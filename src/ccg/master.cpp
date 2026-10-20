// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"

namespace baseline {
namespace {
// Model construction counts against the same allowance as optimization.
void check_time(Clock& t) {
  if (t.expired()) throw std::runtime_error("construction_time_limit");
}
}
Master::Master(const Instance& data, GRBEnv& env, const Options& options,
               Clock& budget,
               const std::vector<int>* fixed)
    : a(data), opt(options), model(env), horizon(time_horizon(data)),
      block_limits(master_block_limits(data, options.master_target_upper)) {
  Clock build;
  if(opt.threads!=1)throw std::runtime_error("MILP master requires one solver thread");
  if (opt.formulation != "block" && opt.formulation != "position")
    throw std::runtime_error("formulation must be block or position");
  configure(model, budget.left());
  model.set(GRB_IntParam_Seed, static_cast<int>(opt.seed));
  // Integer time data allow feasibility bisection in integer tardiness.
  tardiness = model.addVar(0, horizon, 0, GRB_INTEGER, "Tmax");
  x.resize(a.n);
  for (int j = 0; j < a.n; ++j) {
    GRBLinExpr sum = 0;
    for (int k = 0; k < a.m; ++k) {
      int lo = fixed && (*fixed)[j] == k ? 1 : 0;
      int hi = fixed && (*fixed)[j] != k ? 0 : 1;
      x[j].push_back(model.addVar(lo, hi, 0, GRB_BINARY));
      if (!fixed && opt.warm_assignment.size() == static_cast<std::size_t>(a.n))
        x[j].back().set(GRB_DoubleAttr_Start, opt.warm_assignment[j] == k ? 1.0 : 0.0);
      sum += x[j].back();
    }
    model.addConstr(sum == 1);
  }
  // Canonical first use removes only label symmetry among machines whose
  // processing, energy, and capacity data are identical for every job.
  if (!fixed) for (const auto& group : identical_machine_groups(a))
    for (std::size_t q = 1; q < group.size(); ++q)
      for (int j = 0; j < a.n; ++j) {
        GRBLinExpr predecessor = 0;
        for (int i = 0; i < j; ++i) predecessor += x[i][group[q-1]];
        model.addConstr(x[j][group[q]] <= predecessor);
      }
  order.resize(a.n); rank.resize(a.n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
    return a.d[i] < a.d[j];
  });
  for (int h = 0; h < a.n; ++h) rank[order[h]] = h + 1;
  check_time(budget);
  model.update();
  build_seconds += build.elapsed();
}

bool Master::add_scenario(const std::vector<double>& point, Clock& budget) {
  if (point.size() != static_cast<std::size_t>(a.n)) throw std::runtime_error("scenario dimension");
  double used = 0;
  for (double v : point) {
    if (!std::isfinite(v) || v < 0 || v > 1) throw std::runtime_error("scenario coordinate");
    used += v;
  }
  if (used > a.gamma + 1e-10) throw std::runtime_error("scenario budget");
  // Do not merge merely similar scenarios. Tiny differences can change packing.
  if (std::find(points.begin(), points.end(), point) != points.end()) return false;
  Clock build;
  if (opt.formulation == "block") add_blocks(point, budget);
  else add_positions(point, budget);
  model.update(); points.push_back(point); build_seconds += build.elapsed();
  return true;
}

void Master::add_blocks(const std::vector<double>& point, Clock& budget) {
  for (int k = 0; k < a.m; ++k) {
    check_time(budget);
    const int slots = block_limits[k];
    std::vector<std::vector<GRBVar>> z(slots);
    std::vector<GRBVar> used, start;
    for (int b = 0; b < slots; ++b) {
      used.push_back(model.addVar(0, 1, 0, GRB_BINARY));
      start.push_back(model.addVar(0, horizon, 0, GRB_CONTINUOUS));
    }
    for (int j = 0; j < a.n; ++j) {
      check_time(budget);
      GRBLinExpr sum = 0;
      for (int b = 0; b < slots; ++b) {
        z[b].push_back(model.addVar(0, 1, 0, GRB_BINARY));
        sum += z[b].back(); model.addConstr(z[b].back() <= used[b]);
      }
      model.addConstr(sum == x[j][k]);
    }
    model.addConstr(start[0] == 0);
    for (int b = 0; b < slots; ++b) {
      check_time(budget);
      GRBLinExpr size = 0, energy = 0, processing = 0;
      for (int j = 0; j < a.n; ++j) {
        size += z[b][j]; processing += a.p[j][k] * z[b][j];
        energy += (a.nominal[j][k] + a.deviation[j][k]*point[j]) / a.capacity[k] * z[b][j];
      }
      model.addConstr(size >= used[b]);
      model.addConstr(energy <= used[b]);
      // Exact no-idle recurrence, including the unused tail of block positions.
      if (b + 1 < slots) {
        model.addConstr(used[b+1] <= used[b]);
        model.addConstr(start[b+1] == start[b] + processing + a.alpha*used[b+1]);
      }
      GRBLinExpr inside = 0;
      for (int j : order) {
        inside += a.p[j][k]*z[b][j];
        const int big_m = std::max(0, horizon-a.d[j]);
        model.addConstr(start[b] + inside <= a.d[j] + tardiness
                        + big_m*(1-z[b][j]));
      }
    }
  }
}

void Master::add_positions(const std::vector<double>& point, Clock& budget) {
  for (int k = 0; k < a.m; ++k) {
    check_time(budget);
    const int slots = a.n;
    std::vector<std::vector<GRBVar>> z(slots);
    std::vector<GRBVar> used, charge, completion, energy;
    for (int q = 0; q < slots; ++q) {
      used.push_back(model.addVar(0, 1, 0, GRB_BINARY));
      charge.push_back(model.addVar(0, q ? 1 : 0, 0, GRB_BINARY));
      completion.push_back(model.addVar(0, horizon, 0, GRB_CONTINUOUS));
      energy.push_back(model.addVar(0, 1, 0, GRB_CONTINUOUS));
    }
    for (int j = 0; j < a.n; ++j) {
      check_time(budget);
      GRBLinExpr sum = 0;
      for (int q = 0; q < slots; ++q) {
        z[q].push_back(model.addVar(0, 1, 0, GRB_BINARY));
        sum += z[q].back();
      }
      model.addConstr(sum == x[j][k]);
    }
    GRBLinExpr previous_rank = 0;
    for (int q = 0; q < slots; ++q) {
      check_time(budget);
      GRBLinExpr size = 0, demand = 0, processing = 0, current_rank = 0;
      for (int j = 0; j < a.n; ++j) {
        size += z[q][j]; processing += a.p[j][k]*z[q][j];
        demand += (a.nominal[j][k] + a.deviation[j][k]*point[j]) / a.capacity[k] * z[q][j];
        current_rank += rank[j]*z[q][j];
        model.addGenConstrIndicator(z[q][j], 1, completion[q] - tardiness,
                                    GRB_LESS_EQUAL, a.d[j]);
      }
      model.addConstr(size == used[q]);
      model.addConstr(charge[q] <= used[q]);
      model.addConstr(energy[q] >= demand);
      if (!q) {
        model.addConstr(completion[q] == processing);
        model.addConstr(energy[q] == demand);
      } else {
        model.addConstr(used[q] <= used[q-1]);
        model.addConstr(completion[q] == completion[q-1] + processing + a.alpha*charge[q]);
        // Energy since the latest recharge. The second lower bound is essential
        // when a recharge is made, otherwise the current job could consume zero.
        model.addConstr(energy[q] >= energy[q-1] + demand - charge[q]);
        model.addConstr(current_rank >= previous_rank + used[q]
                        - (a.n+1)*(charge[q] + 1 - used[q]));
      }
      previous_rank = current_rank;
    }
  }
}

Master::Answer Master::solve(int target, Clock& budget) {
  Answer out;
  if (budget.expired()) { out.reason = "time_limit"; return out; }
  tardiness.set(GRB_DoubleAttr_LB, target);
  tardiness.set(GRB_DoubleAttr_UB, target);
  configure(model, budget.left());
  model.set(GRB_IntParam_Seed, static_cast<int>(opt.seed));
  model.set(GRB_IntParam_SolutionLimit, 1);
  model.optimize();
  out.nodes = static_cast<long>(model.get(GRB_DoubleAttr_NodeCount));
  int status = model.get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE) { out.status = Status::Infeasible; return out; }
  if (model.get(GRB_IntAttr_SolCount) > 0) {
    out.assignment.assign(a.n, -1);
    for (int j = 0; j < a.n; ++j) for (int k = 0; k < a.m; ++k)
      if (x[j][k].get(GRB_DoubleAttr_X) > .5) out.assignment[j] = k;
    if (std::find(out.assignment.begin(), out.assignment.end(), -1) != out.assignment.end())
      throw std::runtime_error("incomplete master assignment");
    out.status = Status::Feasible;
    return out;
  }
  out.reason = status == GRB_TIME_LIMIT ? "time_limit" :
               status == GRB_MEM_LIMIT ? "memory_limit" : "master_status_" + std::to_string(status);
  return out;
}
int Master::variables() { return model.get(GRB_IntAttr_NumVars); }
int Master::constraints() {
  return model.get(GRB_IntAttr_NumConstrs) + model.get(GRB_IntAttr_NumGenConstrs);
}

#ifndef AGV_WITH_CP
std::unique_ptr<ScenarioMaster> make_master(const Instance& a,GRBEnv& env,
    const Options& options,Clock& budget,const std::vector<int>* fixed) {
  if(options.formulation!="block"&&options.formulation!="position")
    throw std::runtime_error("formulation must be block or position");
  return std::make_unique<Master>(a,env,options,budget,fixed);
}
#endif
}
