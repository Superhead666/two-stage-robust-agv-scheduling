// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "block_mip.hpp"

namespace agv {

DetResult deterministic_block_mip(const Local& a, int target, GRBEnv& env,
                                  Clock& budget, bool robust_static) {
  Clock timing;
  DetResult out;
  auto finish = [&](Status status) {
    out.status = status;
    out.seconds = timing.elapsed();
    return out;
  };
  if (!a.n()) return finish(Status::Feasible);
  if (a.n() > 62 || budget.expired()) return finish(Status::Unknown);

  Local upper = a;
  if (robust_static)
    for (int j = 0; j < a.n(); ++j)
      upper.e[j] = a.nominal[j] + std::min(1.0, a.gamma)*a.deviation[j];
  auto heuristic = edd(upper);
  if (heuristic.lateness <= target) {
    out.schedule = heuristic;
    return finish(Status::Feasible);
  }

  const int total = std::accumulate(a.p.begin(), a.p.end(), 0);
  const int horizon = *std::max_element(a.d.begin(), a.d.end()) + target;
  if (total > horizon) return finish(Status::Infeasible);
  for (int j = 0; j < a.n(); ++j)
    if (a.p[j] > a.d[j] + target || upper.e[j] > a.capacity + ENERGY_EPS)
      return finish(Status::Infeasible);

  const int count = a.alpha
      ? std::min(a.n(), 1 + (horizon-total)/a.alpha)
      : a.n();
  if (count < 1) return finish(Status::Infeasible);

  std::vector<int> order(a.n());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int i, int j) { return a.d[i] < a.d[j]; });

  GRBModel model(env);
  configure(model, budget.left());
  model.set(GRB_IntParam_SolutionLimit, 1);
  std::vector<std::vector<GRBVar>> x(a.n());
  std::vector<GRBVar> used;
  for (int b = 0; b < count; ++b)
    used.push_back(model.addVar(b == 0 ? 1 : 0, 1, 0, GRB_BINARY));

  for (int j = 0; j < a.n(); ++j) {
    GRBLinExpr assignment = 0;
    for (int b = 0; b < count; ++b) {
      const bool possible = a.p[j] + a.alpha*b <= a.d[j] + target;
      x[j].push_back(model.addVar(0, possible ? 1 : 0, 0, GRB_BINARY));
      assignment += x[j].back();
      model.addConstr(x[j].back() <= used[b]);
    }
    model.addConstr(assignment == 1);
  }

  GRBLinExpr preceding = 0;
  for (int b = 0; b < count; ++b) {
    GRBLinExpr energy = 0, jobs = 0, work = 0;
    for (int j = 0; j < a.n(); ++j) {
      energy += (robust_static ? a.nominal[j] : a.e[j])/a.capacity*x[j][b];
      jobs += x[j][b];
      work += a.p[j]*x[j][b];
    }
    if (robust_static && a.gamma > 0) {
      auto price = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
      energy += a.gamma*price;
      for (int j = 0; j < a.n(); ++j) {
        auto excess = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model.addConstr(excess >= a.deviation[j]/a.capacity*x[j][b] - price);
        energy += excess;
      }
    }
    model.addConstr(energy <= used[b]);
    model.addConstr(jobs >= used[b]);
    if (b) model.addConstr(used[b] <= used[b-1]);
    GRBLinExpr partial = 0;
    for (int j : order) {
      partial += a.p[j]*x[j][b];
      const int big_m = std::max(0, total + a.alpha*b - a.d[j] - target);
      model.addConstr(preceding + partial + a.alpha*b <= a.d[j] + target
                      + big_m*(1-x[j][b]));
    }
    preceding += work;
  }

  if (budget.expired()) return finish(Status::Unknown);
  configure(model, budget.left());
  model.optimize();
  out.nodes = static_cast<std::uint64_t>(model.get(GRB_DoubleAttr_NodeCount));
  if (model.get(GRB_IntAttr_SolCount) > 0) {
    std::vector<int> sequence;
    std::vector<bool> cuts;
    for (int b = 0; b < count; ++b) {
      bool first = true;
      for (int j : order)
        if (x[j][b].get(GRB_DoubleAttr_X) > .5) {
          cuts.push_back(first && !sequence.empty());
          sequence.push_back(j);
          first = false;
        }
    }
    out.schedule = evaluate(a, sequence, cuts);
    if (out.schedule.lateness > target) return finish(Status::Unknown);
    if (robust_static)
      for (Mask block : out.schedule.blocks)
        if (worst_block_energy(a, block) > a.capacity + ENERGY_EPS)
          return finish(Status::Unknown);
    return finish(Status::Feasible);
  }
  return finish(model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE
                    ? Status::Infeasible : Status::Unknown);
}

}
