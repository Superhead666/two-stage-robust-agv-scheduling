// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

namespace baseline {
namespace {
struct InitialSolution {
  int value = INF;
  int lower = 0;
  std::vector<int> assignment;
};

InitialSolution read_warm_solution(const Instance& a, const fs::path& path) {
  std::ifstream in(path);
  int n = 0, upper = INF, local_lower = 0;
  if (!(in >> n >> upper >> local_lower) || n != a.n || upper < 0 || upper >= INF)
    throw std::runtime_error("invalid warm solution header");
  InitialSolution out;
  out.value = upper;
  out.lower = local_lower;
  out.assignment.resize(a.n);
  for (int& k : out.assignment)
    if (!(in >> k) || k < 0 || k >= a.m)
      throw std::runtime_error("invalid warm solution assignment");
  return out;
}

InitialSolution singleton_edd_list(const Instance& a) {
  std::vector<int> order(a.n), ready(a.m, 0), count(a.m, 0);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
    return a.d[i] < a.d[j];
  });
  InitialSolution out; out.value = 0; out.assignment.assign(a.n, -1);
  for (int j : order) {
    int selected = -1, completion = INF;
    for (int k = 0; k < a.m; ++k) {
      const int candidate = ready[k] + (count[k] ? a.alpha : 0) + a.p[j][k];
      if (candidate < completion) { completion = candidate; selected = k; }
    }
    if (selected < 0) throw std::runtime_error("initial assignment failed");
    out.assignment[j] = selected; ready[selected] = completion; ++count[selected];
    out.value = std::max(out.value, completion-a.d[j]);
  }
  out.value = std::max(0, out.value);
  return out;
}

int next_target(int lower, int upper, const std::set<int>& unresolved) {
  if (lower >= upper) return -1;
  const int middle=lower+(upper-lower-1)/2;
  for(int distance=0;distance<upper-lower;++distance) {
    const int above=middle+distance;
    if(above<upper&&!unresolved.count(above))return above;
    const int below=middle-distance;
    if(distance&&below>=lower&&!unresolved.count(below))return below;
  }
  return -1;
}

int deadline_prefix_bound(const Local& a) {
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  int best=-INF,time=0;Mask mask=0;
  for(int h=0;h<a.n();++h) {
    const int j=order[h];time+=a.p[j];mask|=Mask(1)<<j;
    if(h+1<a.n()&&a.d[order[h+1]]==a.d[j])continue;
    const int recharges=std::max(0,int(std::ceil((worst_block_energy(a,mask)-1e-9)/a.capacity))-1);
    best=std::max(best,time+a.alpha*recharges-a.d[j]);
  }
  return best;
}
}

std::vector<std::vector<int>> identical_machine_groups(const Instance& a) {
  std::vector<std::vector<int>> groups;
  std::vector<bool> used(a.m, false);
  auto same = [&](int k, int h) {
    if (a.capacity[k] != a.capacity[h]) return false;
    for (int j = 0; j < a.n; ++j)
      if (a.p[j][k] != a.p[j][h] || a.nominal[j][k] != a.nominal[j][h]
          || a.deviation[j][k] != a.deviation[j][h]) return false;
    return true;
  };
  for (int k = 0; k < a.m; ++k) if (!used[k]) {
    std::vector<int> group{k}; used[k] = true;
    for (int h = k+1; h < a.m; ++h) if (!used[h] && same(k,h)) {
      group.push_back(h); used[h] = true;
    }
    if (group.size() > 1) groups.push_back(std::move(group));
  }
  return groups;
}

std::vector<int> master_block_limits(const Instance& a, int exclusive_target_upper) {
  std::vector<int> limits(a.m, a.n);
  if (exclusive_target_upper >= INF/2) return limits;
  if (exclusive_target_upper <= 0) return std::vector<int>(a.m, 1);
  const long long completion_limit = static_cast<long long>(
      *std::max_element(a.d.begin(), a.d.end())) + exclusive_target_upper - 1;
  for (int k = 0; k < a.m; ++k) {
    std::vector<int> processing;
    processing.reserve(a.n);
    for (int j = 0; j < a.n; ++j) processing.push_back(a.p[j][k]);
    std::sort(processing.begin(), processing.end());
    long long total = 0;
    int count = 0;
    for (int q = 1; q <= a.n; ++q) {
      total += processing[q-1];
      if (total + static_cast<long long>(a.alpha)*(q-1) > completion_limit) break;
      count = q;
    }
    limits[k] = std::max(1, count);
  }
  return limits;
}

int time_horizon(const Instance& a) {
  // Every job fits a full battery. At most n-1 recharges suffice for any
  // assignment, with no idle time. This is a data bound, not an incumbent.
  long long horizon = static_cast<long long>(std::max(0, a.n-1))*a.alpha;
  for (const auto& row : a.p) horizon += *std::max_element(row.begin(), row.end());
  if (horizon >= INF/8) throw std::runtime_error("unsafe time horizon");
  return static_cast<int>(horizon);
}
RobustOptions baseline_verifier_options() {
  RobustOptions ro;
  // The baseline and LBBD use the same exact fixed-assignment verifier.
  ro.enhanced = true; ro.shared_blocks = true; ro.early_adversary = true;
  ro.overload_cuts = false; ro.static_policy = false; ro.branch_coverage = false;
  ro.fixed_order = false; ro.order_generation = false; ro.order_bounds = false;
  ro.workload_screen = true; ro.refine = true; ro.reuse_schedules = true;
  ro.screening = false;
  ro.max_schedules = std::numeric_limits<int>::max();
  return ro;
}

void write_solution(const fs::path& path, const Result& r) {
  if (path.empty()) return;
  // Only this run's checkpoint is replaced; existing experiment runs are guarded by the driver.
  fs::path temporary = path; temporary += ".tmp";
  std::ofstream f(temporary);
  if (!f) throw std::runtime_error("cannot write solution checkpoint");
  f << r.assignment.size() << ' ' << r.upper << ' ' << r.lower << '\n';
  for (int k : r.assignment) f << k << ' ';
  f << '\n' << "objective maximum_tardiness\ncertificate " << r.certificate << '\n';
  f.close(); if (!f) throw std::runtime_error("solution checkpoint write failed");
  fs::rename(temporary, path);
}

Result solve(const Instance& a, GRBEnv& env, const Options& opt) {
  if (!(opt.limit > 0) || !std::isfinite(opt.limit)) throw std::runtime_error("positive time limit required");
  if (opt.threads != 1) throw std::runtime_error("fair-comparison baselines require one solver thread");
  env.set(GRB_IntParam_Threads,opt.threads);
  env.set(GRB_IntParam_Seed,static_cast<int>(opt.seed));
  validate_instance(a);
  for (int d : a.d) if (d < 0) throw std::runtime_error("baseline requires nonnegative due dates");
  long long magnitude = static_cast<long long>(a.n)*a.alpha;
  for (int j = 0; j < a.n; ++j) magnitude += *std::max_element(a.p[j].begin(), a.p[j].end());
  if (magnitude + *std::max_element(a.d.begin(), a.d.end()) > INF/8)
    throw std::runtime_error("instance time magnitudes exceed safe range");
  Clock budget(opt.limit); Result out;
  const int horizon = time_horizon(a);
  const auto initial = opt.warm_solution.empty()
      ? singleton_edd_list(a) : read_warm_solution(a, opt.warm_solution);
  if(!opt.warm_solution.empty()&&initial.lower!=opt.initial_lower)
    throw std::runtime_error("warm solution and initial lower bound disagree");
  out.initial = out.upper = initial.value;
  if (opt.initial_lower < 0 || opt.initial_lower > out.upper)
    throw std::runtime_error("invalid initial lower bound");
  out.lower = opt.initial_lower;
  out.assignment = initial.assignment;
  out.certificate = opt.warm_solution.empty() ? "singleton_edd_list" : "certified_hig_warm_start";
  out.time_best = 0;
  std::ofstream trace;
  if (!opt.trace.empty()) {
    trace.open(opt.trace);
    if (!trace) throw std::runtime_error("cannot write trace");
    trace << "elapsed,event,iteration,target,lower,upper,scenarios\n";
  }
  auto record = [&](const std::string& event, int target) {
    if (trace.is_open()) {
      trace << std::setprecision(10) << budget.elapsed() << ',' << event << ',' << out.iterations
            << ',' << target << ',' << out.lower << ',';
      if (out.upper < INF) trace << out.upper;
      trace << ',' << out.scenarios << '\n';
      trace.flush();
    }
    write_solution(opt.solution, out);
  };
  record("initial", out.upper);
  if (opt.formulation != "block" && opt.formulation != "position" && opt.formulation != "cp")
    throw std::runtime_error("formulation must be block, position, or cp");
  std::unique_ptr<ScenarioMaster> master;
  try {
    if (out.lower < out.upper && !budget.expired()) {
      Options master_options = opt;
      master_options.master_target_upper = out.upper;
      master_options.warm_assignment = initial.assignment;
      master = make_master(a, env, master_options, budget);
      master->add_scenario(std::vector<double>(a.n, 0), budget);
      out.scenarios = master->scenarios();
    }
    const auto base_ro = baseline_verifier_options();
    std::vector<std::unordered_map<Mask,std::vector<Schedule>>> pools(a.m);
    std::set<int> unresolved_targets;
    while (master && out.lower < out.upper && !budget.expired()) {
      const int target = next_target(out.lower, out.upper, unresolved_targets);
      if (target < 0) { out.reason = "unresolved_target_set"; break; }
      ++out.iterations; record("master_start", target);
      Clock mt; auto candidate = master->solve(target, budget);
      out.master_seconds += mt.elapsed(); out.master_nodes += candidate.nodes;
      if (candidate.status == Status::Infeasible) {
        out.lower = target+1; out.reason.clear(); record("master_infeasible", target); continue;
      }
      if (candidate.status == Status::Unknown) {
        out.reason = candidate.reason;
        record("master_unresolved", target);
        if (budget.expired()) break;
        unresolved_targets.insert(target); continue;
      }
      bool feasible = true, unresolved = false;
      std::vector<double> adverse;
      std::vector<std::pair<int,int>> ranked;
      for(int k=0;k<a.m;++k) {
        auto local=a.local(k,candidate.assignment);
        if(local.n())ranked.push_back({deadline_prefix_bound(local),k});
      }
      std::stable_sort(ranked.begin(),ranked.end(),[](const auto& x,const auto& y){
        if(x.first!=y.first)return x.first>y.first;return x.second<y.second;
      });
      for (const auto& item:ranked) {
        const int k=item.second;
        Local local = a.local(k, candidate.assignment);
        if (!local.n()) continue;
        if (budget.expired()) { feasible = false; unresolved = true; out.reason = "time_limit"; break; }
        Mask mask=0;for(int j:local.ids)mask|=Mask(1)<<j;
        auto ro=base_ro;ro.screening=local.gamma>=local.n();
        Clock vt; auto r = verify(local, target, env, ro, budget,&pools[k][mask]);
        out.verify_seconds += vt.elapsed(); ++out.verify_calls;
        out.adversary_calls += r.adversary_calls; out.schedules += r.schedules;
        if (r.status == Status::Unknown) {
          feasible = false; unresolved = true; out.reason = r.reason; break;
        }
        if (r.status == Status::Infeasible) {
          feasible = false; adverse.assign(a.n, 0);
          if (r.witness.size() != static_cast<std::size_t>(local.n())) throw std::runtime_error("missing adverse witness");
          for (int q = 0; q < local.n(); ++q) adverse[local.ids[q]] = r.witness[q];
          break;
        }
      }
      if (feasible) {
        if (target < out.upper) {
          out.upper = target; out.assignment = candidate.assignment;
          out.certificate = "adaptive_coverage_solver_tolerance"; out.time_best = budget.elapsed();
        }
        out.reason.clear(); record("robust_feasible", target);
      } else if (unresolved || budget.expired()) {
        record("verification_unresolved", target);
        if (budget.expired()) break;
        unresolved_targets.insert(target); continue;
      } else {
        if (!master->add_scenario(adverse, budget)) {
          out.reason = "duplicate_scenario_numerical_boundary";
          record("duplicate_scenario", target);
          if (budget.expired()) break;
          unresolved_targets.insert(target); continue;
        }
        unresolved_targets.clear(); out.reason.clear();
        out.scenarios = master->scenarios(); record("scenario_added", target);
      }
    }
  } catch (const std::runtime_error& error) {
    if (std::string(error.what()) != "construction_time_limit") throw;
    out.reason = "construction_time_limit";
  }
  if (master) {
    out.scenarios = master->scenarios(); out.build_seconds = master->build_seconds;
    out.extraction_seconds = master->extraction_time();
    out.master_variables = master->variables(); out.master_constraints = master->constraints();
  }
  if (out.lower > out.upper) throw std::runtime_error("global bound contradiction");
  if ((out.upper < INF) != (out.assignment.size() == static_cast<std::size_t>(a.n)))
    throw std::runtime_error("incumbent/certificate mismatch");
  if (out.lower == out.upper) { out.status = "optimal_solver_tolerance"; out.reason = "bounds_match"; }
  else if (budget.expired() || out.reason == "time_limit") { out.status = "time_limit"; if (out.reason.empty()) out.reason = "time_limit"; }
  else if (out.reason == "memory_limit") out.status = "memory_limit";
  else { out.status = "unresolved"; if (out.reason.empty()) out.reason = "verification_unresolved"; }
  out.seconds = budget.elapsed(); record("finished", out.upper < INF ? out.upper : horizon);
  return out;
}

Instance generate_paper(int n, int m, unsigned seed, double gamma,
                        double variation, double scale, double due_factor,
                        int processing_max, int energy_max, double recharge_factor, int recharge_time) {
  if (n < 1 || n > 62 || m < 1 || m > 32 || gamma < 0 || gamma > n ||
      variation < 0 || scale <= 0 || due_factor < 0 || processing_max < 1 ||
      processing_max > 10000 || energy_max < 1 || energy_max > 10000 ||
      !std::isfinite(recharge_factor) || recharge_factor < 0 || recharge_factor > 100 ||
      recharge_time < -1 || recharge_time > 1000000)
    throw std::runtime_error("invalid generator parameters");
  Instance a; a.n=n; a.m=m; a.gamma=gamma; a.id="s"+std::to_string(seed);
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> processing_uniform(1,processing_max),energy_uniform(1,energy_max);
  std::uniform_real_distribution<double> factor(.8,1.2), due(.25,1.1);
  a.p.assign(n,std::vector<int>(m)); a.nominal.assign(n,std::vector<double>(m)); a.deviation=a.nominal;
  a.capacity.assign(m,0); a.d.resize(n); double processing=0;
  for (int j=0;j<n;++j) for (int k=0;k<m;++k) {
    a.p[j][k]=processing_uniform(rng); processing+=a.p[j][k]; a.nominal[j][k]=energy_uniform(rng);
    a.deviation[j][k]=a.nominal[j][k]*variation*factor(rng);
  }
  a.alpha=recharge_time>=0?recharge_time:static_cast<int>(std::round(recharge_factor*processing/(n*m)));
  for (int k=0;k<m;++k) {
    double sum=0,maximum=0;
    for (int j=0;j<n;++j) { double e=a.nominal[j][k]+a.deviation[j][k]; sum+=e; maximum=std::max(maximum,e); }
    a.capacity[k]=std::max(maximum,scale*sum/n);
  }
  double horizon=(processing/(n*m)+a.alpha/4.0)*n/m;
  for (int j=0;j<n;++j) a.d[j]=static_cast<int>(std::round(due_factor*horizon*due(rng)));
  validate_instance(a); return a;
}
}
