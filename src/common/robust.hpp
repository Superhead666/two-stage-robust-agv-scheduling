// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "core.hpp"
#include "gurobi_c++.h"
#include <memory>

namespace agv {
struct RobustOptions {
  bool enhanced=true, screening=true, shared_blocks=true, early_adversary=true;
  bool workload_screen=false, refine=false;
  bool overload_cuts=false, reuse_schedules=false;
  bool static_policy=false;
  bool branch_coverage=false;
  bool fixed_order=false;
  bool order_generation=false;
  bool order_bounds=false;
  DetOptions det;
  int max_schedules=std::numeric_limits<int>::max();
};
struct RobustResult {
  Status status=Status::Unknown; // Feasible = bound is robustly feasible, within solver tolerances.
  std::vector<double> witness;
  int schedules=0, adversary_calls=0, det_calls=0, blocks=0;
  std::uint64_t det_nodes=0;
  double seconds=0, mip_seconds=0, det_seconds=0, margin=0;
  std::string reason;
};
// All verifier submodels inherit one reproducible solver seed. The default
// preserves legacy behaviour; publication drivers set it explicitly.
void set_solver_seed(int seed);
void configure(GRBModel&,double seconds);
RobustResult verify(const Local&,int lambda,GRBEnv&,const RobustOptions&,Clock&,std::vector<Schedule>* pool=nullptr);
RobustResult verify_full_partitions(const Local&,int lambda,GRBEnv&,Clock&);
RobustResult verify_fixed_order(const Local&,const std::vector<int>& order,int lambda,GRBEnv&,Clock&,bool tight=false);
RobustResult verify_order_enumeration(const Local&,const std::vector<int>& order,int lambda,GRBEnv&,Clock&);
struct FleetOptions {
  bool enhanced=true;
  bool robust_workload=false;
  bool primal_search=false, threshold_master=false, conflict_cuts=false, logic_only=false;
  bool insertion_cuts=false;
  bool root_bound=false;
  bool cumulative_workload=false;
  RobustOptions robust;
  int max_iterations=200;
};
struct FleetResult {
  std::string status="unknown";
  int lower=-INF,upper=INF,iterations=0,scenarios=0;
  double seconds=0,master_seconds=0,verify_seconds=0;
  int verify_calls=0,adversary_calls=0,schedules=0,conflicts=0;
  std::vector<int> assignment;
};
FleetResult solve_fleet(const Instance&,GRBEnv&,const FleetOptions&,Clock&);
int enumerate_fleet(const Instance&,GRBEnv&,Clock&);
}
