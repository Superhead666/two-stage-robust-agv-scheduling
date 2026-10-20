// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "evaluator.hpp"
#include "block_mip.hpp"

namespace agv {
DetResult deterministic_dp(const Local&, int, const DetOptions&, Clock&);
}
namespace renewed {
namespace {
thread_local GRBEnv* selected_env=nullptr;
thread_local int selected_kind=0;
thread_local EvalStats stats;
}
void set_evaluator(GRBEnv* env,int kind) { selected_env=env; selected_kind=kind; }
EvalStats evaluator_stats() { return stats; }
void reset_evaluator_stats() { stats={}; }

DetResult block_mip(const Local& a,int bound,GRBEnv& env,Clock& budget,bool robust_static) {
  return agv::deterministic_block_mip(a,bound,env,budget,robust_static);
}
}

namespace agv {
DetResult deterministic(const Local& a,int bound,const DetOptions& options,Clock& budget) {
  using namespace renewed;
  ++stats.calls;
  if(selected_kind==0 || !selected_env || options.dense) {
    ++stats.dp_calls;auto r=deterministic_dp(a,bound,options,budget);stats.dp_seconds+=r.seconds;return r;
  }
  if(selected_kind==2) {
    Clock short_run(std::min(.003,budget.left()));++stats.dp_calls;
    auto r=deterministic_dp(a,bound,options,short_run);stats.dp_seconds+=r.seconds;
    if(r.status!=Status::Unknown || budget.expired())return r;
  }
  ++stats.mip_calls;auto r=block_mip(a,bound,*selected_env,budget);stats.mip_seconds+=r.seconds;return r;
}
}
