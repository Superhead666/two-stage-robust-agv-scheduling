// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"
#include "block_mip.hpp"

namespace baseline {
namespace {
thread_local GRBEnv* environment = nullptr;
thread_local OracleStats stats;
}
void set_oracle(GRBEnv* env) { environment = env; }
void reset_oracle_stats() { stats = {}; }
OracleStats oracle_stats() { return stats; }

DetResult deterministic_mip(const Local& a, int target, GRBEnv& env, Clock& budget) {
  Clock timing;
  ++stats.mip_calls;
  auto out = agv::deterministic_block_mip(a,target,env,budget);
  stats.mip_seconds += timing.elapsed();
  return out;
}
}

namespace agv {
// core.cpp is separately compiled with its DFS renamed deterministic_dp.
// Production robust verification always calls the Gurobi formulation above.
DetResult deterministic(const Local& a, int target, const DetOptions&, Clock& budget) {
  if (!baseline::environment) throw std::runtime_error("deterministic MIP environment not set");
  ++baseline::stats.calls;
  auto r = baseline::deterministic_mip(a, target, *baseline::environment, budget);
  baseline::stats.seconds += r.seconds;
  return r;
}
}
