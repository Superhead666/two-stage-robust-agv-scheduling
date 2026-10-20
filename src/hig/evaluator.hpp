// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "robust.hpp"

namespace renewed {
using namespace agv;
// Only changes the deterministic scheduling routine called by legacy coverage.
// 0 = original DP, 1 = direct block MIP, 2 = short DP then MIP.
void set_evaluator(GRBEnv* env, int kind);
DetResult block_mip(const Local&, int bound, GRBEnv&, Clock&, bool robust_static=false);
struct EvalStats {
  long calls=0, dp_calls=0, mip_calls=0;
  double dp_seconds=0, mip_seconds=0;
};
EvalStats evaluator_stats();
void reset_evaluator_stats();
}
