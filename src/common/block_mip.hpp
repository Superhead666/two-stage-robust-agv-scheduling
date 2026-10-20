// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once

#include "robust.hpp"

namespace agv {

// Shared deterministic block formulation used by HIG, LBBD, and MILP--CCG.
DetResult deterministic_block_mip(const Local&, int target, GRBEnv&, Clock&,
                                  bool robust_static = false);

}
