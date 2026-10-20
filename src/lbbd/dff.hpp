// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include <algorithm>
#include <cmath>
// Conservative rounding DFF. If sum t_i <= 1, then sum f(t_i) <= 1:
// q*sum f(t_i) is an integer strictly smaller than q+1 whenever it is positive.
// Downward adjustment at breakpoints only decreases the valid function.
inline double transformed_size(double t,int q){return std::max(0.0,std::ceil((q+1)*t-1e-9)-1)/q;}
