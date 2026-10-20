// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "core.hpp"
namespace lb {
inline double budget_sum(std::vector<double> sizes,double gamma) {
  std::sort(sizes.rbegin(),sizes.rend());double sum=0;
  for(double v:sizes){double use=std::min(1.0,gamma);sum+=use*v;gamma-=use;if(gamma<=0)break;}return sum;
}
inline std::vector<double> polymatroid_coefficients(const std::vector<double>& deviation,double gamma,const std::vector<int>& order) {
  std::vector<double> coefficient(deviation.size(),0),selected;double previous=0;
  for(int j:order){selected.push_back(deviation[j]);double next=budget_sum(selected,gamma);coefficient[j]=std::max(0.,next-previous);previous=next;}return coefficient;
}
}
