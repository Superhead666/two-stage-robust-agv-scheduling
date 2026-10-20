// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once

#include "core.hpp"
#include "gurobi_c++.h"
#include <cmath>
#include <unordered_set>

namespace lb {

struct PrefixDualResult {
  bool complete=false;
  double value=0;
  double seconds=0;
  int columns=0;
  std::vector<double> weight;
};

// Solve the fractional bin-packing relaxation for one fixed, admissible
// energy realization.  Column generation stops only after an optimal
// knapsack pricing solve proves that the returned dual weights satisfy every
// feasible one-block pattern.  Scaling by the final pricing value protects
// the cut against the pricing solver's feasibility tolerance.
inline PrefixDualResult prefix_dual(const agv::Local& a,GRBEnv& env,agv::Clock& global,double limit) {
  PrefixDualResult out;out.weight.assign(a.n(),0.0);agv::Clock timing;
  if(a.n()==0||global.expired()||limit<=0){out.seconds=timing.elapsed();return out;}
  const double allowed=std::min(limit,global.left());
  if(allowed<=0){out.seconds=timing.elapsed();return out;}
  for(double energy:a.e)if(energy>a.capacity+1e-8||energy<0) {
    out.seconds=timing.elapsed();return out;
  }
  try {
    GRBModel master(env);master.set(GRB_IntParam_OutputFlag,0);master.set(GRB_IntParam_Threads,1);
    master.set(GRB_IntParam_Method,1);master.set(GRB_DoubleParam_TimeLimit,std::max(.001,allowed));
    std::vector<GRBConstr> cover;
    for(int j=0;j<a.n();++j)cover.push_back(master.addConstr(GRBLinExpr()>=1.0));
    std::unordered_set<agv::Mask> patterns;
    auto add_pattern=[&](agv::Mask mask) {
      GRBColumn column;
      for(int j=0;j<a.n();++j)if(mask&(agv::Mask(1)<<j))column.addTerm(1.0,cover[j]);
      master.addVar(0,GRB_INFINITY,1.0,GRB_CONTINUOUS,column);
      patterns.insert(mask);++out.columns;
    };
    for(int j=0;j<a.n();++j)add_pattern(agv::Mask(1)<<j);

    std::vector<double> dual(a.n());double pricing_value=GRB_INFINITY;
    while(!global.expired()&&timing.elapsed()<allowed) {
      master.set(GRB_DoubleParam_TimeLimit,std::max(.001,std::min(global.left(),allowed-timing.elapsed())));
      master.optimize();
      if(master.get(GRB_IntAttr_Status)!=GRB_OPTIMAL)break;
      for(int j=0;j<a.n();++j)dual[j]=std::max(0.0,cover[j].get(GRB_DoubleAttr_Pi));

      GRBModel pricing(env);pricing.set(GRB_IntParam_OutputFlag,0);pricing.set(GRB_IntParam_Threads,1);
      pricing.set(GRB_DoubleParam_MIPGap,0);pricing.set(GRB_DoubleParam_FeasibilityTol,1e-9);
      pricing.set(GRB_DoubleParam_TimeLimit,std::max(.001,std::min(global.left(),allowed-timing.elapsed())));
      std::vector<GRBVar> z;GRBLinExpr energy=0,profit=0;
      for(int j=0;j<a.n();++j) {
        z.push_back(pricing.addVar(0,1,0,GRB_BINARY));energy+=a.e[j]*z.back();profit+=dual[j]*z.back();
      }
      pricing.addConstr(energy<=a.capacity);pricing.setObjective(profit,GRB_MAXIMIZE);pricing.optimize();
      if(pricing.get(GRB_IntAttr_Status)!=GRB_OPTIMAL)break;
      pricing_value=pricing.get(GRB_DoubleAttr_ObjVal);
      agv::Mask mask=0;for(int j=0;j<a.n();++j)if(z[j].get(GRB_DoubleAttr_X)>.5)mask|=agv::Mask(1)<<j;
      if(pricing_value<=1.0+1e-8) {
        const double scale=std::max(1.0,pricing_value);
        for(int j=0;j<a.n();++j)out.weight[j]=dual[j]/scale;
        out.value=0;for(double value:out.weight)out.value+=value;
        out.complete=true;break;
      }
      if(mask==0||patterns.count(mask))break;
      add_pattern(mask);
    }
  } catch(const GRBException&) {
    out.complete=false;
  }
  out.seconds=timing.elapsed();return out;
}

}
