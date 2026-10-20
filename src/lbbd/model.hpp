// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "core.hpp"
#include "gurobi_c++.h"
#include <set>
namespace lb {
using namespace agv;
#ifndef AGV_DUE_DATE_BOUNDS
#define AGV_DUE_DATE_BOUNDS 1
#endif
inline void configure(GRBModel& model,double time,int focus=0,int seed=1) {
  model.set(GRB_IntParam_OutputFlag,0);model.set(GRB_IntParam_Threads,1);
  model.set(GRB_IntParam_Seed,seed);model.set(GRB_IntParam_MIPFocus,focus);
  model.set(GRB_DoubleParam_TimeLimit,std::max(.001,time));
  model.set(GRB_DoubleParam_MIPGap,0);model.set(GRB_DoubleParam_FeasibilityTol,1e-9);
}
struct Relaxation {
  const Instance& a;GRBModel model;GRBVar lambda;
  std::vector<std::vector<GRBVar>> x,r;
  std::vector<int> dates;std::vector<std::vector<int>> prefixes;
  Relaxation(const Instance& data,GRBEnv& env,bool nonnegative=true):a(data),model(env),x(a.n),r(a.m) {
    int low=0;
    if(!nonnegative){low=-INF;for(int j=0;j<a.n;++j)low=std::max(low,*std::min_element(a.p[j].begin(),a.p[j].end())-a.d[j]);}
    lambda=model.addVar(low,GRB_INFINITY,1,GRB_INTEGER,"lambda");
    for(int j=0;j<a.n;++j) {
      GRBLinExpr sum=0;
      for(int k=0;k<a.m;++k){x[j].push_back(model.addVar(0,1,0,GRB_BINARY));sum+=x[j].back();}
      model.addConstr(sum==1);
    }
    std::set<int> ds(a.d.begin(),a.d.end());dates.assign(ds.begin(),ds.end());
    for(int d:dates){std::vector<int> p;for(int j=0;j<a.n;++j)if(a.d[j]<=d)p.push_back(j);prefixes.push_back(p);}
#if AGV_DUE_DATE_BOUNDS
    for(int k=0;k<a.m;++k)for(int h=0;h<int(dates.size());++h) {
      r[k].push_back(model.addVar(0,a.n-1,0,GRB_INTEGER));
      auto price=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);
      GRBLinExpr energy=a.gamma*price,processing=0;
      for(int j:prefixes[h]) {
        energy+=a.nominal[j][k]/a.capacity[k]*x[j][k];processing+=a.p[j][k]*x[j][k];
        auto excess=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);energy+=excess;
        model.addConstr(excess>=a.deviation[j][k]/a.capacity[k]*x[j][k]-price);
      }
      model.addConstr(energy<=1+r[k].back());model.addConstr(processing+a.alpha*r[k].back()<=dates[h]+lambda);
    }
#else
    // Processing-only prefix relaxation used by the component ablation.  The
    // recharge counts remain available to the DFF separator, but the
    // budgeted-volume inequalities are omitted.
    for(int k=0;k<a.m;++k)for(int h=0;h<int(dates.size());++h) {
      r[k].push_back(model.addVar(0,a.n-1,0,GRB_INTEGER));GRBLinExpr processing=0;
      for(int j:prefixes[h])processing+=a.p[j][k]*x[j][k];
      model.addConstr(processing+a.alpha*r[k].back()<=dates[h]+lambda);
    }
#endif
  }
};
inline int integer_bound(double z) {
  if(!std::isfinite(z)||z<=0)return 0;
  if(z>=INF)return INF;
  return static_cast<int>(std::ceil(z-1e-6));
}
}
