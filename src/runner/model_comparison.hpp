// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "long_search.hpp"
#include "evaluator.hpp"
#include "model.hpp"
#include <fstream>

// Model comparison only. The submitted LBBD and HIG executables are unchanged.
namespace model_comparison {
using namespace agv;
struct Solution {
  int lower=0,upper=INF;
  std::vector<int> assignment;
  std::vector<Schedule> schedules;
  std::string status="time_limit";
  double seconds=0,nodes=0;
  long verification_calls=0,deterministic_calls=0,adversary_calls=0;
};
inline Solution read_solution(const Instance& a,const std::string& path) {
  std::ifstream in(path);int n=0;Solution s;
  if(!(in>>n>>s.upper>>s.lower)||n!=a.n||s.lower<0||s.lower>s.upper||s.upper>=INF)
    throw std::runtime_error("invalid reference solution header: "+path);
  s.assignment.resize(n);
  for(int& k:s.assignment)if(!(in>>k)||k<0||k>=a.m)
    throw std::runtime_error("invalid reference assignment: "+path);
  return s;
}
inline Schedule static_schedule(const Local& a) {
  // Processing times are deterministic. Zero energies recover timing and masks;
  // battery feasibility is checked with the budgeted block energies below.
  Local zero=a;std::fill(zero.e.begin(),zero.e.end(),0);
  Schedule best;best.lateness=INF;
  for(const auto& blocks:{longrun::static_reference(a),longrun::insertion_reference(a)}) {
    std::vector<int> order;std::vector<bool> cuts;
    for(const auto& block:blocks)for(std::size_t q=0;q<block.size();++q) {
      cuts.push_back(!order.empty()&&q==0);order.push_back(block[q]);
    }
    auto candidate=evaluate(zero,order,cuts);if(candidate.lateness<best.lateness)best=std::move(candidate);
  }
  return best;
}
inline void validate_static(const Instance& a,const Solution& s) {
  if(s.assignment.size()!=std::size_t(a.n)||s.schedules.size()!=std::size_t(a.m))
    throw std::runtime_error("static certificate dimensions");
  for(int k:s.assignment)if(k<0||k>=a.m)throw std::runtime_error("invalid static assignment");
  int value=0;
  for(int k=0;k<a.m;++k) {
    Local local=a.local(k,s.assignment);const auto& schedule=s.schedules[k];
    Local zero=local;std::fill(zero.e.begin(),zero.e.end(),0);
    if(!valid_schedule(zero,schedule,s.upper))throw std::runtime_error("static timing certificate");
    for(Mask block:schedule.blocks)
      if(worst_block_energy(local,block)>local.capacity+ENERGY_EPS)
        throw std::runtime_error("static battery certificate");
    value=std::max(value,schedule.lateness);
  }
  if(value!=s.upper)throw std::runtime_error("static objective certificate");
}
inline Solution static_incumbent(const Instance& a,const std::vector<int>& assignment) {
  Solution s;s.assignment=assignment;s.upper=0;
  for(int k=0;k<a.m;++k) {
    s.schedules.push_back(static_schedule(a.local(k,assignment)));
    s.upper=std::max(s.upper,s.schedules.back().lateness);
  }
  validate_static(a,s);return s;
}

// Exact single-stage robust MILP: assignment, job order, and recharge positions
// are all fixed before energy is observed. Each block must satisfy its own
// budgeted robust counterpart under the SAME global uncertainty set. This does
// not require allocating separate uncertainty budgets to simultaneous blocks.
inline Solution solve_static(const Instance& a,GRBEnv& env,Clock& clock,int seed,
                             const std::vector<std::vector<int>>& seeds,int robust_lower) {
  Solution out;
  for(const auto& x:seeds) {
    auto s=static_incumbent(a,x);if(s.upper<out.upper)out=std::move(s);
  }
  if(out.upper==INF)throw std::runtime_error("static solve needs an initial assignment");
  out.lower=std::max(0,robust_lower);
  if(out.lower>out.upper)throw std::runtime_error("reference bound contradicts static certificate");
  if(out.lower==out.upper){out.status="optimal_certified";out.seconds=clock.elapsed();return out;}
  if(clock.expired()){out.seconds=clock.elapsed();return out;}
  // Use the manuscript's aggregate due-date rows as valid strengthening.
  lb::Relaxation f(a,env);f.lambda.set(GRB_DoubleAttr_LB,out.lower);
  f.lambda.set(GRB_DoubleAttr_UB,out.upper);
  f.lambda.set(GRB_DoubleAttr_Start,out.upper);
  for(int j=0;j<a.n;++j)for(int k=0;k<a.m;++k)
    f.x[j][k].set(GRB_DoubleAttr_Start,out.assignment[j]==k);
  long long horizon=static_cast<long long>(a.n-1)*a.alpha;
  for(const auto& p:a.p)horizon+=*std::max_element(p.begin(),p.end());
  if(horizon>=INF/8)throw std::runtime_error("unsafe static horizon");
  std::vector<int> order(a.n);std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  using Matrix=std::vector<std::vector<GRBVar>>;
  std::vector<Matrix> z(a.m);
  std::vector<std::vector<GRBVar>> used(a.m),start(a.m);
  for(int k=0;k<a.m;++k) {
    if(clock.expired()){out.seconds=clock.elapsed();return out;}
    // b nonempty blocks require b distinct jobs and b-1 recharges. This
    // incumbent-dependent limit preserves every schedule with value <= UB.
    auto processing=std::vector<int>();for(const auto& p:a.p)processing.push_back(p[k]);
    std::sort(processing.begin(),processing.end());long long work=0;int slots=0;
    const long long deadline=static_cast<long long>(*std::max_element(a.d.begin(),a.d.end()))+out.upper;
    for(int b=1;b<=a.n;++b) {
      work+=processing[b-1];if(work+static_cast<long long>(a.alpha)*(b-1)>deadline)break;
      slots=b;
    }
    slots=std::max(1,slots);z[k].resize(slots);
    Local warm=a.local(k,out.assignment);const auto& ws=out.schedules[k];
    std::vector<int> warm_block(a.n,-1);
    for(int b=0;b<int(ws.blocks.size());++b)for(int j=0;j<warm.n();++j)
      if(ws.blocks[b]&(Mask(1)<<j))warm_block[warm.ids[j]]=b;
    if(int(ws.blocks.size())>slots)throw std::runtime_error("block limit removed warm certificate");
    for(int b=0;b<slots;++b) {
      used[k].push_back(f.model.addVar(0,1,0,GRB_BINARY));
      start[k].push_back(f.model.addVar(0,horizon,0,GRB_CONTINUOUS));
      used[k].back().set(GRB_DoubleAttr_Start,b<int(ws.blocks.size()));
      for(int j=0;j<a.n;++j) {
        z[k][b].push_back(f.model.addVar(0,1,0,GRB_BINARY));
        z[k][b].back().set(GRB_DoubleAttr_Start,warm_block[j]==b);
        f.model.addConstr(z[k][b][j]<=used[k][b]);
      }
    }
    for(int j=0;j<a.n;++j) {
      GRBLinExpr sum=0;for(int b=0;b<slots;++b)sum+=z[k][b][j];
      f.model.addConstr(sum==f.x[j][k]);
    }
    f.model.addConstr(start[k][0]==0);
    for(int b=0;b<slots;++b) {
      if(clock.expired()){out.seconds=clock.elapsed();return out;}
      auto price=f.model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);
      GRBLinExpr energy=a.gamma*price,size=0,time=0;
      for(int j=0;j<a.n;++j) {
        size+=z[k][b][j];time+=a.p[j][k]*z[k][b][j];
        energy+=a.nominal[j][k]/a.capacity[k]*z[k][b][j];
        auto excess=f.model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);energy+=excess;
        f.model.addConstr(excess>=a.deviation[j][k]/a.capacity[k]*z[k][b][j]-price);
      }
      f.model.addConstr(size>=used[k][b]);f.model.addConstr(energy<=used[k][b]);
      if(b+1<slots) {
        f.model.addConstr(used[k][b+1]<=used[k][b]);
        f.model.addConstr(start[k][b+1]==start[k][b]+time+a.alpha*used[k][b+1]);
      }
      GRBLinExpr inside=0;
      for(int j:order) {
        inside+=a.p[j][k]*z[k][b][j];
        f.model.addConstr(start[k][b]+inside<=a.d[j]+f.lambda+
                          std::max(0LL,horizon-a.d[j])*(1-z[k][b][j]));
      }
    }
  }
  if(clock.expired()){out.seconds=clock.elapsed();return out;}
  lb::configure(f.model,clock.left(),0,seed);f.model.optimize();
  const int status=f.model.get(GRB_IntAttr_Status);
  if(status==GRB_INFEASIBLE||status==GRB_UNBOUNDED||status==GRB_INF_OR_UNBD)
    throw std::runtime_error("static model contradicts its explicit feasible certificate");
  out.lower=std::max(out.lower,lb::integer_bound(f.model.get(GRB_DoubleAttr_ObjBound)));
  out.nodes=f.model.get(GRB_DoubleAttr_NodeCount);
  if(f.model.get(GRB_IntAttr_SolCount)>0) {
    Solution candidate;candidate.assignment.assign(a.n,-1);candidate.upper=0;
    for(int j=0;j<a.n;++j)for(int k=0;k<a.m;++k)
      if(f.x[j][k].get(GRB_DoubleAttr_X)>.5)candidate.assignment[j]=k;
    for(int k=0;k<a.m;++k) {
      Local local=a.local(k,candidate.assignment);Local zero=local;std::fill(zero.e.begin(),zero.e.end(),0);
      std::vector<int> sequence;std::vector<bool> cuts;
      for(const auto& block:z[k]) {
        bool first=true;
        for(int j:order)if(block[j].get(GRB_DoubleAttr_X)>.5) {
          cuts.push_back(!sequence.empty()&&first);first=false;
          sequence.push_back(int(std::find(local.ids.begin(),local.ids.end(),j)-local.ids.begin()));
        }
      }
      candidate.schedules.push_back(evaluate(zero,sequence,cuts));
      candidate.upper=std::max(candidate.upper,candidate.schedules.back().lateness);
    }
    validate_static(a,candidate);
    if(candidate.upper<out.upper) {
      out.upper=candidate.upper;out.assignment=candidate.assignment;out.schedules=candidate.schedules;
    }
  }
  if(out.lower>out.upper)throw std::runtime_error("invalid static bound interval");
  if(status==GRB_OPTIMAL&&out.lower!=out.upper)throw std::runtime_error("static optimal status without matching certificate");
  out.status=out.lower==out.upper?"optimal_certified":"time_limit";
  out.seconds=clock.elapsed();return out;
}

// Fix only the nominal-model assignment; order and recharge positions remain
// adjustable. Global budget concentration gives max_k of the local robust
// optima, each with the original Gamma (never Gamma/m).
inline Solution evaluate_assignment(const Instance& a,const std::vector<int>& x,
                                    GRBEnv& env,Clock& clock,int robust_lower) {
  Solution out=static_incumbent(a,x);out.lower=std::max(0,robust_lower);
  std::vector<int> lower(a.m),upper(a.m);std::vector<std::vector<Schedule>> pools(a.m);
  for(int k=0;k<a.m;++k) {
    Local local=a.local(k,x);lower[k]=std::max(0,longrun::workload_lower(local));
    upper[k]=std::max(0,out.schedules[k].lateness);pools[k].push_back(out.schedules[k]);
  }
  auto refresh=[&]() {
    out.lower=std::max(out.lower,*std::max_element(lower.begin(),lower.end()));
    out.upper=*std::max_element(upper.begin(),upper.end());
    if(out.lower>out.upper)throw std::runtime_error("invalid assignment evaluation bounds");
  };
  refresh();std::string unresolved;
  while(out.lower<out.upper&&!clock.expired()) {
    const int target=out.lower+(out.upper-out.lower-1)/2;
    bool all=true;
    std::vector<int> machines(a.m);std::iota(machines.begin(),machines.end(),0);
    std::stable_sort(machines.begin(),machines.end(),[&](int k,int h){return upper[k]>upper[h];});
    for(int k:machines) {
      if(upper[k]<=target)continue;
      if(lower[k]>target){all=false;break;}
      if(clock.expired()){all=false;break;}
      Local local=a.local(k,x);
      RobustOptions options;options.workload_screen=true;options.refine=true;options.reuse_schedules=true;
      options.screening=local.gamma>=local.n();
      auto r=verify(local,target,env,options,clock,&pools[k]);
      ++out.verification_calls;out.deterministic_calls+=r.det_calls;out.adversary_calls+=r.adversary_calls;
      if(r.status==Status::Feasible)upper[k]=target;
      else if(r.status==Status::Infeasible){lower[k]=target+1;all=false;break;}
      else {unresolved=r.reason.empty()?"oracle_unknown":r.reason;all=false;break;}
    }
    refresh();if(all)out.upper=std::min(out.upper,target);
    if(!unresolved.empty())break;
  }
  out.status=out.lower==out.upper?"optimal_certified":
    (unresolved.empty()?"time_limit":"verification_unresolved_"+unresolved);
  // An adaptive UB is a verifier certificate, not the timing of the initial
  // static schedule. Do not export that initial schedule as an adaptive policy.
  out.schedules.clear();out.seconds=clock.elapsed();return out;
}
}
