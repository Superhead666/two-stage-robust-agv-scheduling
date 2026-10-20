// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "../src/runner/model_comparison.hpp"
#include <iostream>
using namespace agv;
using namespace model_comparison;

void require(bool test,const char* message) {if(!test)throw std::runtime_error(message);}
int exhaustive_static_local(const Local& a) {
  if(!a.n())return 0;int best=INF;
  for(const auto& s:enumerate_partitions(a)) {
    bool feasible=true;
    for(auto block:s.blocks)if(worst_block_energy(a,block)>a.capacity+ENERGY_EPS)feasible=false;
    if(feasible)best=std::min(best,std::max(0,s.lateness));
  }
  return best;
}
int exhaustive_adaptive_local(const Local& a,GRBEnv& env) {
  if(!a.n())return 0;int lower=0,upper=exhaustive_static_local(a);Clock clock(20);
  while(lower<upper) {
    int target=lower+(upper-lower-1)/2;auto r=verify_full_partitions(a,target,env,clock);
    require(r.status!=Status::Unknown,"complete partition verification unresolved");
    if(r.status==Status::Feasible)upper=target;else lower=target+1;
  }
  return lower;
}
int exhaustive_static_fleet(const Instance& a) {
  int best=INF;std::vector<int> assignment(a.n);
  std::function<void(int)> recurse=[&](int j) {
    if(j<a.n){for(int k=0;k<a.m;++k){assignment[j]=k;recurse(j+1);}return;}
    int value=0;for(int k=0;k<a.m;++k)value=std::max(value,exhaustive_static_local(a.local(k,assignment)));
    best=std::min(best,value);
  };recurse(0);return best;
}
int main() {try {
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.start();set_solver_seed(1);renewed::set_evaluator(&env,1);
  // Recourse covers the uncertainty set with several schedules; no one fixed
  // two-block partition is robust. Also test Gamma=0 and full-box equivalence.
  Instance trap;trap.n=3;trap.m=2;trap.alpha=10;trap.capacity={10,10};trap.d={3,3,3};
  trap.p={{1,100},{1,100},{1,100}};trap.nominal={{4,4},{4,4},{4,4}};
  trap.deviation={{3,3},{3,3},{3,3}};std::vector<int> x(3,0);
  for(double gamma:{0.,.5,1.,3.}) {
    trap.gamma=gamma;validate_instance(trap);
    const int exact_static=exhaustive_static_fleet(trap);
    const int exact_adaptive=exhaustive_adaptive_local(trap.local(0,x),env);
    Clock static_clock(20);auto s=solve_static(trap,env,static_clock,1,{x},0);
    require(s.lower==exact_static&&s.upper==exact_static,"static counterexample optimum mismatch");
    Clock adaptive_clock(20);auto r=evaluate_assignment(trap,x,env,adaptive_clock,0);
    require(r.lower==exact_adaptive&&r.upper==exact_adaptive,"adaptive counterexample mismatch");
    if(gamma==1)require(s.upper==20&&r.upper==10,"strict recourse benefit lost");
    if(gamma==0||gamma==3)require(s.upper==r.upper,"endpoint model equivalence failed");
    Clock interrupted(0);auto partial=solve_static(trap,env,interrupted,1,{x},0);
    require(partial.lower<=exact_static&&exact_static<=partial.upper,"interrupted static interval invalid");
    Clock interrupted_eval(0);partial=evaluate_assignment(trap,x,env,interrupted_eval,0);
    require(partial.lower<=exact_adaptive&&exact_adaptive<=partial.upper,"interrupted evaluation interval invalid");
  }
  // Unrelated vehicles, fractional budgets, duplicate due dates, zero recharge
  // duration, empty machines, and differing assignment choices.
  for(int test=0;test<12;++test) {
    Instance a=generate(4,2,52310+test,1,.3,1.4,.2);
    a.gamma=std::vector<double>{0,.5,1,1.5,4}[test%5];a.alpha=test%3==0?0:7;
    if(test%2==0)a.d[1]=a.d[0];validate_instance(a);
    std::vector<int> seed(4,0);const int expected=exhaustive_static_fleet(a);
    Clock static_clock(20);auto s=solve_static(a,env,static_clock,1,{seed},0);
    require(s.lower==expected&&s.upper==expected,"unrelated-fleet static MILP mismatch");
    seed={test%2,0,1,0};int fixed_expected=0;
    for(int k=0;k<a.m;++k) {
      auto local=a.local(k,seed);auto blocks=longrun::insertion_reference(local);
      std::vector<int> order;std::vector<bool> cuts;
      for(const auto& block:blocks)for(std::size_t q=0;q<block.size();++q) {
        cuts.push_back(!order.empty()&&q==0);order.push_back(block[q]);
      }
      Local zero=local;std::fill(zero.e.begin(),zero.e.end(),0);
      require(evaluate(zero,order,cuts).lateness==longrun::insertion_upper(local),"production insertion decoder certificate mismatch");
      fixed_expected=std::max(fixed_expected,exhaustive_adaptive_local(local,env));
    }
    Clock eval_clock(20);auto r=evaluate_assignment(a,seed,env,eval_clock,0);
    require(r.lower==fixed_expected&&r.upper==fixed_expected,"fixed-assignment verifier mismatch");
  }
  std::cout<<"PASS 16 exhaustive model comparisons, strict recourse example, endpoints, and interrupted bounds\n";
}catch(const GRBException& e){std::cerr<<e.getMessage()<<'\n';return 1;}
 catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
