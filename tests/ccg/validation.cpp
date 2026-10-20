// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"
#include <iostream>

namespace baseline {
namespace {
void require(bool okay, const std::string& message) {
  if (!okay) throw std::runtime_error("VALIDATION FAILED: " + message);
}
int full_local(const Local& local, GRBEnv& env) {
  int lower = 0, upper = std::accumulate(local.p.begin(), local.p.end(), 0)
      + std::max(0,local.n()-1)*local.alpha;
  Clock budget(20);
  while (lower < upper) {
    int target = lower + (upper-lower)/2;
    auto r = verify_full_partitions(local, target, env, budget);
    require(r.status != Status::Unknown, "partition-enumeration timeout");
    if (r.status == Status::Feasible) upper = target;
    else lower = target+1;
  }
  return upper;
}
int full_fleet(const Instance& a, GRBEnv& env, bool robust) {
  std::vector<std::map<Mask,int>> cache(a.m);
  std::vector<int> x(a.n,0); int best=INF; bool done=false;
  while (!done) {
    int value = 0;
    for (int k=0;k<a.m;++k) {
      Mask mask=0;for(int j=0;j<a.n;++j)if(x[j]==k)mask|=Mask(1)<<j;
      auto it=cache[k].find(mask);
      if (it==cache[k].end()) {
        auto local=a.local(k,x);
        int cost=robust?full_local(local,env):std::max(0,brute_optimum(local));
        it=cache[k].emplace(mask,cost).first;
      }
      value=std::max(value,it->second);
    }
    best=std::min(best,value);
    int j=0;while(j<a.n && ++x[j]==a.m){x[j]=0;++j;}done=j==a.n;
  }
  return best;
}
Instance trap(bool adaptive) {
  Instance a; a.n=3;a.m=1;a.alpha=10;a.gamma=1;a.capacity={10};
  a.d=adaptive?std::vector<int>{3,3,3}:std::vector<int>{1,2,3};
  a.p={{1},{1},{1}};a.nominal={{4},{4},{4}};
  a.deviation=adaptive?std::vector<std::vector<double>>{{3},{3},{3}}:
                       std::vector<std::vector<double>>{{4},{4},{4}};
  return a;
}
std::vector<std::string> formulations() {
  std::vector<std::string> result={"block","position"};
#ifdef AGV_WITH_CP
  result.push_back("cp");
#endif
  return result;
}
}
void validate(GRBEnv& env, bool quick) {
  Clock total;std::mt19937 rng(260831);int deterministic_checks=0,master_checks=0,robust_checks=0;
  env.set(GRB_IntParam_Threads,1);env.set(GRB_IntParam_Seed,19);
  set_solver_seed(19);
  {
    GRBModel probe(env);configure(probe,.1);
    require(probe.get(GRB_IntParam_Threads)==1,"every subproblem uses one thread");
    require(probe.get(GRB_IntParam_Seed)==19,"subproblem configure preserves seed setting");
  }
  env.set(GRB_IntParam_Threads,1);env.set(GRB_IntParam_Seed,1);set_solver_seed(1);
  std::cout<<"PASS single-thread solver configuration and seed inheritance\n";
  const auto policy=baseline_verifier_options();
  require(policy.enhanced&&!policy.screening&&policy.shared_blocks&&policy.early_adversary
          &&policy.workload_screen&&policy.refine&&!policy.overload_cuts&&policy.reuse_schedules
          &&!policy.static_policy&&!policy.branch_coverage&&!policy.fixed_order
          &&!policy.order_generation&&!policy.order_bounds,"CCG and LBBD verifier policies differ");
  std::cout<<"PASS common exact-verifier policy\n";
  for(int c=0;c<(quick?12:60);++c) {
    auto a=generate_paper(1+c%6,1,951000+c,0,.3,2,.4);
    auto local=a.local(0,std::vector<int>(a.n,0));
    local.capacity=10;local.alpha=c%7==0?0:5+c%4;
    for(int j=0;j<local.n();++j) {
      local.e[j]=(rng()%41)/4.0;local.nominal[j]=local.e[j];local.deviation[j]=0;
      local.p[j]=1+rng()%7;local.d[j]=rng()%45;
    }
    int truth=brute_optimum(local);
    for(int target:{truth-1,truth,truth+1}) {
      Clock t(10);auto r=deterministic_mip(local,target,env,t);
      require(r.status!=Status::Unknown,"deterministic MIP timeout");
      require((r.status==Status::Feasible)==(target>=truth),"deterministic MILP versus all permutations/recharges");
      if(r.status==Status::Feasible)require(valid_schedule(local,r.schedule,target),"deterministic schedule validation");
      ++deterministic_checks;
    }
  }
  std::cout<<"PASS deterministic checks "<<deterministic_checks<<std::endl;
  for(int c=0;c<(quick?4:16);++c) {
    auto a=generate_paper(3+c%3,2,952000+c,0,0,1.5+c%3,.3);
    if(c%7==0)a.alpha=0;
    if(c%5==0)for(auto& row:a.nominal)for(double& e:row)e=0;
    if(c%6==0)for(int& d:a.d)d+=100;
    int truth=full_fleet(a,env,false);
    for(const std::string& method:formulations()) {
      Options opt;opt.formulation=method;opt.master_target_upper=truth+3;
      Clock budget(20);auto master=make_master(a,env,opt,budget);
      master->add_scenario(std::vector<double>(a.n,0),budget);
      // Repeatedly change the target to exercise both tightening and relaxation.
      double first_extraction=-1;
      for(int target:{truth+1,truth,std::max(0,truth-1),truth+2}) {
        auto r=master->solve(target,budget);
        require(r.status!=Status::Unknown,"finite master timeout");
        require((r.status==Status::Feasible)==(target>=truth),"master formulation "+method+" at target "+std::to_string(target));
        if(r.status==Status::Feasible) {
          int value=0;for(int k=0;k<a.m;++k)value=std::max(value,brute_optimum(a.local(k,r.assignment)));
          require(value<=target,"master assignment independent feasibility");
        }
        if(method=="cp") {
          if(first_extraction<0)first_extraction=master->extraction_time();
          else require(master->extraction_time()==first_extraction,"CP model re-extracted without a new scenario");
        }
        ++master_checks;
      }
    }
  }
  std::cout<<"PASS finite-scenario master checks "<<master_checks<<std::endl;
  int scenario_checks=0;
  for(double gamma:{0.,.5,1.,1.5,3.,4.})for(int mask=0;mask<(quick?4:16);++mask) {
    auto a=generate_paper(4,2,954001,gamma,.6,1.5,0);
    a.alpha=7;std::vector<int> fixed(a.n);
    for(int j=0;j<a.n;++j)fixed[j]=(mask>>j)&1;
    // Two different, budget-feasible scenarios: same assignment, independent
    // sequencing/recharging. Exhaustive permutations are the reference.
    std::vector<std::vector<double>> points(2,std::vector<double>(a.n,0));
    double remaining=gamma;
    for(int j=0;j<a.n;++j){points[0][j]=std::min(1.,remaining);remaining-=points[0][j];}
    for(int j=0;j<a.n;++j)points[1][j]=gamma/a.n;
    int truth=0;
    for(const auto& point:points)for(int k=0;k<a.m;++k) {
      auto local=a.local(k,fixed);std::vector<double> xi;
      for(int j:local.ids)xi.push_back(point[j]);
      truth=std::max(truth,brute_optimum(realize(local,xi)));
    }
    for(const std::string& method:formulations()) {
      Options opt;opt.formulation=method;Clock budget(10);
      auto master=make_master(a,env,opt,budget,&fixed);
      for(const auto& point:points)master->add_scenario(point,budget);
      for(int target:{truth-1,truth}) {
        auto result=master->solve(target,budget);
        require(result.status!=Status::Unknown,"multi-scenario master timeout");
        require((result.status==Status::Feasible)==(target>=truth),"multi-scenario "+method+" versus exhaustive schedules");
        ++scenario_checks;
      }
    }
  }
  std::cout<<"PASS fixed-assignment multi-scenario checks "<<scenario_checks<<std::endl;
  {
    auto identical=generate_paper(6,3,955000,0,.3,2,.4);
    for(int k=1;k<identical.m;++k) {
      identical.capacity[k]=identical.capacity[0];
      for(int j=0;j<identical.n;++j) {
        identical.p[j][k]=identical.p[j][0];
        identical.nominal[j][k]=identical.nominal[j][0];
        identical.deviation[j][k]=identical.deviation[j][0];
      }
    }
    const int truth=full_fleet(identical,env,false);
    for(const std::string& method:formulations()) {
      Options opt;opt.formulation=method;Clock budget(20);
      auto master=make_master(identical,env,opt,budget);
      master->add_scenario(std::vector<double>(identical.n,0),budget);
      auto result=master->solve(truth,budget);
      require(result.status==Status::Feasible,"identical-machine symmetry master feasibility");
      for(int q=1;q<identical.m;++q)for(int j=0;j<identical.n;++j)if(result.assignment[j]==q) {
        bool predecessor=false;for(int i=0;i<j;++i)predecessor|=result.assignment[i]==q-1;
        require(predecessor,"first-use ordering for identical machines");
      }
    }
    std::cout<<"PASS identical-machine first-use symmetry constraints\n";
  }
  for(int c=0;c<(quick?4:10);++c) {
    Instance a;
    if(c<2)a=trap(c==1);
    else {
      int n=4+c%2;
      double gamma=c%4==0?0:c%4==1?.5:c%4==2?1.5:n;
      a=generate_paper(n,2,953000+c,gamma,c%6==0?0:.6,1.5+c%2,.3);
      if(c==8)for(int& d:a.d)d+=200;
      if(c==9)a.alpha=0;
    }
    int truth=full_fleet(a,env,true);
    if(c<2)require(truth==(c==0?20:10),"analytic nonvertex/adaptive-recourse example");
    for(const std::string& method:formulations()) {
      Options opt;opt.formulation=method;opt.limit=30;
      auto r=solve(a,env,opt);
      require(r.initial<INF&&r.certificate!="none","certified singleton initial solution");
      require(r.lower==truth && r.upper==truth,"full robust "+method+" case "+std::to_string(c)+" truth "+std::to_string(truth)+" got ["+std::to_string(r.lower)+","+std::to_string(r.upper)+"] reason "+r.reason);
      for(int k=0;k<a.m;++k)require(full_local(a.local(k,r.assignment),env)<=r.upper,"returned robust assignment");
      ++robust_checks;
    }
  }
  // A vanishing allowance must return honest bounds, never an optimality claim by timeout.
  auto a=trap(false);Options short_run;short_run.limit=1e-9;
  auto limited=solve(a,env,short_run);
  require(limited.lower<=20 && limited.upper>=20,"timeout bound validity");
  require(limited.status=="time_limit","timeout status");
  require(limited.upper==20&&limited.assignment.size()==3&&limited.certificate=="singleton_edd_list",
          "vanishing-time run retains only its certified initial solution");
  std::cout<<"PASS robust solves "<<robust_checks<<", including nonvertex uncertainty, adaptive advantage, zero tardiness, fractional budgets, and timeout handling\n"
           <<"ALL VALIDATION PASSED in "<<total.elapsed()<<" seconds\n";
}
}
