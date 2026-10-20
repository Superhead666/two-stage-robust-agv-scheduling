// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "long_search.hpp"
#include "evaluator.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>
using namespace agv;

void require(bool test,const char* message) {if(!test)throw std::runtime_error(message);}
int exact_local(const Local& a,GRBEnv& env) {
  if(!a.n())return -INF;
  int lo=-INF,hi=std::accumulate(a.p.begin(),a.p.end(),0)+(a.n()-1)*a.alpha-*std::min_element(a.d.begin(),a.d.end());
  for(int j=0;j<a.n();++j)lo=std::max(lo,a.p[j]-a.d[j]);
  Clock clock(10);
  while(lo<hi) {
    int mid=lo+(hi-lo)/2;auto r=verify_full_partitions(a,mid,env,clock);
    require(r.status!=Status::Unknown,"complete-partition reference timed out");
    if(r.status==Status::Feasible)hi=mid;else lo=mid+1;
  }
  return lo;
}
std::vector<int> clipped(std::vector<int> v) {for(int& x:v)x=std::max(0,x);return v;}
int main(int argc,char** argv) {try {
  std::ofstream expected;
  if(argc>1) {
    std::filesystem::create_directories(argv[1]);
    expected.open(std::string(argv[1])+"/expected.csv");expected<<"instance,optimum\n";
  }
  std::vector<std::vector<int>> vectors;
  for(int i=-2;i<=2;++i)for(int j=-2;j<=i;++j)for(int k=-2;k<=j;++k)vectors.push_back({i,j,k});
  int pairs=0;
  for(const auto& x:vectors)for(const auto& y:vectors) {
    auto tx=clipped(x),ty=clipped(y);
    require(!(x<y) || tx<=ty,"positive-part projection reverses ranking");
    require((x<y)==(std::make_pair(tx,x)<std::make_pair(ty,y)),"signed tie-breaking equivalence");++pairs;
  }
  std::cout<<"PASS "<<pairs<<" sorted-vector comparisons\n";
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.start();renewed::set_evaluator(&env,2);
  {
    int comparisons=0,portfolios=0;
    for(int c=0;c<1200;++c) {
      int n=1+c%50;auto a=longrun::generate_paper(n,2,780000+c,0,.3,2,.45);
      a.gamma=std::vector<double>{0,.5,std::min(2.5,double(n)),double(n),n-.5}[c%5];
      auto b=a.local(c%2,std::vector<int>(n,c%2));
      require(longrun::workload_lower(b)==longrun::incremental_lower(b),"incremental bound mismatch");++comparisons;
      for(const auto& s:longrun::nominal_portfolio(b)){require(valid_schedule(b,s,s.lateness),"invalid nominal seed");++portfolios;}
    }
    std::cout<<"PASS "<<comparisons<<" incremental bounds and "<<portfolios<<" nominal schedules\n";
    Local b;b.p={1,1,1};b.d={14,14,14};b.ids={0,1,2};b.e=b.nominal={4,4,4};b.deviation={3,3,3};b.capacity=10;b.gamma=1;b.alpha=10;
    auto pool=longrun::nominal_portfolio(b);require(pool.size()==3,"symmetric portfolio should contain all three pairs");
    require(longrun::static_upper(b)==9,"static portfolio example");
    RobustOptions ro;ro.reuse_schedules=true;Clock clock(5);auto r=verify(b,-1,env,ro,clock,&pool);
    require(r.status==Status::Feasible,"portfolio coverage example");
    std::cout<<"PASS nominal portfolio certifies -1 while static UB is 9\n";
    b.d={10,11,12};b.e=b.nominal={5,6,5};b.deviation={0,0,0};b.gamma=0;b.alpha=5;
    require(brute_optimum(b)==-3,"noncontiguous EDD example");
    Local moderate;moderate.p.assign(10,1);moderate.d.assign(10,0);moderate.ids.resize(10);std::iota(moderate.ids.begin(),moderate.ids.end(),0);
    moderate.e=moderate.nominal=std::vector<double>(10,1);moderate.deviation.assign(10,.3);moderate.capacity=2.15;moderate.alpha=10;moderate.gamma=1;
    auto policies=longrun::budget_repair_portfolio(moderate);int value=-INF;for(const auto& s:policies)value=std::max(value,s.lateness);
    require(policies.size()==5 && value==60 && longrun::static_upper(moderate)==100,"budget certificate counterexample");
    std::cout<<"PASS five half-budget repair policies certify 60 versus static 100\n";
  }
  {
    Instance a;a.n=2;a.m=2;a.alpha=10;a.gamma=1;a.capacity={100,100};a.d={10,11};
    a.p={{9,1},{10,100}};a.nominal={{1,1},{1,1}};a.deviation={{0,0},{0,0}};
    for(int mode=0;mode<3;++mode) {
      auto o=longrun::options_for("staged");o.limit=.1;o.root_mode=0;o.stop_zero=true;
      o.zero_aware_refinement=mode>0;o.tardiness_score=mode==2;
      longrun::Search search(a,env,o);auto r=search.run();
      require(std::max(0,r.initial)==(mode==2?8:0),"greedy signed-guidance example");
      require(r.upper<=0,"small greedy example not repaired");
    }
    std::cout<<"PASS signed construction Tmax=0; clipped construction Tmax=8 (later repaired)\n";
  }
  {
    Instance a;a.n=3;a.m=1;a.alpha=10;a.gamma=1;a.capacity={10};a.d={14,14,14};
    a.p={{1},{1},{1}};a.nominal={{4},{4},{4}};a.deviation={{3},{3},{3}};
    auto local=a.local(0,{0,0,0});require(exact_local(local,env)==-1,"adaptive on-time reference");
    require(longrun::static_upper(local)==9,"static reference bound");
    for(int mode=0;mode<3;++mode) {
      auto o=longrun::options_for("staged");o.limit=2;o.zero_aware_refinement=mode>0;o.tardiness_score=mode==2;
      longrun::Search search(a,env,o);auto r=search.audit({0,0,0});
      require(r.assignment_lower<=-1 && r.upper>=-1,"invalid signed interval");
      require(r.upper==(mode?0:-1),"zero-aware refinement should stop at zero");
      if(mode)require(r.stats.negative_targets==0 && r.assignment_lower<r.upper,"falsely resolved negative lateness");
    }
    std::cout<<"PASS zero-aware recourse certifies Tmax=0 without claiming exact negative lateness\n";
  }
  int assignments=0,checks=0;
  for(int c=0;c<6;++c) {
    auto a=longrun::generate_paper(4+c%2,2,318900+c,.5+.5*(c%4),.6,1.5+.5*(c%3),.4);
    if(c==4)for(int& d:a.d)d+=20;
    if(c==5)a.alpha=0;
    std::map<std::vector<int>,int> values;int optimum=INF;
    for(int bits=0;bits<(1<<a.n);++bits) {
      std::vector<int> x(a.n);for(int j=0;j<a.n;++j)x[j]=(bits>>j)&1;
      int value=-INF;for(int k=0;k<a.m;++k){auto local=a.local(k,x);int exact=exact_local(local,env);value=std::max(value,exact);int certificate=-INF;for(const auto& s:longrun::budget_repair_portfolio(local))certificate=std::max(certificate,s.lateness);require(certificate>=exact,"invalid budget repair certificate");int improved=-INF;for(const auto& s:longrun::budget_repair_portfolio(local,true))improved=std::max(improved,s.lateness);require(improved>=exact && improved<=certificate,"invalid improved budget certificate");}
      values[x]=value;optimum=std::min(optimum,value);++assignments;
    }
    if(expected.is_open()) {
      const std::string name="c"+std::to_string(c)+".agv";
      save_instance(a,std::string(argv[1])+"/"+name);
      expected<<name<<','<<std::max(0,optimum)<<'\n';expected.flush();
    }
    for(int mode=0;mode<10;++mode)for(unsigned seed:{1u,2u}) {
      auto o=longrun::options_for("staged");o.limit=.2;o.root_mode=0;o.seed=seed;o.stop_zero=true;
      o.managed=true;o.fast_screen=true;o.polishing=true;o.polish_reserve=.02;o.warm_fraction=.1;
      o.focus=true;o.witness_memory=true;o.probe_warm=true;o.replacement_bounds=true;
      o.block_removal=true;o.recourse_reference=true;o.ejection_chains=true;o.chain_after=true;
      o.incremental_lower=mode==1 || mode==6;o.seed_portfolio=mode==2 || mode==6;
      o.critical_prefix=mode==3;o.lex_acceptance=mode==4;o.interval_targets=mode==5;
      if(mode==7){o.incremental_lower=true;o.memo_lower=true;}
      o.budget_portfolio=mode>=8;o.budget_best=mode==9;
      int observed=INF;double observed_time=-1;
      o.progress=[&](double t,int u,int l,const longrun::Stats&,const std::vector<int>& x){
        require(t>=observed_time,"backward progress time");
        require(std::max(0,u)<=observed,"increasing progress upper bound");
        require(std::max(0,l)<=std::max(0,optimum),"invalid progress lower bound");
        require(std::max(0,u)>=std::max(0,values.at(x)),"invalid progress upper bound");
        observed=std::max(0,u);observed_time=t;
      };
      longrun::Search search(a,env,o);auto r=search.run();int exact=values.at(r.assignment);
      require(observed==std::max(0,r.upper),"final progress does not match result");
      require(r.assignment_lower<=exact && r.upper>=exact && r.lower<=optimum,"invalid search certificate");
      ++checks;
    }
  }
  std::cout<<"PASS "<<checks<<" search certificates against "<<assignments<<" exhaustive assignments\n";
  return 0;
}catch(const GRBException& e){std::cerr<<e.getMessage()<<'\n';return 2;}
 catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
