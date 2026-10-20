// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "long_search.hpp"
#include "evaluator.hpp"
#include <clocale>
#include <iostream>

using namespace agv;
void require(bool ok,const char* message) {
  if(!ok)throw std::runtime_error(message);
}

namespace longrun {
// Narrow test access to the production comparison and perturbation procedures.
struct SearchValidation {
  static std::vector<int> compare(Search& s,const std::vector<int>& x,const std::vector<int>& incumbent) {
    return s.score(x,true,&incumbent);
  }
  static bool screen(Search& s,const std::vector<int>& x,const std::vector<int>& incumbent) {
    return s.cheap_screen(x,incumbent);
  }
  static Stats stats(const Search& s){return s.stats;}
  static void best_bound(Search& s,const std::vector<int>& bound){s.best_score=bound;}
  static std::vector<int> perturb(Search& s,const std::vector<int>& x) {
    s.remember(x);return s.block_shake(x);
  }
};
}

Instance symmetric() {
  Instance a;a.n=3;a.m=2;a.alpha=10;a.gamma=1;a.capacity={10,10};a.d={14,14,14};
  a.p={{1,1},{1,1},{1,1}};a.nominal={{4,4},{4,4},{4,4}};a.deviation={{3,3},{3,3},{3,3}};
  return a;
}
int exact_local(const Local& a,GRBEnv& env) {
  if(!a.n())return -INF;
  int lo=-INF,hi=std::accumulate(a.p.begin(),a.p.end(),0)+(a.n()-1)*a.alpha-*std::min_element(a.d.begin(),a.d.end());
  for(int j=0;j<a.n();++j)lo=std::max(lo,a.p[j]-a.d[j]);
  Clock clock(20);
  while(lo<hi) {
    int target=lo+(hi-lo)/2;auto r=verify_full_partitions(a,target,env,clock);
    require(r.status!=Status::Unknown,"full-partition reference timed out");
    if(r.status==Status::Feasible)hi=target;else lo=target+1;
  }
  return lo;
}

int main() {try {
  std::setlocale(LC_ALL,"C");
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.set(GRB_IntParam_Threads,1);env.start();
  renewed::set_evaluator(&env,1);
  auto options=longrun::options_for("manuscript");
  require(options.local_limit==.1,"submission HIG verification allowance must remain 0.1 seconds");
  // The longer allowance below is only for the exact certificate checks.
  options.limit=2;options.local_limit=.5;
  require(options.manuscript && options.adaptive && options.targeted,"comparison verification is inactive");
  require(!options.staged && !options.managed && !options.polishing && !options.probe_warm && !options.recombination,"unreported search stages are active");
  require(!options.inheritance && !options.replacement_bounds,"unreported bound transfers are active");
  {
    auto a=symmetric();std::vector<int> x(3,0);
    longrun::Search s(a,env,options);
    require(longrun::SearchValidation::screen(s,x,{-1,-INF}),"workload screening missed a nonimproving candidate");
    require(longrun::SearchValidation::stats(s).entries==0 && longrun::SearchValidation::stats(s).verifications==0,"screening constructed schedules or invoked Verify");
    auto& e=s.entry(0,{0,1,2});require(e.lower==-1 && e.upper==9,"symmetric reference interval changed");
    auto before=longrun::SearchValidation::stats(s).verifications;
    auto value=longrun::SearchValidation::compare(s,x,{-1,-INF});
    require(value.front()==9 && longrun::SearchValidation::stats(s).verifications==before,"screened candidate reached Verify");
    // A certified incumbent bound of 1 requires target 0, even though the
    // intermediate assignment being compared has a leading bound of 5.
    longrun::SearchValidation::best_bound(s,{1,-INF});
    value=longrun::SearchValidation::compare(s,x,{5,-INF});
    require(value.front()==0 && e.upper==0,"verification did not target one below the best incumbent bound");
    require(e.lower==-1 && longrun::SearchValidation::stats(s).verifications==before+1,"feasible target changed the lower bound or used multiple targets");
  }
  {
    auto a=symmetric();a.d={0,0,0};std::vector<int> x(3,0);longrun::Search s(a,env,options);
    auto& e=s.entry(0,{0,1,2});require(e.lower==13 && e.upper==23,"positive target reference changed");
    e.lower=11; // A valid weaker bound isolates the infeasible-return update.
    longrun::SearchValidation::best_bound(s,{13,-INF});
    auto value=longrun::SearchValidation::compare(s,x,{22,-INF});
    require(e.lower==13 && e.upper==23 && value.front()==23,"infeasible best-incumbent target corrupted retained bounds");
    require(longrun::SearchValidation::stats(s).verifications==1,"infeasibility did not stop candidate verification");
    longrun::SearchValidation::compare(s,x,{22,-INF});
    require(longrun::SearchValidation::stats(s).verifications==1,"a retained bound excluding incumbent improvement reached Verify again");
    int lo=e.lower,hi=e.upper;auto calls=longrun::SearchValidation::stats(s).verifications;
    s.refine(e,0,0);
    require(e.lower==lo && e.upper==hi && longrun::SearchValidation::stats(s).verifications==calls,"exhausted call allowance changed bounds");
  }
  std::cout<<"PASS screening and strict best-incumbent verification targets\n";
  {
    // With one machine there is no alternative assignment to compare. The
    // manuscript flow must not refine the incumbent at an unrelated target.
    auto a=symmetric();a.m=1;a.capacity.resize(1);
    for(int j=0;j<a.n;++j){a.p[j].resize(1);a.nominal[j].resize(1);a.deviation[j].resize(1);}
    auto o=longrun::options_for("manuscript");o.limit=.03;o.root_mode=0;
    longrun::Search s(a,env,o);auto r=s.run();
    require(r.stats.verifications==0 && r.upper==9,"incumbent-only verification ran without a candidate comparison");
  }
  std::cout<<"PASS no verification outside candidate comparisons\n";
  {
    Instance a;a.n=8;a.m=2;a.alpha=10;a.gamma=1;a.capacity={15,15};a.d={16,66,0,10,49,16,26,58};
    std::vector<int> p={12,5,19,3,6,4,12,19},energy={6,6,8,10,5,9,4,2};
    for(int j=0;j<a.n;++j){a.p.push_back({p[j],p[j]});a.nominal.push_back({double(energy[j]),double(energy[j])});a.deviation.push_back({.3*energy[j],.3*energy[j]});}
    std::vector<int> x(a.n,0);auto local=a.local(0,x);
    require(longrun::insertion_upper(local)==57,"ordered construction reference changed");
    auto blocks=longrun::insertion_reference(local);
    require(blocks==std::vector<std::vector<int>>{{2},{3},{5},{0,6},{4,7,1}},"unexpected ordered construction schedule");
    longrun::Search s(a,env,options);auto next=longrun::SearchValidation::perturb(s,x);
    require(next[7]==1 && next[5]==0,"perturbation used an unordered construction reference or returned the pivot to its source");
    for(int j:{0,2,3,5,6})require(next[j]==0,"perturbation removed a job outside the selected ordered block");
    require(longrun::SearchValidation::stats(s).verifications==0,"constructive reinsertion invoked Verify");
  }
  std::cout<<"PASS ordered fallback block and pivot exclusion\n";
  int checks=0,assignments=0;
  for(int c=0;c<6;++c) {
    auto a=longrun::read_instance("tests/data/hig/c"+std::to_string(c)+".agv");
    std::map<std::vector<int>,int> values;int optimum=INF;
    for(int bits=0;bits<(1<<a.n);++bits) {
      std::vector<int> x(a.n);for(int j=0;j<a.n;++j)x[j]=(bits>>j)&1;
      int value=-INF;for(int k=0;k<a.m;++k)value=std::max(value,exact_local(a.local(k,x),env));
      values[x]=value;optimum=std::min(optimum,value);++assignments;
    }
    for(unsigned seed:{1u,2u}) {
      auto o=longrun::options_for("manuscript");o.limit=.5;o.seed=seed;
      int observed=INF,root_bound=-INF;double previous_time=-1;bool root_seen=false;
      o.progress=[&](double t,int upper,int lower,const longrun::Stats& stats,const std::vector<int>& x) {
        require(t>=previous_time,"backward progress time");previous_time=t;
        require(std::max(0,upper)<=observed,"increasing upper bound");observed=std::max(0,upper);
        require(upper>=values.at(x) && lower<=optimum,"invalid progress certificate");
        if(stats.mip_seconds>0){if(root_seen)require(root_bound==lower,"global bound changed after initialization");root_bound=lower;root_seen=true;}
        require(stats.pulses==0 && stats.pulse_improvements==0 && stats.pool_calls==0,"unreported incumbent refinement or recombination ran");
        require(stats.replacement_upper==0 && stats.replacement_lower==0 && stats.deletion_updates==0,"unreported bound transfer ran");
      };
      longrun::Search s(a,env,o);auto r=s.run();
      require(r.lower<=optimum && r.assignment_lower<=values.at(r.assignment) && r.upper>=values.at(r.assignment),"invalid final certificate");
      for(int k=0;k<a.m;++k){int value=exact_local(a.local(k,r.assignment),env);require(r.local_lower[k]<=value && value<=r.local_upper[k],"invalid local interval");}
      require(r.stats.negative_targets+r.stats.zero_targets<=r.stats.verifications,"invalid verification counters");
      ++checks;
    }
  }
  std::cout<<"PASS "<<checks<<" manuscript-mode search certificates against "<<assignments<<" exhaustive assignments\n";
  return 0;
}catch(const GRBException& e){std::cerr<<e.getMessage()<<'\n';return 2;}
 catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
