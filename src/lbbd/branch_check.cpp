// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "model.hpp"
#include "evaluator.hpp"
#include "robust.hpp"
#include "prefix_dual.hpp"
#define BOUNDS_LIBRARY
#include "bounds.cpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <unordered_set>

using namespace lb;

#ifndef AGV_CERTIFICATE_MINING
#define AGV_CERTIFICATE_MINING 0
#endif

#ifndef AGV_MINING_SECONDS
#define AGV_MINING_SECONDS 1.0
#endif

#ifndef AGV_EXPLORATORY
#define AGV_EXPLORATORY 0
#endif

#ifndef AGV_VERIFY_SECONDS
#define AGV_VERIFY_SECONDS 1e9
#endif

#ifndef AGV_OVERLOAD_CUTS
#define AGV_OVERLOAD_CUTS 0
#endif

#ifndef AGV_BRANCH_COVERAGE
#define AGV_BRANCH_COVERAGE 0
#endif

#ifndef AGV_DFF_CUTS
#define AGV_DFF_CUTS 1
#endif

#ifndef AGV_PREFIX_REDUCTION
#define AGV_PREFIX_REDUCTION 1
#endif

#ifndef AGV_LIFT_CONFLICTS
#define AGV_LIFT_CONFLICTS 1
#endif

#ifndef AGV_PREFIX_DUAL_CUTS
#define AGV_PREFIX_DUAL_CUTS 0
#endif

#ifndef AGV_PREFIX_DUAL_SECONDS
#define AGV_PREFIX_DUAL_SECONDS 0.25
#endif

#ifndef AGV_MASTER_FOCUS
#define AGV_MASTER_FOCUS 0
#endif

#ifndef AGV_FIXED_ORDER_SCREEN
#define AGV_FIXED_ORDER_SCREEN 0
#endif

#ifndef AGV_SHARED_BLOCKS
#define AGV_SHARED_BLOCKS 1
#endif

#ifndef AGV_ITERATIVE_RESOLVE
#define AGV_ITERATIVE_RESOLVE 0
#endif

namespace {
Local due_prefix(const Local& a,const std::vector<int>& order,int count) {
  Local b;b.capacity=a.capacity;b.gamma=a.gamma;b.alpha=a.alpha;
  for(int h=0;h<count;++h) {
    int j=order[h];b.ids.push_back(a.ids[j]);b.p.push_back(a.p[j]);b.d.push_back(a.d[j]);
    b.e.push_back(a.e[j]);b.nominal.push_back(a.nominal[j]);b.deviation.push_back(a.deviation[j]);
  }
  return b;
}

Local without_job(const Local& a,int removed) {
  Local b;b.capacity=a.capacity;b.gamma=a.gamma;b.alpha=a.alpha;
  for(int j=0;j<a.n();++j)if(j!=removed) {
    b.ids.push_back(a.ids[j]);b.p.push_back(a.p[j]);b.d.push_back(a.d[j]);b.e.push_back(a.e[j]);
    b.nominal.push_back(a.nominal[j]);b.deviation.push_back(a.deviation[j]);
  }
  return b;
}

int deadline_prefix_bound(const Local& a) {
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  int best=-INF,time=0;Mask mask=0;
  for(int h=0;h<a.n();++h) {
    int j=order[h];time+=a.p[j];mask|=Mask(1)<<j;
    if(h+1<a.n()&&a.d[order[h+1]]==a.d[j])continue;
    int r=std::max(0,int(std::ceil((worst_block_energy(a,mask)-1e-9)/a.capacity))-1);
    int value=time+a.alpha*r-a.d[j];
    best=std::max(best,value);
  }
  return best;
}

struct BranchStats {
  int candidates=0,certified=0,conflicts=0,unknown=0,det_calls=0,robust_calls=0;
  int lifted_conflicts=0,lifted_coefficients=0,lift_gain=0;
  int iterations=0,conflict_jobs_before=0,conflict_jobs_after=0;
  int mined_cuts=0,mined_jobs_before=0,mined_jobs_after=0,mined_q_gain=0;
  int dual_calls=0,dual_complete=0,dual_cuts=0,dual_columns=0;
  double master_seconds=0,verify_seconds=0,mining_seconds=0,dual_seconds=0;
  double nodes=0;
};

struct ConflictRecord {
  int machine=0,q=0,base=0;
  std::vector<int> jobs;
};

struct BranchCallback final:Callback {
  GRBEnv& verify_env;
  int initial_lower;
  int& certified_upper;
  std::vector<int>& best_assignment;
  BranchStats& branch;
  std::vector<std::unordered_map<Mask,std::vector<Schedule>>> pools;
  std::unordered_set<std::string> certified;
  std::ofstream progress,cut_file,unresolved_file;
  int last_progress_lower=-1,last_progress_upper=INF;
  double last_progress_time=-1;
  bool unresolved=false;
  std::string unresolved_reason;

  BranchCallback(Relaxation& form,Clock& timer,const std::string& output_prefix,GRBEnv& sub_env,
                 int lower,int& upper,std::vector<int>& best,BranchStats& stats)
    :Callback(form,timer,AGV_DFF_CUTS,output_prefix+"_packing_trace.csv"),verify_env(sub_env),initial_lower(lower),
     certified_upper(upper),best_assignment(best),branch(stats),pools(form.a.m),
     progress(output_prefix+"_progress.csv"),cut_file(output_prefix+"_cuts.csv"),
     unresolved_file(output_prefix+"_unresolved.csv") {
    scale=true;
    progress<<"seconds,lower,upper,event,candidates,certified,conflicts,unknown\n";
    cut_file<<"machine,q,base,jobs,energies\n";
    unresolved_file<<"target,machine,assignment\n";
    record_progress("initial",true);
  }

  void record_progress(const char* event,bool force=false) {
    const int lower=std::max(initial_lower,integer_bound(bound));
    const double now=clock.elapsed();
    if(force||lower>last_progress_lower||certified_upper<last_progress_upper||now-last_progress_time>=1) {
      progress<<std::setprecision(12)<<now<<','<<lower<<','<<certified_upper<<','<<event<<','
        <<branch.candidates<<','<<branch.certified<<','<<branch.conflicts<<','<<branch.unknown<<'\n';
      progress.flush();last_progress_lower=lower;last_progress_upper=certified_upper;last_progress_time=now;
    }
  }

  std::string key(const std::vector<int>& assignment,int target) const {
    std::string s=std::to_string(target);s.push_back(':');
    for(int k:assignment){s+=std::to_string(k);s.push_back(',');}
    return s;
  }

  Status deterministic_status(const Local& a,int target,Clock& mining) {
    if(mining.expired()||clock.expired())return Status::Unknown;
    Clock limited(std::min(mining.left(),clock.left()));
    auto result=renewed::block_mip(a,target,verify_env,limited);++branch.det_calls;
    return result.status;
  }

  Local minimize_conflict(Local conflict,int target,Clock& mining) {
    std::vector<int> order(conflict.n());std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){
      const int wi=conflict.p[i]+int(std::ceil(conflict.e[i]/std::max(1.0,conflict.capacity)))*conflict.alpha;
      const int wj=conflict.p[j]+int(std::ceil(conflict.e[j]/std::max(1.0,conflict.capacity)))*conflict.alpha;
      if(wi!=wj)return wi>wj;return conflict.ids[i]<conflict.ids[j];
    });
    std::vector<int> ids;for(int i:order)ids.push_back(conflict.ids[i]);
    for(int id:ids) {
      if(conflict.n()<=1||mining.expired()||clock.expired())break;
      int position=-1;for(int j=0;j<conflict.n();++j)if(conflict.ids[j]==id){position=j;break;}
      if(position<0)continue;
      auto trial=without_job(conflict,position);
      if(deterministic_status(trial,target,mining)==Status::Infeasible)conflict=std::move(trial);
    }
    return conflict;
  }

  int strengthen_value(const Local& conflict,int target,Clock& mining) {
    int lower=target+1,upper=edd(conflict).lateness;
    if(upper<lower)return lower;
    while(lower<upper&&!mining.expired()&&!clock.expired()) {
      const int middle=lower+(upper-lower)/2;
      const auto status=deterministic_status(conflict,middle,mining);
      if(status==Status::Infeasible)lower=middle+1;
      else if(status==Status::Feasible)upper=middle;
      else break;
    }
    return lower;
  }

  void add_conflict_cut(const Local& conflict,int k,int q,int current_lower) {
    const int basic=q-current_lower;
    GRBLinExpr rhs=q;
    bool safe_lift=true;
    for(int j=0;j<conflict.n();++j)if(conflict.p[j]-conflict.d[j]>current_lower)safe_lift=false;
    bool strengthened=false;
    for(int j=0;j<conflict.n();++j) {
      int coefficient=basic;
      if(AGV_LIFT_CONFLICTS&&safe_lift)coefficient=std::min(coefficient,conflict.p[j]+f.a.alpha);
      if(coefficient<basic) {
        strengthened=true;++branch.lifted_coefficients;branch.lift_gain+=basic-coefficient;
      }
      rhs-=coefficient*(1-f.x[conflict.ids[j]][k]);
    }
    if(strengthened)++branch.lifted_conflicts;
    addLazy(f.lambda>=rhs);
    ++branch.conflicts;
    cut_file<<k<<','<<q<<','<<current_lower<<",\"";
    for(int j:conflict.ids)cut_file<<j<<' ';
    cut_file<<"\",\"";for(double e:conflict.e)cut_file<<std::setprecision(17)<<e<<' ';
    cut_file<<"\"\n";cut_file.flush();
  }

  void add_prefix_dual_cut(const Local& conflict,int k) {
#if AGV_PREFIX_DUAL_CUTS
    if(conflict.n()==0||clock.expired())return;
    ++branch.dual_calls;
    auto dual=prefix_dual(conflict,verify_env,clock,AGV_PREFIX_DUAL_SECONDS);
    branch.dual_seconds+=dual.seconds;branch.dual_columns+=dual.columns;
    if(!dual.complete)return;
    ++branch.dual_complete;
    const int deadline=*std::max_element(conflict.d.begin(),conflict.d.end());
    const auto position=std::lower_bound(f.dates.begin(),f.dates.end(),deadline);
    if(position==f.dates.end()||*position!=deadline)return;
    const int h=int(position-f.dates.begin());
    double actual=0;GRBLinExpr lhs=0;
    for(int j=0;j<conflict.n();++j) {
      actual+=dual.weight[j]*getSolution(f.x[conflict.ids[j]][k]);
      lhs+=dual.weight[j]*f.x[conflict.ids[j]][k];
    }
    const double rhs=1+getSolution(f.r[k][h]);
    if(actual>rhs+1e-6){addLazy(lhs<=1+f.r[k][h]);++branch.dual_cuts;}
#else
    (void)conflict;(void)k;
#endif
  }

  void callback() override {try {
    const long packing_before=cuts;
    Callback::callback();
    record_progress("master");
    if(failed||where!=GRB_CB_MIPSOL||cuts>packing_before||clock.expired())return;

    ++branch.candidates;
    const int target=int(std::llround(getSolution(f.lambda)));
    std::vector<int> assignment(f.a.n,-1);
    for(int j=0;j<f.a.n;++j)for(int k=0;k<f.a.m;++k)
      if(getSolution(f.x[j][k])>.5)assignment[j]=k;
    const auto candidate_key=key(assignment,target);
    if(certified.count(candidate_key))return;

    std::vector<std::pair<int,int>> ranked_machines;
    for(int k=0;k<f.a.m;++k) {
      auto loc=f.a.local(k,assignment);
      if(loc.n())ranked_machines.push_back({deadline_prefix_bound(loc),k});
    }
    std::stable_sort(ranked_machines.begin(),ranked_machines.end(),[](const auto& x,const auto& y){
      if(x.first!=y.first)return x.first>y.first;return x.second<y.second;
    });
    bool all=true;
    for(const auto& ranked:ranked_machines) {
      if(clock.expired())break;
      const int k=ranked.second;
      Local loc=f.a.local(k,assignment);if(!loc.n())continue;
      Mask mask=0;for(int j:loc.ids)mask|=Mask(1)<<j;
      RobustOptions options;options.workload_screen=true;options.refine=true;options.reuse_schedules=true;
      options.overload_cuts=AGV_OVERLOAD_CUTS;options.branch_coverage=AGV_BRANCH_COVERAGE;
      options.fixed_order=AGV_FIXED_ORDER_SCREEN;options.shared_blocks=AGV_SHARED_BLOCKS;
      options.screening=loc.gamma>=loc.n();
      Clock limited(std::min(clock.left(),double(AGV_VERIFY_SECONDS)));
      Clock timing;auto robust=verify(loc,target,verify_env,options,limited,&pools[k][mask]);
      branch.verify_seconds+=timing.elapsed();++branch.robust_calls;Status state=robust.status;
      Local realized;if(state==Status::Infeasible)realized=realize(loc,robust.witness);
      if(state==Status::Unknown)unresolved_reason=robust.reason;
      if(state==Status::Unknown) {
        ++branch.unknown;
        unresolved_file<<target<<','<<k<<",\"";for(int machine:assignment)unresolved_file<<machine<<' ';
        unresolved_file<<"\"\n";unresolved_file.flush();
#if AGV_EXPLORATORY
        GRBLinExpr selected=0;for(int j=0;j<f.a.n;++j)selected+=f.x[j][assignment[j]];
        addLazy(selected<=f.a.n-1);record_progress("temporary_exclusion",true);return;
#else
        unresolved=true;abort();return;
#endif
      }
      if(state==Status::Feasible)continue;
      all=false;

      std::vector<int> order(loc.n());std::iota(order.begin(),order.end(),0);
      std::stable_sort(order.begin(),order.end(),[&](int i,int j){return loc.d[i]<loc.d[j];});
      int hi=loc.n();
#if AGV_PREFIX_REDUCTION
      int lo=0;
      while(hi-lo>1&&!clock.expired()) {
        int mid=lo+(hi-lo)/2;auto sub=due_prefix(realized,order,mid);
        Clock limited(clock.left());
        auto result=agv::deterministic(sub,target,DetOptions{},limited);++branch.det_calls;
        if(result.status==Status::Infeasible)hi=mid;
        else if(result.status==Status::Feasible)lo=mid;
        else break;
      }
#endif
      Local conflict=due_prefix(realized,order,hi);
      branch.conflict_jobs_before+=loc.n();branch.conflict_jobs_after+=conflict.n();
      add_prefix_dual_cut(conflict,k);
      const int current_lower=AGV_EXPLORATORY?initial_lower:
        std::max(initial_lower,integer_bound(getDoubleInfo(GRB_CB_MIPSOL_OBJBND)));
#if AGV_CERTIFICATE_MINING > 0
      const int original_jobs=conflict.n();
      Clock mining(std::min(double(AGV_MINING_SECONDS),clock.left()));Clock mining_timing;
      auto small=minimize_conflict(conflict,target,mining);
      add_conflict_cut(small,k,target+1,current_lower);
      ++branch.mined_cuts;branch.mined_jobs_before+=original_jobs;branch.mined_jobs_after+=small.n();
#if AGV_CERTIFICATE_MINING > 1
      if(!mining.expired()&&!clock.expired()) {
        const int strengthened=strengthen_value(conflict,target,mining);
        if(strengthened>target+1) {
          auto strong=minimize_conflict(conflict,strengthened-1,mining);
          if(strengthened!=target+1||strong.ids!=small.ids) {
            add_conflict_cut(strong,k,strengthened,current_lower);++branch.mined_cuts;
            branch.mined_jobs_before+=original_jobs;branch.mined_jobs_after+=strong.n();
            branch.mined_q_gain+=strengthened-(target+1);
          }
        }
      }
#endif
      branch.mining_seconds+=mining_timing.elapsed();
#else
      add_conflict_cut(conflict,k,target+1,current_lower);
#endif
      record_progress("conflict",true);
    }
    if(all&&!clock.expired()) {
      certified.insert(candidate_key);++branch.certified;
      if(target<certified_upper){certified_upper=target;best_assignment=assignment;record_progress("incumbent",true);}
    }
  } catch(const std::exception& error) {
    failed=true;unresolved=true;unresolved_reason=std::string("callback_exception_")+error.what();abort();
  } catch(...) {failed=true;unresolved=true;unresolved_reason="callback_exception_unknown";abort();}}
};

enum class IterativeOutcome {Feasible,Cut,Unknown};

struct IterativeSolver {
  Relaxation& master;Clock& clock;GRBEnv& verify_env;int initial_lower;
  int& upper;std::vector<int>& best;BranchStats& stats;
  std::vector<std::unordered_map<Mask,std::vector<Schedule>>> pools;
  Callback packing;
  std::ofstream progress,cut_file,unresolved_file;
  bool unresolved=false;std::string unresolved_reason;
  int first_root_lower=-1;

  IterativeSolver(Relaxation& form,Clock& timer,const std::string& output_prefix,
                  GRBEnv& sub_env,int lower,int& incumbent,std::vector<int>& assignment,
                  BranchStats& statistics)
    :master(form),clock(timer),verify_env(sub_env),initial_lower(lower),upper(incumbent),best(assignment),
     stats(statistics),pools(form.a.m),packing(form,timer,AGV_DFF_CUTS,output_prefix+"_packing_trace.csv"),
     progress(output_prefix+"_progress.csv"),cut_file(output_prefix+"_cuts.csv"),
     unresolved_file(output_prefix+"_unresolved.csv") {
    packing.scale=true;
    progress<<"seconds,lower,upper,event,candidates,certified,conflicts,unknown,iterations\n";
    cut_file<<"machine,q,base,jobs,energies\n";
    unresolved_file<<"target,machine,assignment\n";
  }

  void record(int lower,const char* event) {
    progress<<std::setprecision(12)<<clock.elapsed()<<','<<lower<<','<<upper<<','<<event<<','
      <<stats.candidates<<','<<stats.certified<<','<<stats.conflicts<<','<<stats.unknown<<','
      <<stats.iterations<<'\n';progress.flush();
  }

  Status deterministic_status(const Local& a,int target) {
    if(clock.expired())return Status::Unknown;
    Clock limited(clock.left());auto result=agv::deterministic(a,target,DetOptions{},limited);
    ++stats.det_calls;return result.status;
  }

  void add_conflict(const Local& conflict,int k,int q,int current_lower) {
    const int basic=q-current_lower;GRBLinExpr rhs=q;bool safe_lift=true,strengthened=false;
    for(int j=0;j<conflict.n();++j)if(conflict.p[j]-conflict.d[j]>current_lower)safe_lift=false;
    for(int j=0;j<conflict.n();++j) {
      int coefficient=basic;
      if(AGV_LIFT_CONFLICTS&&safe_lift)coefficient=std::min(coefficient,conflict.p[j]+master.a.alpha);
      if(coefficient<basic){strengthened=true;++stats.lifted_coefficients;stats.lift_gain+=basic-coefficient;}
      rhs-=coefficient*(1-master.x[conflict.ids[j]][k]);
    }
    if(strengthened)++stats.lifted_conflicts;
    master.model.addConstr(master.lambda>=rhs);master.model.update();++stats.conflicts;
    cut_file<<k<<','<<q<<','<<current_lower<<",\"";
    for(int j:conflict.ids)cut_file<<j<<' ';
    cut_file<<"\",\"";for(double e:conflict.e)cut_file<<std::setprecision(17)<<e<<' ';
    cut_file<<"\"\n";cut_file.flush();
  }

  IterativeOutcome verify_candidate(const std::vector<int>& assignment,int target,int current_lower) {
    ++stats.candidates;std::vector<std::pair<int,int>> ranked_machines;
    for(int k=0;k<master.a.m;++k) {
      auto loc=master.a.local(k,assignment);
      if(loc.n())ranked_machines.push_back({deadline_prefix_bound(loc),k});
    }
    std::stable_sort(ranked_machines.begin(),ranked_machines.end(),[](const auto& x,const auto& y){
      if(x.first!=y.first)return x.first>y.first;return x.second<y.second;
    });
    for(const auto& ranked:ranked_machines) {
      if(clock.expired()){unresolved=true;unresolved_reason="time_limit";return IterativeOutcome::Unknown;}
      const int k=ranked.second;Local loc=master.a.local(k,assignment);if(!loc.n())continue;
      Mask mask=0;for(int j:loc.ids)mask|=Mask(1)<<j;
      RobustOptions options;options.workload_screen=true;options.refine=true;options.reuse_schedules=true;
      options.overload_cuts=AGV_OVERLOAD_CUTS;options.branch_coverage=AGV_BRANCH_COVERAGE;
      options.fixed_order=AGV_FIXED_ORDER_SCREEN;options.shared_blocks=AGV_SHARED_BLOCKS;
      options.screening=loc.gamma>=loc.n();Clock limited(std::min(clock.left(),double(AGV_VERIFY_SECONDS)));
      Clock timing;auto robust=verify(loc,target,verify_env,options,limited,&pools[k][mask]);
      stats.verify_seconds+=timing.elapsed();++stats.robust_calls;
      if(robust.status==Status::Unknown) {
        ++stats.unknown;unresolved=true;unresolved_reason=robust.reason;
        unresolved_file<<target<<','<<k<<",\"";for(int machine:assignment)unresolved_file<<machine<<' ';
        unresolved_file<<"\"\n";unresolved_file.flush();return IterativeOutcome::Unknown;
      }
      if(robust.status==Status::Feasible)continue;
      Local realized=realize(loc,robust.witness);
      std::vector<int> order(loc.n());std::iota(order.begin(),order.end(),0);
      std::stable_sort(order.begin(),order.end(),[&](int i,int j){return loc.d[i]<loc.d[j];});
      int hi=loc.n();
#if AGV_PREFIX_REDUCTION
      int lo=0;
      while(hi-lo>1&&!clock.expired()) {
        const int mid=lo+(hi-lo)/2;auto sub=due_prefix(realized,order,mid);
        const auto state=deterministic_status(sub,target);
        if(state==Status::Infeasible)hi=mid;
        else if(state==Status::Feasible)lo=mid;
        else break;
      }
#endif
      Local conflict=due_prefix(realized,order,hi);
      stats.conflict_jobs_before+=loc.n();stats.conflict_jobs_after+=conflict.n();
      add_conflict(conflict,k,target+1,current_lower);return IterativeOutcome::Cut;
    }
    ++stats.certified;upper=target;best=assignment;return IterativeOutcome::Feasible;
  }

  std::pair<int,std::string> solve(int lower,int solver_seed) {
    record(lower,"initial");
    while(lower<upper&&!clock.expired()) {
      ++stats.iterations;lb::configure(master.model,clock.left(),AGV_MASTER_FOCUS,solver_seed);
      master.model.set(GRB_IntParam_LazyConstraints,1);master.model.setCallback(&packing);
      Clock timing;master.model.optimize();stats.master_seconds+=timing.elapsed();
      master.model.setCallback(nullptr);stats.nodes+=master.model.get(GRB_DoubleAttr_NodeCount);
      if(stats.iterations==1) {
        first_root_lower=std::max(initial_lower,integer_bound(packing.root_bound));
        packing.root_open=false;
      }
      if(packing.failed){unresolved=true;unresolved_reason="dff_callback_exception";break;}
      const int status=master.model.get(GRB_IntAttr_Status);
      const double dual=master.model.get(GRB_DoubleAttr_ObjBound);
      if(std::isfinite(dual)&&dual>-INF&&dual<INF)lower=std::max(lower,integer_bound(dual));
      record(lower,"master");
      if(status==GRB_INFEASIBLE){lower=upper;return {lower,"optimal_solver_tolerance"};}
      if(status!=GRB_OPTIMAL||master.model.get(GRB_IntAttr_SolCount)==0||clock.expired())break;
      const int target=int(std::llround(master.model.get(GRB_DoubleAttr_ObjVal)));
      lower=std::max(lower,target);std::vector<int> assignment(master.a.n,-1);
      for(int j=0;j<master.a.n;++j)for(int k=0;k<master.a.m;++k)
        if(master.x[j][k].get(GRB_DoubleAttr_X)>.5)assignment[j]=k;
      const auto outcome=verify_candidate(assignment,target,lower);
      if(outcome==IterativeOutcome::Feasible){lower=target;record(lower,"incumbent");return {lower,"optimal_solver_tolerance"};}
      if(outcome==IterativeOutcome::Unknown)break;
      record(lower,"conflict");
    }
    if(lower>upper)lower=upper;
    return {lower,unresolved?"verification_unresolved_"+unresolved_reason:(clock.expired()?"time_limit":"master_interrupted")};
  }
};
}

int main(int argc,char** argv) {
  try {
    if(argc!=7&&argc!=8)throw std::runtime_error("agv_lbbd INSTANCE TIME OUTPUT_PREFIX WARM_SOLUTION INITIAL_LB SOLVER_SEED [FIXED_TARGET]");
    const auto a=load_instance(argv[1]);const double limit=std::stod(argv[2]);
    const std::string out=argv[3];const int initial_lower=std::stoi(argv[5]);
    const int solver_seed=std::stoi(argv[6]);
    const std::optional<int> fixed_target=argc>=8&&std::string(argv[7])!="-"?std::optional<int>(std::stoi(argv[7])):std::nullopt;
    if(!(limit>0)||!std::isfinite(limit))throw std::runtime_error("time limit must be positive and finite");

    std::ifstream warm(argv[4]);int n=0,file_lower=0,upper=INF;
    if(!(warm>>n>>upper>>file_lower)||n!=a.n||upper<0||upper>=INF)
      throw std::runtime_error("invalid certified warm solution");
    std::vector<int> best(a.n);
    for(int& k:best)if(!(warm>>k)||k<0||k>=a.m)throw std::runtime_error("invalid warm assignment");
    if(initial_lower<0||initial_lower>upper||file_lower!=initial_lower)throw std::runtime_error("invalid or inconsistent initial lower bound");
    if(solver_seed<0||solver_seed>2000000000)throw std::runtime_error("solver seed must be in [0,2000000000]");

    GRBEnv master_env(true),verify_env(true);
    for(GRBEnv* env:{&master_env,&verify_env}) {
      env->set(GRB_IntParam_OutputFlag,0);env->set(GRB_IntParam_Threads,1);env->set(GRB_IntParam_Seed,solver_seed);env->start();
    }
    agv::set_solver_seed(solver_seed);renewed::set_evaluator(&verify_env,1);
    Clock clock(limit);Relaxation master(a,master_env);
    master.lambda.set(GRB_DoubleAttr_LB,initial_lower);
    if(fixed_target) {
      if(*fixed_target<initial_lower||*fixed_target>=upper)
        throw std::runtime_error("fixed target must lie in the open incumbent interval");
      master.lambda.set(GRB_DoubleAttr_LB,*fixed_target);master.lambda.set(GRB_DoubleAttr_UB,*fixed_target);
      master.model.set(GRB_IntParam_SolutionLimit,1);
    } else master.lambda.set(GRB_DoubleAttr_UB,upper-1);
    for(int j=0;j<a.n;++j)for(int k=0;k<a.m;++k)master.x[j][k].set(GRB_DoubleAttr_Start,best[j]==k?1.0:0.0);
    BranchStats stats;
#if AGV_ITERATIVE_RESOLVE
    IterativeSolver iterative(master,clock,out,verify_env,initial_lower,upper,best,stats);
    const auto outcome=iterative.solve(initial_lower,solver_seed);int lower=outcome.first;
    const std::string reason=outcome.second;
    const int root_lower=iterative.first_root_lower<0?initial_lower:iterative.first_root_lower;
#else
    master.model.set(GRB_IntParam_LazyConstraints,1);
    stats.iterations=1;
    BranchCallback callback(master,clock,out,verify_env,initial_lower,upper,best,stats);
    lb::configure(master.model,clock.left(),AGV_MASTER_FOCUS,solver_seed);master.model.setCallback(&callback);
    Clock master_timing;master.model.optimize();stats.master_seconds=master_timing.elapsed();
    master.model.setCallback(nullptr);stats.nodes=master.model.get(GRB_DoubleAttr_NodeCount);
    if(callback.failed&&!callback.unresolved)throw std::runtime_error("branch-and-check callback failed");

    int lower=initial_lower;const int status=master.model.get(GRB_IntAttr_Status);
    const double dual=master.model.get(GRB_DoubleAttr_ObjBound);
    // A dual bound from a model in which lambda is fixed at a trial value is
    // not a global problem bound. It becomes informative only when that
    // target is proved infeasible.
    if(!AGV_EXPLORATORY&&!fixed_target&&std::isfinite(dual)&&dual>-INF&&dual<INF)
      lower=std::max(lower,integer_bound(dual));
    std::string reason;
    if(AGV_EXPLORATORY)reason="exploratory_only";
    else if(callback.unresolved)reason="verification_unresolved_"+callback.unresolved_reason;
    else if(status==GRB_INFEASIBLE) {
      if(fixed_target)lower=std::max(lower,*fixed_target+1);
      else lower=upper;
      reason=fixed_target?"target_infeasible":"optimal_solver_tolerance";
    }
    else if(fixed_target&&upper<=*fixed_target)reason="target_robust_feasible";
    else if(lower>=upper)reason="bounds_match";
    else reason=clock.expired()?"time_limit":"master_status_"+std::to_string(status);
    if(lower>upper)lower=upper;
    callback.record_progress("final",true);
    const int root_lower=std::max(initial_lower,integer_bound(callback.root_bound));
#endif

    std::ofstream result(out+".csv");
    result<<"lower,upper,reason,seconds,architecture,iterations,root_lower,master_seconds,candidates,certified,conflicts,conflict_jobs_before,conflict_jobs_after,unknown,det_calls,robust_calls,verify_seconds,packing_cuts,packing_seconds,nodes,mined_cuts,mined_jobs_before,mined_jobs_after,mined_q_gain,mining_seconds,dual_calls,dual_complete,dual_cuts,dual_columns,dual_seconds,lifted_conflicts,lifted_coefficients,lift_gain\n";
    result<<std::setprecision(12)<<lower<<','<<upper<<','<<reason<<','<<clock.elapsed()<<','
      <<(AGV_ITERATIVE_RESOLVE?"iterative_resolve":"single_tree")<<','<<stats.iterations<<','<<root_lower<<','
      <<(AGV_ITERATIVE_RESOLVE?stats.master_seconds:std::max(0.0,stats.master_seconds-stats.verify_seconds))<<','<<stats.candidates<<','
      <<stats.certified<<','<<stats.conflicts<<','<<stats.conflict_jobs_before<<','<<stats.conflict_jobs_after<<','
      <<stats.unknown<<','<<stats.det_calls<<','<<stats.robust_calls<<','<<stats.verify_seconds<<','
#if AGV_ITERATIVE_RESOLVE
      <<iterative.packing.cuts<<','<<iterative.packing.sep_seconds<<','<<stats.nodes<<','
#else
      <<callback.cuts<<','<<callback.sep_seconds<<','<<stats.nodes<<','
#endif
      <<stats.mined_cuts<<','<<stats.mined_jobs_before<<','<<stats.mined_jobs_after<<','<<stats.mined_q_gain<<','<<stats.mining_seconds<<','
      <<stats.dual_calls<<','<<stats.dual_complete<<','<<stats.dual_cuts<<','<<stats.dual_columns<<','<<stats.dual_seconds<<','
      <<stats.lifted_conflicts<<','<<stats.lifted_coefficients<<','<<stats.lift_gain<<'\n';
    std::ofstream solution(out+".sol");solution<<a.n<<' '<<upper<<' '<<lower<<'\n';
    for(int k:best)solution<<k<<' ';solution<<'\n';
    std::cout<<(AGV_ITERATIVE_RESOLVE?"iterative_lbbd":"branch_check")<<" LB="<<lower<<" UB="<<upper<<" reason="<<reason<<" candidates="<<stats.candidates
      <<" conflicts="<<stats.conflicts<<" seconds="<<clock.elapsed()<<std::endl;
    return 0;
  } catch(const GRBException& e) {
    std::cerr<<e.getMessage()<<'\n';return 1;
  } catch(const std::exception& e) {
    std::cerr<<e.what()<<'\n';return 1;
  }
}
