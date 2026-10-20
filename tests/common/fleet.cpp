// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "robust.hpp"
#include <iostream>

namespace agv {
static std::pair<int,std::vector<int>> initial_assignment(const Instance& a){
  std::vector<int> order(a.n),assignment(a.n,-1),t(a.m,0);std::vector<double> battery=a.capacity;
  std::iota(order.begin(),order.end(),0);std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  int bound=-INF;
  for(int j:order){int best=0,bestfinish=INF;for(int k=0;k<a.m;++k){double e=a.nominal[j][k]+a.deviation[j][k];int f=t[k]+a.p[j][k]+(e>battery[k]+ENERGY_EPS?a.alpha:0);if(f<bestfinish){bestfinish=f;best=k;}}
    double e=a.nominal[j][best]+a.deviation[j][best];if(e>battery[best]+ENERGY_EPS)battery[best]=a.capacity[best];battery[best]-=e;t[best]=bestfinish;assignment[j]=best;bound=std::max(bound,bestfinish-a.d[j]);}
  return {bound,assignment};
}
static std::pair<int,std::vector<int>> improve_assignment(const Instance& a,std::vector<int> assignment,double seconds){
  Clock budget(seconds);std::vector<std::unordered_map<Mask,int>> cache(a.m);
  auto score=[&](const std::vector<int>& x){std::vector<int> values;
    for(int k=0;k<a.m;++k){Mask set=0;for(int j=0;j<a.n;++j)if(x[j]==k)set|=Mask(1)<<j;
      auto found=cache[k].find(set);if(found!=cache[k].end()){values.push_back(found->second);continue;}
      Local b=a.local(k,x);std::vector<int> order(b.n());std::iota(order.begin(),order.end(),0);std::stable_sort(order.begin(),order.end(),[&](int i,int j){return b.d[i]<b.d[j];});
      Mask block=0;int t=0,value=-INF;for(int j:order){Mask bit=Mask(1)<<j;if(worst_block_energy(b,block|bit)>b.capacity+ENERGY_EPS){block=0;t+=b.alpha;}block|=bit;t+=b.p[j];value=std::max(value,t-b.d[j]);}
      cache[k][set]=value;values.push_back(value);}
    std::sort(values.rbegin(),values.rend());return values;};
  auto best=score(assignment);
  while(!budget.expired()){auto next=assignment;auto candidate=best;
    for(int j=0;j<a.n&&!budget.expired();++j)for(int k=0;k<a.m;++k)if(k!=assignment[j]){auto x=assignment;x[j]=k;auto value=score(x);if(value<candidate){candidate=value;next=x;}}
    for(int i=0;i<a.n&&!budget.expired();++i)for(int j=i+1;j<a.n;++j)if(assignment[i]!=assignment[j]){auto x=assignment;std::swap(x[i],x[j]);auto value=score(x);if(value<candidate){candidate=value;next=x;}}
    if(candidate>=best)break;best=candidate;assignment=next;
  }
  return {best.front(),assignment};
}
static Local shrink_conflict(const Local& a,const std::vector<double>& witness,int lambda,double seconds){
  Clock budget(seconds);Local current=realize(a,witness);std::vector<int> candidates=current.ids;
  std::stable_sort(candidates.begin(),candidates.end(),[&](int i,int j){auto first=std::find(a.ids.begin(),a.ids.end(),i)-a.ids.begin(),second=std::find(a.ids.begin(),a.ids.end(),j)-a.ids.begin();return a.d[first]>a.d[second];});
  for(int remove:candidates){if(budget.expired()||current.n()<=1)break;Local next;next.capacity=a.capacity;next.gamma=a.gamma;next.alpha=a.alpha;
    for(int j=0;j<current.n();++j)if(current.ids[j]!=remove){next.ids.push_back(current.ids[j]);next.p.push_back(current.p[j]);next.d.push_back(current.d[j]);next.nominal.push_back(current.nominal[j]);next.deviation.push_back(current.deviation[j]);next.e.push_back(current.e[j]);}
    auto result=deterministic(next,lambda,DetOptions{},budget);if(result.status==Status::Infeasible)current=next;
  }return current;
}
class ScenarioMaster {
public:
  const Instance& a;GRBModel model;GRBVar lambda;std::vector<std::vector<GRBVar>> x;int upper;bool enhanced;int scenario_count=0;
  ScenarioMaster(const Instance& data,GRBEnv& env,int ub,bool improve,bool robust_workload,bool cumulative,const std::vector<int>& start):a(data),model(env),upper(ub),enhanced(improve){
    configure(model,1e9);int lower=-INF;for(int j=0;j<a.n;++j)lower=std::max(lower,*std::min_element(a.p[j].begin(),a.p[j].end())-a.d[j]);
    lambda=model.addVar(lower,upper,1,GRB_INTEGER,"lateness");x.resize(a.n);
    for(int j=0;j<a.n;++j){GRBLinExpr assign=0;for(int k=0;k<a.m;++k){auto v=model.addVar(0,1,0,GRB_BINARY);v.set(GRB_DoubleAttr_Start,start[j]==k);x[j].push_back(v);assign+=v;}model.addConstr(assign==1);}
    if(enhanced){std::vector<int> ds=a.d;std::sort(ds.begin(),ds.end());ds.erase(std::unique(ds.begin(),ds.end()),ds.end());std::vector<GRBVar> previous(a.m);std::vector<bool> has_previous(a.m,false);std::vector<int> previous_due(a.m,-INF);
      for(int k=0;k<a.m;++k)for(int d:ds){GRBLinExpr work=0;for(int j=0;j<a.n;++j)if(a.d[j]<=d)work+=a.p[j][k]*x[j][k];model.addConstr(work<=d+lambda);
        if(robust_workload){
          int eligible=0;GRBLinExpr count=0;for(int j=0;j<a.n;++j)if(a.d[j]<=d){++eligible;count+=x[j][k];}
          auto recharges=model.addVar(0,cumulative?eligible-1:a.n-1,0,GRB_INTEGER);auto price=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);GRBLinExpr energy=a.gamma*price;
          for(int j=0;j<a.n;++j)if(a.d[j]<=d){auto excess=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);model.addConstr(excess>=a.deviation[j][k]*x[j][k]-price);energy+=a.nominal[j][k]*x[j][k]+excess;}
          model.addConstr(energy<=a.capacity[k]*(1+recharges));model.addConstr(work+a.alpha*recharges<=d+lambda);
          if(cumulative){for(int j=0;j<a.n;++j)if(a.d[j]<=d)model.addConstr(recharges<=count-x[j][k]);
            if(has_previous[k]){GRBLinExpr added=0;for(int j=0;j<a.n;++j)if(a.d[j]>previous_due[k]&&a.d[j]<=d)added+=x[j][k];model.addConstr(recharges>=previous[k]);model.addConstr(recharges<=previous[k]+added);}
            previous[k]=recharges;previous_due[k]=d;has_previous[k]=true;}
        }
      }}
  }
  void add_scenario(const std::vector<double>& xi){
    ++scenario_count;int maxd=*std::max_element(a.d.begin(),a.d.end());
    std::vector<int> order(a.n);std::iota(order.begin(),order.end(),0);std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
    std::vector<GRBVar> completion;for(int j=0;j<a.n;++j){completion.push_back(model.addVar(0,upper+a.d[j],0,GRB_CONTINUOUS));model.addConstr(completion.back()<=a.d[j]+lambda);}
    for(int k=0;k<a.m;++k){int minp=INF;for(int j=0;j<a.n;++j)minp=std::min(minp,a.p[j][k]);
      int count=enhanced?std::min(a.n,(upper+maxd+a.alpha)/(minp+a.alpha)):a.n;count=std::max(1,count);
      std::vector<std::vector<GRBVar>> y(count);std::vector<GRBVar> z,begin;
      for(int b=0;b<count;++b){z.push_back(model.addVar(0,1,0,GRB_BINARY));begin.push_back(model.addVar(0,std::max(0,upper+maxd),0,GRB_CONTINUOUS));
        for(int j=0;j<a.n;++j){bool impossible=enhanced && b*(minp+a.alpha)+a.p[j][k]>a.d[j]+upper;y[b].push_back(model.addVar(0,impossible?0:1,0,GRB_BINARY));}}
      model.addConstr(begin[0]==0);
      for(int j=0;j<a.n;++j){GRBLinExpr assigned=0;for(int b=0;b<count;++b)assigned+=y[b][j];model.addConstr(assigned==x[j][k]);}
      for(int b=0;b<count;++b){GRBLinExpr size=0,energy=0,work=0;
        for(int j=0;j<a.n;++j){size+=y[b][j];energy+=(a.nominal[j][k]+a.deviation[j][k]*xi[j])*y[b][j];work+=a.p[j][k]*y[b][j];}
        model.addConstr(size>=z[b]);model.addConstr(size<=a.n*z[b]);model.addConstr(energy<=a.capacity[k]*z[b]);
        if(b+1<count){model.addConstr(z[b+1]<=z[b]);model.addConstr(begin[b+1]==begin[b]+work+a.alpha*z[b+1]);}
        GRBLinExpr internal=0;
        for(int j:order){internal+=a.p[j][k]*y[b][j];model.addGenConstrIndicator(y[b][j],1,completion[j]-begin[b]-internal,GRB_GREATER_EQUAL,0);}
      }
    }
  }
  void add_conflict(const Local& local,int k,int tested,int lower){
    // If this entire job set remains on k, its robust optimum exceeds tested.
    // Deleting jobs cannot worsen lateness or battery feasibility.
    GRBLinExpr missing=static_cast<int>(local.ids.size());for(int j:local.ids)missing-=x[j][k];
    model.addConstr(lambda>=tested+1-(tested+1-lower)*missing);
  }
  void add_insertion_cut(const Local& local,int k,int value,int lower){
    GRBLinExpr rhs=value;for(int q=0;q<local.n();++q)rhs-=std::min(value-lower,local.p[q]+a.alpha)*(1-x[local.ids[q]][k]);model.addConstr(lambda>=rhs);
  }
};
FleetResult solve_fleet(const Instance& a,GRBEnv& env,const FleetOptions& opt,Clock& budget){
  Clock timing;FleetResult out;auto initial=initial_assignment(a);out.upper=initial.first;out.assignment=initial.second;
  int universal_lower=-INF;for(int j=0;j<a.n;++j)universal_lower=std::max(universal_lower,*std::min_element(a.p[j].begin(),a.p[j].end())-a.d[j]);out.lower=universal_lower;
  if(opt.primal_search){auto better=improve_assignment(a,out.assignment,std::min(.15,budget.left()*.05));if(better.first<out.upper){out.upper=better.first;out.assignment=better.second;}}
  ScenarioMaster master(a,env,out.upper,opt.enhanced,opt.robust_workload,opt.cumulative_workload,out.assignment);
  if(opt.root_bound){Clock root;configure(master.model,std::min(.5,budget.left()*.1));master.model.optimize();out.master_seconds+=root.elapsed();if(master.model.get(GRB_IntAttr_Status)==GRB_INFEASIBLE)throw std::runtime_error("workload relaxation contradicts incumbent");double bound=master.model.get(GRB_DoubleAttr_ObjBound);if(bound>-INF&&bound<INF)out.lower=std::max(out.lower,static_cast<int>(std::ceil(bound-1e-6)));}
  if(!opt.logic_only)master.add_scenario(std::vector<double>(a.n,0));
  if(opt.threshold_master){master.lambda.set(GRB_DoubleAttr_Obj,0);master.model.set(GRB_IntParam_SolutionLimit,1);}
  std::vector<std::unordered_map<Mask,std::vector<Schedule>>> pools(a.m);
  for(int iter=0;iter<opt.max_iterations && !budget.expired();++iter){
    if(out.lower>=out.upper){out.status="optimal_solver_tolerance";break;}
    int target=out.lower+(out.upper-out.lower)/2;
    if(opt.threshold_master)master.lambda.set(GRB_DoubleAttr_UB,target);
    ++out.iterations;Clock mt;configure(master.model,budget.left());master.model.optimize();out.master_seconds+=mt.elapsed();
    int status=master.model.get(GRB_IntAttr_Status);if(status==GRB_INFEASIBLE){if(opt.threshold_master){out.lower=target+1;continue;}throw std::runtime_error("scenario master contradicts robust upper bound");}
    if(!opt.threshold_master){double bound=master.model.get(GRB_DoubleAttr_ObjBound);if(bound>-INF && bound<INF)out.lower=std::max(out.lower,static_cast<int>(std::ceil(bound-1e-6)));}
    if(out.lower>=out.upper){out.status="optimal_solver_tolerance";break;}
    if(master.model.get(GRB_IntAttr_SolCount)==0)break;
    int lambda=opt.threshold_master?target:static_cast<int>(std::llround(master.lambda.get(GRB_DoubleAttr_X)));std::vector<int> assignment(a.n,-1);
    for(int j=0;j<a.n;++j)for(int k=0;k<a.m;++k)if(master.x[j][k].get(GRB_DoubleAttr_X)>.5)assignment[j]=k;
    bool all=true,unknown=false;std::vector<std::vector<double>> counterexamples;
    for(int k=0;k<a.m && !budget.expired();++k){Local local=a.local(k,assignment);Mask assigned=0;for(int j:local.ids)assigned|=Mask(1)<<j;
      auto r=verify(local,lambda,env,opt.robust,budget,opt.robust.reuse_schedules?&pools[k][assigned]:nullptr);
      ++out.verify_calls;out.verify_seconds+=r.seconds;out.adversary_calls+=r.adversary_calls;out.schedules+=r.schedules;
      if(r.status==Status::Unknown){unknown=true;all=false;break;}
      if(r.status==Status::Infeasible){all=false;if(opt.conflict_cuts){auto conflict=shrink_conflict(local,r.witness,lambda,std::min(.05,budget.left()));master.add_conflict(conflict,k,lambda,universal_lower);++out.conflicts;
          if(opt.insertion_cuts){Clock limited(std::min(.025,budget.left()));int value=deterministic_optimum(realize(local,r.witness),DetOptions{},limited);if(value<INF){master.add_insertion_cut(local,k,value,universal_lower);++out.conflicts;}}}
        std::vector<double> point(a.n,0);for(int q=0;q<local.n();++q)point[local.ids[q]]=r.witness[q];counterexamples.push_back(point);if(!opt.enhanced)break;}
    }
    if(budget.expired())all=false;
    if(all){out.upper=std::min(out.upper,lambda);out.assignment=assignment;if(out.lower>=out.upper){out.status="optimal_solver_tolerance";break;}master.lambda.set(GRB_DoubleAttr_UB,out.upper);}
    if(unknown || budget.expired())break;
    if(!opt.logic_only)for(const auto& point:counterexamples)master.add_scenario(point);
  }
  out.scenarios=master.scenario_count;out.seconds=timing.elapsed();if(out.lower>=out.upper)out.status="optimal_solver_tolerance";if(out.status=="unknown")out.status=budget.expired()?"time_limit":"unresolved";return out;
}
int enumerate_fleet(const Instance& a,GRBEnv& env,Clock& budget){
  if(a.n>7)throw std::runtime_error("enumeration fleet size");int best=INF;std::vector<int> assignment(a.n,0);bool done=false;
  while(!done){int value=-INF;for(int k=0;k<a.m;++k){Local b=a.local(k,assignment);if(!b.n())continue;
      int lo=-INF,hi=edd(realize(b,std::vector<double>(b.n(),1))).lateness;for(int j=0;j<b.n();++j)lo=std::max(lo,b.p[j]-b.d[j]);
      while(lo<hi){int mid=lo+(hi-lo)/2;auto r=verify_full_partitions(b,mid,env,budget);if(r.status==Status::Unknown)return INF;if(r.status==Status::Feasible)hi=mid;else lo=mid+1;}value=std::max(value,lo);}
    best=std::min(best,value);int j=0;while(j<a.n && ++assignment[j]==a.m){assignment[j]=0;++j;}done=j==a.n;if(budget.expired())return INF;
  }return best;
}
}
