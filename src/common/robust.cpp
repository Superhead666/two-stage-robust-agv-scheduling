// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "robust.hpp"
#include <set>

namespace agv {
namespace {thread_local int configured_solver_seed=1;}
void set_solver_seed(int seed){configured_solver_seed=seed;}
void configure(GRBModel& m,double seconds){m.set(GRB_IntParam_OutputFlag,0);m.set(GRB_IntParam_Threads,1);m.set(GRB_DoubleParam_TimeLimit,std::max(.001,seconds));m.set(GRB_DoubleParam_MIPGap,0);m.set(GRB_DoubleParam_MIPGapAbs,1e-9);m.set(GRB_DoubleParam_FeasibilityTol,1e-9);m.set(GRB_DoubleParam_IntFeasTol,1e-9);m.set(GRB_IntParam_Seed,configured_solver_seed);}
class PositiveMargin:public GRBCallback {
  GRBVar margin;
  void callback() override {try{if(where==GRB_CB_MIPSOL && getSolution(margin)>1e-6)abort();}catch(const GRBException&) {}}
public:explicit PositiveMargin(GRBVar v):margin(v){}
};
class Coverage {
public:
  const Local& a;GRBModel model;std::vector<GRBVar> xi;GRBVar margin;
  std::unordered_map<Mask,GRBVar> shared;
  std::set<std::vector<Mask>> seen;
  bool enhanced,share,cost_cuts;int count=0,block_count=0;
  std::unordered_map<Mask,double> costs;
  Coverage(const Local& data,GRBEnv& env,bool strengthen,bool reuse,bool cuts=false):a(data),model(env),enhanced(strengthen),share(reuse),cost_cuts(cuts){
    configure(model,1e9);GRBLinExpr budget=0;
    for(int j=0;j<a.n();++j){xi.push_back(model.addVar(0,1,0,GRB_CONTINUOUS));budget+=xi.back();}
    model.addConstr(budget<=a.gamma);margin=model.addVar(0,1,1,GRB_CONTINUOUS);model.set(GRB_IntAttr_ModelSense,GRB_MAXIMIZE);
  }
  bool add(const Schedule& s){
    std::vector<Mask> blocks;
    for(Mask b:s.blocks)if(!enhanced || worst_block_energy(a,b)>a.capacity+ENERGY_EPS)blocks.push_back(b);
    if(blocks.empty())return false;
    std::sort(blocks.begin(),blocks.end());
    if(!seen.insert(blocks).second)return true;
    GRBLinExpr cut=0;
    for(Mask b:blocks){
      if(share && shared.count(b)){cut+=shared.at(b);continue;}
      GRBVar z=model.addVar(0,1,0,GRB_BINARY);GRBLinExpr energy=0;double nominal=0;
      for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j)){nominal+=a.nominal[j];energy+=a.deviation[j]/a.capacity*xi[j];}
      energy+=nominal/a.capacity;
      double bigM=std::max(0.0,2-nominal/a.capacity);
      model.addConstr(energy-margin>=1-bigM*(1-z));
      if(cost_cuts){
        std::vector<double> ds;GRBLinExpr effort=0;for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j)){ds.push_back(a.deviation[j]);effort+=xi[j];}
        std::sort(ds.rbegin(),ds.rend());double deficit=std::max(0.0,a.capacity-nominal),cost=0;
        for(double d:ds){if(deficit<=0)break;if(d<=0)continue;double fraction=std::min(1.0,deficit/d);cost+=fraction;deficit-=fraction*d;}
        if(deficit>ENERGY_EPS)model.addConstr(z==0);
        else {
          model.addConstr(effort>=std::max(0.0,cost-1e-10)*z);
          model.addConstr(energy-nominal/a.capacity>=std::max(0.0,1-nominal/a.capacity)*z);
          if(share)for(const auto& previous:shared)if(!(previous.first&b) && costs.at(previous.first)+cost>a.gamma+1e-8)
            model.addConstr(previous.second+z<=1);
        }
        costs[b]=cost;
      }
      if(share)shared.emplace(b,z);cut+=z;++block_count;
    }
    model.addConstr(cut>=1);++count;return true;
  }
  // 1 = a strictly uncovered candidate, 0 = numerical coverage certificate, -1 = unresolved.
  int candidate(Clock& budget,bool early,std::vector<double>& point,double& value){
    configure(model,budget.left());PositiveMargin callback(margin);if(early)model.setCallback(&callback);
    model.optimize();model.setCallback(nullptr);int status=model.get(GRB_IntAttr_Status);
    if(status==GRB_INFEASIBLE)return 0;
    if(model.get(GRB_IntAttr_SolCount)>0){value=margin.get(GRB_DoubleAttr_X);
      if(value>1e-8){point.clear();double sum=0;for(auto v:xi){point.push_back(std::max(0.0,std::min(1.0,v.get(GRB_DoubleAttr_X))));sum+=point.back();}
        if(sum>a.gamma && sum>0)for(double& v:point)v*=a.gamma/sum;
        // Independently check every stored schedule after correcting solver roundoff.
        for(const auto& schedule:seen){double violation=-1e100;
          for(Mask b:schedule){double e=0;for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j))e+=a.nominal[j]+a.deviation[j]*point[j];violation=std::max(violation,e-a.capacity);}
          if(violation<=ENERGY_EPS)return -1;
        }return 1;
      }
    }
    if(status==GRB_OPTIMAL && model.get(GRB_DoubleAttr_ObjBound)<=1e-9)return 0;
    return -1;
  }
  int candidate_lp(Clock& budget,std::vector<double>& point,double& value){
    GRBModel lp(model.getEnv());configure(lp,budget.left());lp.set(GRB_IntAttr_ModelSense,GRB_MAXIMIZE);
    std::vector<GRBVar> variables;GRBLinExpr used=0;for(int j=0;j<a.n();++j){variables.push_back(lp.addVar(0,1,0,GRB_CONTINUOUS));used+=variables.back();}lp.addConstr(used<=a.gamma);
    auto delta=lp.addVar(0,1,1,GRB_CONTINUOUS);std::vector<Mask> selected;std::set<std::vector<Mask>> exhausted;
    auto demand=[&](Mask b,const std::vector<double>& x){double total=0;for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j))total+=a.nominal[j]+a.deviation[j]*x[j];return total;};
    std::function<int()> search=[&](){
      if(budget.expired()||selected.size()>512)return -1;
      auto key=selected;std::sort(key.begin(),key.end());if(exhausted.count(key))return 0;
      configure(lp,budget.left());lp.optimize();int status=lp.get(GRB_IntAttr_Status);if(status==GRB_INFEASIBLE){exhausted.insert(key);return 0;}if(status!=GRB_OPTIMAL)return -1;
      if(lp.get(GRB_DoubleAttr_ObjVal)<=1e-9){exhausted.insert(key);return 0;}
      std::vector<double> x;double sum=0;for(auto v:variables){x.push_back(std::clamp(v.get(GRB_DoubleAttr_X),0.0,1.0));sum+=x.back();}if(sum>a.gamma&&sum>0)for(double& v:x)v*=a.gamma/sum;
      for(Mask b:selected)if(demand(b,x)<=a.capacity+ENERGY_EPS)return -1;
      const std::vector<Mask>* branch=nullptr;
      for(const auto& schedule:seen){bool hit=false;for(Mask b:schedule)if(demand(b,x)>a.capacity+ENERGY_EPS){hit=true;break;}
        if(!hit&&(!branch||schedule.size()<branch->size()))branch=&schedule;}
      if(!branch){point=x;value=delta.get(GRB_DoubleAttr_X);return 1;}
      auto choices=*branch;std::stable_sort(choices.begin(),choices.end(),[&](Mask b,Mask c){return demand(b,x)>demand(c,x);});
      for(Mask b:choices){GRBLinExpr energy=0;for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j))energy+=a.nominal[j]/a.capacity+a.deviation[j]/a.capacity*variables[j];
        auto constraint=lp.addConstr(energy-delta>=1);selected.push_back(b);int result=search();selected.pop_back();lp.remove(constraint);if(result!=0)return result;
      }
      exhausted.insert(key);return 0;
    };
    return search();
  }
};
static bool workload_witness(const Local& a,int lambda,std::vector<double>& xi){
  for(int due:a.d){Mask jobs=0;int work=0;std::vector<int> order;
    for(int j=0;j<a.n();++j)if(a.d[j]<=due){jobs|=Mask(1)<<j;work+=a.p[j];order.push_back(j);}
    int charges=std::max(0,static_cast<int>(std::ceil((worst_block_energy(a,jobs)-ENERGY_EPS)/a.capacity))-1);
    if(work+charges*a.alpha<=due+lambda)continue;
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.deviation[i]>a.deviation[j];});
    xi.assign(a.n(),0);double left=a.gamma;for(int j:order){xi[j]=std::min(1.0,std::max(0.0,left));left-=xi[j];}return true;
  }return false;
}
static Schedule refine_schedule(const Local& a,Schedule s,int lambda){
  // Splitting a block enlarges its feasible energy region. Added recharges are
  // accepted only if the original deadline bound remains satisfied.
  while(true){int position=-1;double best=-1;std::vector<bool> selected;
    for(int q=1;q<a.n();++q)if(!s.recharge[q]){
      auto cuts=s.recharge;cuts[q]=true;auto next=evaluate(a,s.order,cuts);if(next.lateness>lambda)continue;
      double before=0,after=0;for(Mask b:s.blocks)before+=std::max(0.0,worst_block_energy(a,b)-a.capacity);
      for(Mask b:next.blocks)after+=std::max(0.0,worst_block_energy(a,b)-a.capacity);
      double gain=before-after;if(gain>best+1e-10){best=gain;position=q;selected=cuts;}
    }
    if(position<0)break;s=evaluate(a,s.order,selected);
  }return s;
}
static bool find_static_policy(const Local& a,int lambda,double seconds,Schedule& solution){
  Clock budget(seconds);int n=a.n();std::vector<int> priority(n),order;std::vector<bool> cuts;
  std::iota(priority.begin(),priority.end(),0);std::stable_sort(priority.begin(),priority.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  struct Key{Mask done,block;bool operator==(const Key& other)const{return done==other.done&&block==other.block;}};
  struct Hash{std::size_t operator()(const Key& key)const{return std::hash<Mask>{}(key.done)^(std::hash<Mask>{}(key.block)+0x9e3779b97f4a7c15ULL+(key.done<<6)+(key.done>>2));}};
  std::unordered_map<Key,int,Hash> visited;std::unordered_map<Mask,double> energies;std::uint64_t nodes=0;bool expired=false;
  auto energy=[&](Mask b){auto it=energies.find(b);if(it!=energies.end())return it->second;double e=worst_block_energy(a,b);energies.emplace(b,e);return e;};
  Mask full=(Mask(1)<<n)-1;
  std::function<bool(Mask,Mask,int,int)> dfs=[&](Mask done,Mask block,int r,int t){
    if((++nodes&63)==0&&budget.expired()){expired=true;return false;}
    if(done==full){solution=evaluate(a,order,cuts);return true;}
    Key key{done,block};auto it=visited.find(key);if(it!=visited.end()&&it->second<=r)return false;visited[key]=r;
    int work=0;double demand=0,battery=a.capacity-energy(block);
    for(int j:priority)if(!(done&(Mask(1)<<j))){work+=a.p[j];demand+=a.nominal[j];int count=static_cast<int>(std::max(0.0,std::ceil((demand-battery-ENERGY_EPS)/a.capacity)));if(t+work+a.alpha*count>a.d[j]+lambda)return false;}
    for(int j:priority)if(!(done&(Mask(1)<<j))){Mask bit=Mask(1)<<j;bool charge=energy(block|bit)>a.capacity+ENERGY_EPS;int finish=t+a.p[j]+charge*a.alpha;
      if(finish>a.d[j]+lambda)continue;order.push_back(j);cuts.push_back(charge);
      if(dfs(done|bit,charge?bit:block|bit,r+charge,finish))return true;order.pop_back();cuts.pop_back();if(expired)return false;
    }return false;
  };
  return dfs(0,0,0,0);
}
static Schedule ordered_greedy(const Local& a,const std::vector<int>& order){
  double battery=a.capacity;std::vector<bool> recharge;for(int j:order){bool charge=a.e[j]>battery+ENERGY_EPS;recharge.push_back(charge);if(charge)battery=a.capacity;battery-=a.e[j];}return evaluate(a,order,recharge);
}
class OrderCoverage {
public:
  const Local& a;int lambda;bool tight;GRBModel model;GRBVar margin;std::vector<GRBVar> xi;
  std::set<std::vector<int>> orders;std::unordered_map<Mask,GRBVar> shared;
  OrderCoverage(const Local& data,int bound,GRBEnv& env,bool improve=false):a(data),lambda(bound),tight(improve),model(env){
    configure(model,1e9);model.set(GRB_IntAttr_ModelSense,GRB_MAXIMIZE);margin=model.addVar(0,1,1,GRB_CONTINUOUS);
    GRBLinExpr used=0;for(int j=0;j<a.n();++j){xi.push_back(model.addVar(0,1,0,GRB_CONTINUOUS));used+=xi.back();}model.addConstr(used<=a.gamma);
  }
  // Zero means this single order covers all scenarios. Minus one is a duplicate.
  int add(const std::vector<int>& order){
    if(!orders.insert(order).second)return -1;if(a.alpha==0)return 0;
    if(ordered_greedy(a,order).lateness>lambda)throw std::runtime_error("nominally infeasible coverage order");
    std::vector<int> lower(a.n()+1,0),upper(a.n()+1);std::iota(upper.begin(),upper.end(),0);
    if(tight){double low_battery=a.capacity,high_battery=a.capacity;int low_count=1,high_count=1;
      for(int h=1;h<=a.n();++h){int j=order[h-1];double lo=a.nominal[j],hi=lo+std::min(1.0,a.gamma)*a.deviation[j];if(lo>low_battery+ENERGY_EPS){++low_count;low_battery=a.capacity;}if(hi>high_battery+ENERGY_EPS){++high_count;high_battery=a.capacity;}low_battery-=lo;high_battery-=hi;lower[h]=low_count;upper[h]=high_count;}}
    std::vector<GRBVar> potential;for(int h=0;h<=a.n();++h){potential.push_back(model.addVar(lower[h],upper[h],0,GRB_CONTINUOUS));if(tight&&h)model.addConstr(potential[h]>=potential[h-1]);}
    GRBLinExpr choose=0;int work=0,targets=0;
    for(int h=1;h<=a.n();++h){int j=order[h-1];work+=a.p[j];int required=(a.d[j]+lambda-work)/a.alpha+2;
      if(required<=upper[h]){auto target=model.addVar(0,1,0,GRB_BINARY);choose+=target;++targets;model.addConstr(potential[h]>=lower[h]+(required-lower[h])*target);}}
    if(!targets)return 0;if(tight)model.addConstr(choose==1);else model.addConstr(choose>=1);
    for(int u=0;u<a.n();++u){Mask block=0;double nominal=0;GRBLinExpr energy=0;
      for(int v=u+1;v<=a.n();++v){int j=order[v-1];block|=Mask(1)<<j;nominal+=a.nominal[j];energy+=a.nominal[j]/a.capacity+a.deviation[j]/a.capacity*xi[j];
        if(nominal>a.capacity+ENERGY_EPS)break;
        if(worst_block_energy(a,block)<=a.capacity+ENERGY_EPS)model.addConstr(potential[v]<=potential[u]+1);
        else {int bigM=tight?std::max(0,std::min(v-u-1,upper[v]-lower[u]-1)):a.n();if(bigM==0)continue;
          auto found=shared.find(block);GRBVar inactive;
          if(found!=shared.end())inactive=found->second;
          else {inactive=model.addVar(0,1,0,GRB_BINARY);if(tight)for(const auto& previous:shared){if((previous.first&block)==previous.first)model.addConstr(previous.second<=inactive);else if((previous.first&block)==block)model.addConstr(inactive<=previous.second);}shared.emplace(block,inactive);model.addConstr(energy-margin>=1-std::max(0.0,2-nominal/a.capacity)*(1-inactive));}
          model.addConstr(potential[v]<=potential[u]+1+bigM*inactive);
        }
      }
    }return 1;
  }
  int candidate(Clock& budget,std::vector<double>& point,double& value){
    PositiveMargin callback(margin);model.setCallback(&callback);configure(model,budget.left());model.optimize();model.setCallback(nullptr);int status=model.get(GRB_IntAttr_Status);
    if(status==GRB_INFEASIBLE)return 0;
    if(model.get(GRB_IntAttr_SolCount)>0 && margin.get(GRB_DoubleAttr_X)>1e-8){value=margin.get(GRB_DoubleAttr_X);double sum=0;point.clear();for(auto v:xi){point.push_back(std::clamp(v.get(GRB_DoubleAttr_X),0.0,1.0));sum+=point.back();}if(sum>a.gamma&&sum>0)for(double& x:point)x*=a.gamma/sum;
      Local realized=realize(a,point);for(const auto& order:orders)if(ordered_greedy(realized,order).lateness<=lambda)return -1;return 1;
    }
    if(status==GRB_OPTIMAL&&model.get(GRB_DoubleAttr_ObjBound)<=1e-9)return 0;return -1;
  }
};
RobustResult verify_fixed_order(const Local& a,const std::vector<int>& order,int lambda,GRBEnv& env,Clock& budget,bool tight){
  Clock timing;RobustResult out;auto finish=[&](Status status,const std::string& reason){out.status=status;out.reason=reason;out.seconds=timing.elapsed();return out;};
  if(ordered_greedy(a,order).lateness>lambda){out.witness.assign(a.n(),0);return finish(Status::Infeasible,"order_nominal");}
  OrderCoverage coverage(a,lambda,env,tight);if(coverage.add(order)==0)return finish(Status::Feasible,"order_timing");
  int state=coverage.candidate(budget,out.witness,out.margin);return finish(state==0?Status::Feasible:state==1?Status::Infeasible:Status::Unknown,state==0?"fixed_order_certificate":state==1?"order_adverse_scenario":"order_verification_limit");
}
RobustResult verify_order_enumeration(const Local& a,const std::vector<int>& order,int lambda,GRBEnv& env,Clock& budget){
  if(a.n()>10)throw std::runtime_error("order enumeration size");RobustResult out;Coverage coverage(a,env,true,true);Local zero=a;std::fill(zero.e.begin(),zero.e.end(),0);bool any=false;
  for(Mask mask=0;mask<(Mask(1)<<std::max(0,a.n()-1));++mask){std::vector<bool> recharge(a.n());for(int q=1;q<a.n();++q)recharge[q]=(mask>>(q-1))&1;auto schedule=evaluate(zero,order,recharge);
    if(schedule.lateness<=lambda){any=true;if(!coverage.add(schedule)){out.status=Status::Feasible;return out;}}}
  if(!any){out.status=Status::Infeasible;return out;}int state=coverage.candidate(budget,false,out.witness,out.margin);out.status=state==0?Status::Feasible:state==1?Status::Infeasible:Status::Unknown;return out;
}
RobustResult verify(const Local& a,int lambda,GRBEnv& env,const RobustOptions& opt,Clock& budget,std::vector<Schedule>* pool){
  Clock timing;RobustResult out;if(a.n()==0){out.status=Status::Feasible;out.reason="empty";return out;}
  auto solve=[&](const Local& b){auto r=deterministic(b,lambda,opt.det,budget);++out.det_calls;out.det_nodes+=r.nodes;out.det_seconds+=r.seconds;return r;};
  auto finish=[&](Status status,const std::string& why){out.status=status;out.reason=why;out.seconds=timing.elapsed();return out;};
  if(opt.workload_screen && workload_witness(a,lambda,out.witness))return finish(Status::Infeasible,"robust_workload");
  auto nominal=solve(a);if(nominal.status==Status::Unknown)return finish(Status::Unknown,"scheduling_limit");
  if(nominal.status==Status::Infeasible){out.witness.assign(a.n(),0);return finish(Status::Infeasible,"nominal");}
  if(a.gamma<=0)return finish(Status::Feasible,"zero_budget");
  if(opt.screening){
    Local box=a;for(int j=0;j<a.n();++j)box.e[j]=a.nominal[j]+std::min(1.0,a.gamma)*a.deviation[j];
    auto upper=solve(box);if(upper.status==Status::Feasible)return finish(Status::Feasible,"box_schedule");
    if(budget.expired())return finish(Status::Unknown,"scheduling_limit");
    if(upper.status==Status::Infeasible && a.gamma>=a.n()){out.witness.assign(a.n(),1);return finish(Status::Infeasible,"full_budget_box");}
  }
  if(opt.static_policy){Schedule policy;Clock pt;bool found=find_static_policy(a,lambda,std::min(.025,budget.left()),policy);out.det_seconds+=pt.elapsed();
    if(found){for(Mask b:policy.blocks)if(worst_block_energy(a,b)>a.capacity+ENERGY_EPS)throw std::runtime_error("invalid static policy certificate");if(policy.lateness>lambda)throw std::runtime_error("invalid static policy timing");return finish(Status::Feasible,"static_policy_certificate");}}
  if(opt.order_generation){OrderCoverage coverage(a,lambda,env,opt.order_bounds);if(coverage.add(nominal.schedule.order)==0)return finish(Status::Feasible,"order_timing");
    if(pool && opt.reuse_schedules)for(const auto& schedule:*pool)if(ordered_greedy(a,schedule.order).lateness<=lambda)if(coverage.add(schedule.order)==0)return finish(Status::Feasible,"order_timing");
    while(!budget.expired()&&coverage.orders.size()<static_cast<std::size_t>(opt.max_schedules)){
      Clock mt;std::vector<double> point;int state=coverage.candidate(budget,point,out.margin);out.mip_seconds+=mt.elapsed();++out.adversary_calls;out.schedules=coverage.orders.size();out.blocks=coverage.shared.size();
      if(state==0)return finish(Status::Feasible,"order_coverage_certificate");if(state<0)return finish(Status::Unknown,budget.expired()?"adversary_limit":"numerical_boundary");
      auto next=solve(realize(a,point));if(next.status==Status::Unknown)return finish(Status::Unknown,"scheduling_limit");if(next.status==Status::Infeasible){out.witness=point;return finish(Status::Infeasible,"adverse_scenario");}
      int added=coverage.add(next.schedule.order);if(added==0)return finish(Status::Feasible,"order_timing");if(added<0)return finish(Status::Unknown,"duplicate_order_boundary");
      if(pool&&opt.reuse_schedules)pool->push_back(next.schedule);
    }return finish(Status::Unknown,"verification_limit");
  }
  std::vector<double> order_point;
  if(opt.fixed_order){Clock limited(std::min(.2,budget.left()));auto r=verify_fixed_order(a,nominal.schedule.order,lambda,env,limited);out.mip_seconds+=r.seconds;++out.adversary_calls;
    if(r.status==Status::Feasible)return finish(Status::Feasible,r.reason);if(r.status==Status::Infeasible)order_point=r.witness;}
  Coverage coverage(a,env,opt.enhanced,opt.shared_blocks,opt.overload_cuts);
  if(pool && opt.reuse_schedules)for(const auto& s:*pool)if(s.lateness<=lambda)if(!coverage.add(s))return finish(Status::Feasible,"cached_robust_schedule");
  auto add=[&](Schedule s){
    if(opt.refine)s=refine_schedule(a,std::move(s),lambda);
    if(pool && opt.reuse_schedules){bool duplicate=false;for(const auto& old:*pool)if(old.blocks==s.blocks && old.lateness<=s.lateness){duplicate=true;break;}if(!duplicate)pool->push_back(s);}
    return coverage.add(s);
  };
  if(!add(nominal.schedule))return finish(Status::Feasible,"robust_schedule");
  if(!order_point.empty()){auto next=solve(realize(a,order_point));if(next.status==Status::Infeasible){out.witness=order_point;return finish(Status::Infeasible,"order_witness_full_infeasibility");}if(next.status==Status::Feasible&&!add(next.schedule))return finish(Status::Feasible,"robust_schedule");if(budget.expired())return finish(Status::Unknown,"scheduling_limit");}
  while(!budget.expired() && coverage.count<opt.max_schedules){
    std::vector<double> point;Clock mt;int state=opt.branch_coverage?coverage.candidate_lp(budget,point,out.margin):coverage.candidate(budget,opt.early_adversary,point,out.margin);
    out.mip_seconds+=mt.elapsed();++out.adversary_calls;out.schedules=coverage.count;out.blocks=coverage.block_count;
    if(state==0)return finish(Status::Feasible,"coverage_solver_tolerance");
    if(state<0)return finish(Status::Unknown,budget.expired()?"adversary_limit":"numerical_boundary");
    auto next=solve(realize(a,point));
    if(next.status==Status::Unknown)return finish(Status::Unknown,"scheduling_limit");
    if(next.status==Status::Infeasible){out.witness=point;return finish(Status::Infeasible,"adverse_scenario");}
    if(!add(next.schedule))return finish(Status::Feasible,"robust_schedule");
  }return finish(Status::Unknown,"verification_limit");
}
RobustResult verify_full_partitions(const Local& a,int lambda,GRBEnv& env,Clock& budget){
  Clock timing;RobustResult out;Coverage coverage(a,env,true,true);
  auto schedules=enumerate_partitions(a);bool any=false;
  for(const auto& s:schedules)if(s.lateness<=lambda){any=true;if(!coverage.add(s)){out.status=Status::Feasible;out.reason="full_robust_schedule";out.seconds=timing.elapsed();return out;}}
  if(!any){out.status=Status::Infeasible;out.witness.assign(a.n(),0);out.reason="no_deadline_schedule";return out;}
  int state=coverage.candidate(budget,false,out.witness,out.margin);
  out.status=state==0?Status::Feasible:state==1?Status::Infeasible:Status::Unknown;out.schedules=coverage.count;out.blocks=coverage.block_count;out.seconds=timing.elapsed();return out;
}
}
