// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "long_search.hpp"
#include <set>

namespace longrun {
void configure(GRBModel& m,double seconds){m.set(GRB_IntParam_OutputFlag,0);m.set(GRB_IntParam_Threads,1);m.set(GRB_DoubleParam_TimeLimit,std::max(.001,seconds));m.set(GRB_DoubleParam_MIPGap,0);m.set(GRB_DoubleParam_MIPGapAbs,1e-9);m.set(GRB_DoubleParam_FeasibilityTol,1e-9);m.set(GRB_DoubleParam_IntFeasTol,1e-9);}
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
  double margin_upper=1;
  struct OverloadRow {GRBConstr row;GRBVar selector;double nominal;};
  std::vector<OverloadRow> rows;
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
    // An adverse realization must violate this schedule. Its common margin
    // cannot exceed the largest individually possible block overload.
    double cap=0;for(Mask b:blocks)cap=std::max(cap,worst_block_energy(a,b)/a.capacity-1);
    if(cap<margin_upper) {
      if(!rows.empty())model.update();
      margin_upper=std::max(0.0,cap);margin.set(GRB_DoubleAttr_UB,margin_upper);
      for(auto& row:rows) {
        double M=std::max(0.0,1+margin_upper-row.nominal);
        // Gurobi has moved the constant nominal energy to the row RHS.
        model.chgCoeff(row.row,row.selector,-M);row.row.set(GRB_DoubleAttr_RHS,1-M-row.nominal);
      }
    }
    if(!seen.insert(blocks).second)return true;
    GRBLinExpr cut=0;
    for(Mask b:blocks){
      if(share && shared.count(b)){cut+=shared.at(b);continue;}
      GRBVar z=model.addVar(0,1,0,GRB_BINARY);GRBLinExpr energy=0;double nominal=0;
      for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j)){nominal+=a.nominal[j];energy+=a.deviation[j]/a.capacity*xi[j];}
      energy+=nominal/a.capacity;
      double bigM=std::max(0.0,1+margin_upper-nominal/a.capacity);
      auto row=model.addConstr(energy-margin-bigM*z>=1-bigM);
      rows.push_back({row,z,nominal/a.capacity});
      if(cost_cuts){
        std::vector<double> ds;GRBLinExpr effort=0;for(int j=0;j<a.n();++j)if(b&(Mask(1)<<j)){ds.push_back(a.deviation[j]);effort+=xi[j];}
        std::sort(ds.rbegin(),ds.rend());double deficit=std::max(0.0,a.capacity-nominal),cost=0;
        for(double d:ds){if(deficit<=0)break;if(d<=0)continue;double fraction=std::min(1.0,deficit/d);cost+=fraction;deficit-=fraction*d;}
        if(deficit>ENERGY_EPS)model.addConstr(z==0);
        else {
          model.addConstr(effort>=std::max(0.0,cost-1e-10)*z);
          model.addConstr(energy-nominal/a.capacity>=std::max(0.0,1-nominal/a.capacity)*z);
          if(share)for(const auto& previous:shared)if(!(previous.first&b) && costs.at(previous.first)+cost>a.gamma+1e-8)model.addConstr(previous.second+z<=1);
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
RobustResult tight_verify(const Local& a,int lambda,GRBEnv& env,const RobustOptions& opt,Clock& budget,std::vector<Schedule>* pool) {
  Clock timing;RobustResult out;
  auto finish=[&](Status status,const std::string& reason){out.status=status;out.reason=reason;out.seconds=timing.elapsed();return out;};
  if(a.n()==0)return finish(Status::Feasible,"empty");
  auto solve=[&](const Local& b){auto r=deterministic(b,lambda,opt.det,budget);++out.det_calls;out.det_nodes+=r.nodes;out.det_seconds+=r.seconds;return r;};
  if(opt.workload_screen && workload_witness(a,lambda,out.witness))return finish(Status::Infeasible,"robust_workload");
  auto nominal=solve(a);if(nominal.status==Status::Unknown)return finish(Status::Unknown,"scheduling_limit");
  if(nominal.status==Status::Infeasible){out.witness.assign(a.n(),0);return finish(Status::Infeasible,"nominal");}
  if(a.gamma<=0)return finish(Status::Feasible,"zero_budget");
  if(opt.screening){Local box=a;for(int j=0;j<a.n();++j)box.e[j]+=std::min(1.0,a.gamma)*a.deviation[j];
    auto upper=solve(box);if(upper.status==Status::Feasible)return finish(Status::Feasible,"box_schedule");
    if(budget.expired())return finish(Status::Unknown,"scheduling_limit");
    if(upper.status==Status::Infeasible && a.gamma>=a.n()){out.witness.assign(a.n(),1);return finish(Status::Infeasible,"full_budget_box");}}
  Coverage coverage(a,env,opt.enhanced,opt.shared_blocks,opt.overload_cuts);
  if(pool && opt.reuse_schedules)for(const auto& s:*pool)if(s.lateness<=lambda)if(!coverage.add(s))return finish(Status::Feasible,"cached_robust_schedule");
  auto add=[&](Schedule s){if(opt.refine)s=refine_schedule(a,s,lambda);
    if(pool && opt.reuse_schedules){bool duplicate=false;for(const auto& old:*pool)if(old.blocks==s.blocks && old.lateness<=s.lateness){duplicate=true;break;}if(!duplicate)pool->push_back(s);}return coverage.add(s);};
  if(!add(nominal.schedule))return finish(Status::Feasible,"robust_schedule");
  while(!budget.expired() && coverage.count<opt.max_schedules){std::vector<double> point;Clock mt;
    int state=coverage.candidate(budget,opt.early_adversary,point,out.margin);out.mip_seconds+=mt.elapsed();++out.adversary_calls;out.schedules=coverage.count;out.blocks=coverage.block_count;
    if(state==0)return finish(Status::Feasible,"coverage_solver_tolerance");
    if(state<0)return finish(Status::Unknown,budget.expired()?"adversary_limit":"numerical_boundary");
    auto next=solve(realize(a,point));if(next.status==Status::Unknown)return finish(Status::Unknown,"scheduling_limit");
    if(next.status==Status::Infeasible){out.witness=point;return finish(Status::Infeasible,"adverse_scenario");}
    if(!add(next.schedule))return finish(Status::Feasible,"robust_schedule");
  }return finish(Status::Unknown,"verification_limit");
}
}
