// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "model.hpp"
#include "packing.hpp"
#include "polymatroid.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace lb;
struct Callback:GRBCallback {
  Relaxation& f;Clock& clock;bool pack,poly=false,root_only=false,scale=false;std::ofstream trace;
  std::vector<std::map<Mask,PackResult>> cache;
  double bound=0,root_bound=0,incumbent=GRB_INFINITY,last_time=-1,sep_seconds=0;
  long calls=0,cuts=0,skipped=0,hits=0;int last_bound=-1;bool failed=false;
  bool root_open=true;
  Callback(Relaxation& form,Clock& timer,bool packing,const std::string& path):f(form),clock(timer),pack(packing),trace(path),cache(f.a.m) {
    trace<<"seconds,lower,raw_bound,relaxation_incumbent,cuts,sep_seconds\n";record(true);
  }
  void record(bool force=false) {
    const double t=clock.elapsed();const int b=integer_bound(bound);
    if(force || b>last_bound || t-last_time>=1) {
      trace<<std::setprecision(12)<<t<<','<<b<<','<<bound<<',';
      if(incumbent<INF)trace<<incumbent;
      trace<<','<<cuts<<','<<sep_seconds<<'\n';trace.flush();last_bound=b;last_time=t;
    }
  }
  void update(double b,double u) {
    if(std::isfinite(b)&&b>-INF&&b<INF)bound=std::max(bound,b);
    if(std::isfinite(u)&&u>-INF&&u<INF)incumbent=std::min(incumbent,u);
  }
  void callback() override {try {
    if(where==GRB_CB_MIP) {
      const double current_bound=getDoubleInfo(GRB_CB_MIP_OBJBND);
      const double nodes=getDoubleInfo(GRB_CB_MIP_NODCNT);
      update(current_bound,getDoubleInfo(GRB_CB_MIP_OBJBST));
      if(root_open&&nodes<.5&&std::isfinite(current_bound))root_bound=std::max(root_bound,current_bound);
      else if(nodes>=.5)root_open=false;
    }
    if(poly&&where==GRB_CB_MIPNODE&&getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL&&(!root_only||getDoubleInfo(GRB_CB_MIPNODE_NODCNT)<.5)) {
      Clock timer;++calls;
      for(int k=0;k<f.a.m;++k) {
        std::vector<double> value(f.a.n),deviation(f.a.n);for(int j=0;j<f.a.n;++j){value[j]=getNodeRel(f.x[j][k]);deviation[j]=f.a.deviation[j][k]/f.a.capacity[k];}
        for(int h=0;h<int(f.dates.size());++h) {
          if(clock.expired()){++skipped;break;}
          auto order=f.prefixes[h];std::stable_sort(order.begin(),order.end(),[&](int i,int j){return value[i]>value[j];});
          auto coefficient=polymatroid_coefficients(deviation,f.a.gamma,order);GRBLinExpr lhs=0;double actual=0;
          for(int j:order){double weight=f.a.nominal[j][k]/f.a.capacity[k]+coefficient[j];lhs+=weight*f.x[j][k];actual+=weight*value[j];}
          if(actual>1+getNodeRel(f.r[k][h])+1e-6){addCut(lhs<=1+f.r[k][h]);++cuts;}
        }
      }
      sep_seconds+=timer.elapsed();
    }
    if(where==GRB_CB_MIPSOL) {
      // MIPSOL candidates may be rejected by our lazy constraints. Do not
      // report their objective as a feasible strengthened-relaxation incumbent.
      update(getDoubleInfo(GRB_CB_MIPSOL_OBJBND),GRB_INFINITY);
      if(pack) {
        Clock timer;++calls;
        std::vector<std::vector<double>> value(f.a.n,std::vector<double>(f.a.m));
        for(int j=0;j<f.a.n;++j)for(int k=0;k<f.a.m;++k)value[j][k]=getSolution(f.x[j][k]);
        for(int k=0;k<f.a.m;++k)for(int h=0;h<int(f.dates.size());++h) {
          if(clock.expired()){++skipped;continue;}
          std::vector<int> jobs;Mask mask=0;
          for(int j:f.prefixes[h])if(value[j][k]>.5){jobs.push_back(j);mask|=Mask(1)<<j;}
          if(jobs.empty())continue;
          auto it=cache[k].find(mask);
          if(it==cache[k].end())it=cache[k].emplace(mask,packing_bound(f.a,k,jobs,clock,scale)).first;
          else ++hits;
          const auto& b=it->second;if(!b.complete)++skipped;
          const double rhs=1+getSolution(f.r[k][h]);
          if(!b.q || b.value<=rhs+1e-6)continue;
          GRBLinExpr lhs=0;double actual=0;
          for(int j:f.prefixes[h]) {
            const double w=transformed_size((f.a.nominal[j][k]+f.a.deviation[j][k]*b.xi[j])/f.a.capacity[k],b.q);
            lhs+=w*f.x[j][k];actual+=w*value[j][k];
          }
          if(actual>rhs+1e-6){addLazy(lhs<=1+f.r[k][h]);++cuts;}
        }
        sep_seconds+=timer.elapsed();
      }
    }
    record();
  }catch(...){failed=true;abort();}}
};
#ifndef BOUNDS_LIBRARY
int main(int argc,char** argv){try {
  if(argc!=5)throw std::runtime_error("bounds INSTANCE original|volume|focus|pack|packscale|packfocus|polyroot|polytree|cumulative TIME OUTPUT_PREFIX");
  const auto a=load_instance(argv[1]);const std::string mode=argv[2],out=argv[4];const double limit=std::stod(argv[3]);
  const bool pack=mode=="pack"||mode=="packfocus"||mode=="packscale";const int focus=mode=="focus"||mode=="packfocus"?3:0;
  const bool poly=mode=="polyroot"||mode=="polytree";
  if(mode!="original"&&mode!="volume"&&mode!="focus"&&mode!="cumulative"&&!pack&&!poly)throw std::runtime_error("unknown mode");
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.start();Clock clock(limit);
  Relaxation f(a,env,mode!="original");Callback cb(f,clock,pack,out+"_trace.csv");
  cb.scale=mode=="packscale";
  if(mode=="cumulative")for(int k=0;k<a.m;++k)for(int h=0;h<int(f.dates.size());++h) {
    GRBLinExpr count=0,added=0;
    for(int j:f.prefixes[h]){count+=f.x[j][k];if(h==0||a.d[j]>f.dates[h-1])added+=f.x[j][k];}
    for(int j:f.prefixes[h])f.model.addConstr(f.r[k][h]<=count-f.x[j][k]);
    if(h){f.model.addConstr(f.r[k][h]>=f.r[k][h-1]);f.model.addConstr(f.r[k][h]<=f.r[k][h-1]+added);}
  }
  cb.poly=poly;cb.root_only=mode=="polyroot";if(poly)f.model.set(GRB_IntParam_PreCrush,1);
  if(pack)f.model.set(GRB_IntParam_LazyConstraints,1);
  configure(f.model,clock.left(),focus);f.model.setCallback(&cb);f.model.optimize();
  if(cb.failed)throw std::runtime_error("separation callback failed");
  const int status=f.model.get(GRB_IntAttr_Status),count=f.model.get(GRB_IntAttr_SolCount);
  if(status==GRB_INFEASIBLE||status==GRB_UNBOUNDED||status==GRB_INF_OR_UNBD)throw std::runtime_error("invalid bound model");
  cb.update(f.model.get(GRB_DoubleAttr_ObjBound),count?f.model.get(GRB_DoubleAttr_ObjVal):GRB_INFINITY);cb.record(true);
  if(count) {
    // This is a relaxation assignment, not a robust feasible schedule/UB.
    std::ofstream sol(out+".sol");sol<<a.n<<' '<<INF<<' '<<integer_bound(cb.bound)<<'\n';
    for(int j=0;j<a.n;++j){int chosen=-1;for(int k=0;k<a.m;++k)if(f.x[j][k].get(GRB_DoubleAttr_X)>.5)chosen=k;sol<<chosen<<' ';}
    sol<<"\nRELAXATION_ASSIGNMENT_NOT_A_ROBUST_UPPER_BOUND\n";
  }
  std::ofstream result(out+".csv");
  result<<"mode,lower,relaxation_incumbent,status,seconds,cuts,sep_calls,sep_seconds,sep_skipped,cache_hits,nodes\n";
  result<<std::setprecision(12)<<mode<<','<<integer_bound(cb.bound)<<',';
  if(count)result<<f.model.get(GRB_DoubleAttr_ObjVal);
  result<<','<<status<<','<<clock.elapsed()<<','<<cb.cuts<<','<<cb.calls<<','<<cb.sep_seconds<<','<<cb.skipped<<','<<cb.hits<<','<<f.model.get(GRB_DoubleAttr_NodeCount)<<'\n';
  std::cout<<mode<<" LB="<<integer_bound(cb.bound)<<" incumbent="<<cb.incumbent<<" cuts="<<cb.cuts<<" seconds="<<clock.elapsed()<<std::endl;
}catch(const GRBException& e){std::cerr<<e.getMessage()<<'\n';return 1;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
