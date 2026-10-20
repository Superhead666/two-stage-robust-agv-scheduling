// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "model.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace lb;
struct Path {double profit=0;std::vector<int> count;bool complete=true;};
// Non-elementary paths may repeat a job. This is a relaxation, not a feasible
// AGV schedule. Integer nominal energies are represented exactly; other
// energies are rounded DOWN, preserving the lower-bound direction.
Path price(const Instance& a,int k,int target,const std::vector<double>& weights,Clock& clock,Mask tracked=0) {
  const int horizon=*std::max_element(a.d.begin(),a.d.end())+target;
  const int capacity=a.alpha?int(std::floor(a.capacity[k])):0,width=capacity+1;
  Path out;out.count.assign(a.n,0);if(horizon<0)return out;
  std::vector<int> bit(a.n,-1);int tracked_count=0;for(int j=0;j<a.n;++j)if((tracked>>j)&1)bit[j]=tracked_count++;
  if(tracked_count>=30){out.complete=false;return out;}
  const int labels=1<<tracked_count;
  const std::size_t states=std::size_t(horizon+1)*width*labels;
  if(states>10000000){out.complete=false;return out;}
  std::vector<int> energy(a.n);for(int j=0;j<a.n;++j)energy[j]=a.alpha?int(std::floor(a.nominal[j][k])):0;
  std::vector<double> dp(states,-1);std::vector<int> parent(states,-1),job(states,-1);dp[0]=0;int best=0;
  for(int t=0;t<=horizon;++t) {
    if(clock.expired()){out.complete=false;return out;}
    for(int e=0;e<=capacity;++e)for(int mask=0;mask<labels;++mask) {
      const int s=(t*width+e)*labels+mask;if(dp[s]<0)continue;
      auto update=[&](int next,int j,double profit){if(profit>dp[next]+1e-12){dp[next]=profit;parent[next]=s;job[next]=j;if(profit>dp[best])best=next;}};
      if(e>0&&a.alpha>0&&t+a.alpha<=horizon)update(((t+a.alpha)*width)*labels+mask,-1,dp[s]);
      for(int j=0;j<a.n;++j)if(weights[j]>0&&t+a.p[j][k]<=a.d[j]+target&&e+energy[j]<=capacity&&(bit[j]<0||!(mask&(1<<bit[j]))))
        update(((t+a.p[j][k])*width+e+energy[j])*labels+(bit[j]<0?mask:mask|(1<<bit[j])),j,dp[s]+weights[j]);
    }
  }
  out.profit=dp[best];for(int s=best;parent[s]>=0;s=parent[s])if(job[s]>=0)++out.count[job[s]];
  return out;
}
struct Totals {int targets=0,iterations=0,calls=0,columns=0,refinements=0,max_tracked=0;double pricing=0,master=0;};
struct PathProof {std::vector<double> weights;double upper=0,margin=0;};
int test(const Instance& a,int target,GRBEnv& env,Clock& clock,Totals& totals,PathProof* proof=nullptr,bool refine=false) {
  GRBModel master(env);std::vector<GRBConstr> cover,machine;
  for(int j=0;j<a.n;++j){auto slack=master.addVar(0,GRB_INFINITY,1,GRB_CONTINUOUS);cover.push_back(master.addConstr(slack>=1));}
  for(int k=0;k<a.m;++k)machine.push_back(master.addConstr(GRBLinExpr(0)<=1));
  std::vector<std::set<std::vector<int>>> seen(a.m);
  struct Column {int vehicle;std::vector<int> counts;GRBVar var;};std::vector<Column> columns;std::vector<Mask> tracked(a.m,0);
  auto add=[&](int k,const std::vector<int>& counts){if(!seen[k].insert(counts).second)return false;GRBColumn col;col.addTerm(1,machine[k]);for(int j=0;j<a.n;++j)if(counts[j])col.addTerm(counts[j],cover[j]);auto z=master.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS,col);columns.push_back({k,counts,z});++totals.columns;return true;};
  for(int k=0;k<a.m;++k)for(int j=0;j<a.n;++j)if(a.p[j][k]<=a.d[j]+target){std::vector<int> c(a.n,0);c[j]=1;add(k,c);}
  while(!clock.expired()) {
    lb::configure(master,clock.left());Clock mt;master.optimize();totals.master+=mt.elapsed();++totals.iterations;
    if(master.get(GRB_IntAttr_Status)!=GRB_OPTIMAL)return 0;
    if(master.get(GRB_DoubleAttr_ObjVal)<=1e-7) {
      if(!refine)return 1;
      std::vector<std::vector<double>> score(a.m,std::vector<double>(a.n,0));
      for(const auto& col:columns){double z=col.var.get(GRB_DoubleAttr_X);if(z>1e-9)for(int j=0;j<a.n;++j)score[col.vehicle][j]+=z*std::max(0,col.counts[j]-1);}
      int vehicle=-1,chosen=-1;double best=1e-7;
      for(int k=0;k<a.m;++k)for(int j=0;j<a.n;++j)if(!((tracked[k]>>j)&1)&&score[k][j]>best){vehicle=k;chosen=j;best=score[k][j];}
      if(chosen<0)return 1;
      tracked[vehicle]|=Mask(1)<<chosen;++totals.refinements;totals.max_tracked=std::max(totals.max_tracked,__builtin_popcountll(tracked[vehicle]));
      for(auto& col:columns)if(col.vehicle==vehicle&&col.counts[chosen]>1)col.var.set(GRB_DoubleAttr_UB,0);
      continue;
    }
    std::vector<double> weights(a.n),v(a.m);double certificate=0;
    for(int j=0;j<a.n;++j){weights[j]=std::clamp(cover[j].get(GRB_DoubleAttr_Pi),0.0,1.0);certificate+=weights[j];}
    for(int k=0;k<a.m;++k)v[k]=machine[k].get(GRB_DoubleAttr_Pi);
    int added=0;double total_upper=0;
    for(int k=0;k<a.m;++k) {
      Clock pt;auto path=price(a,k,target,weights,clock,tracked[k]);totals.pricing+=pt.elapsed();++totals.calls;if(!path.complete)return 0;
      // Upward numerical guard on a maximization DP used as a dual bound.
      double upper=path.profit+1e-8*(1+a.n+*std::max_element(a.d.begin(),a.d.end())+target);certificate-=upper;total_upper+=upper;
      if(path.profit+v[k]>1e-7)added+=add(k,path.count);
    }
    if(certificate>1e-6){if(proof){proof->weights=weights;proof->upper=total_upper;proof->margin=certificate;}return -1;}
    if(!added)return 0;
  }
  return 0;
}
#ifndef PATH_BOUND_LIBRARY
int main(int argc,char** argv){try {
  if(argc!=4&&argc!=5)throw std::runtime_error("path_bound INSTANCE TIME OUTPUT_PREFIX [refine]");
  const auto a=load_instance(argv[1]);const double limit=std::stod(argv[2]);const std::string out=argv[3];
  const bool refine=argc==5&&std::string(argv[4])=="refine";
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.start();Clock clock(limit);Totals totals;
  int lower=0;for(int j=0;j<a.n;++j)lower=std::max(lower,*std::min_element(a.p[j].begin(),a.p[j].end())-a.d[j]);
  int upper=0;for(int j=0;j<a.n;++j)upper+=a.p[j][0]+a.alpha;upper-=*std::min_element(a.d.begin(),a.d.end());upper=std::max(upper,lower);
  {Relaxation root(a,env);lb::configure(root.model,std::min(1.0,clock.left()));root.model.optimize();double dual=root.model.get(GRB_DoubleAttr_ObjBound);if(std::isfinite(dual)&&dual>-INF&&dual<INF)lower=std::max(lower,integer_bound(dual));}
  std::ofstream trace(out+"_trace.csv");trace<<"seconds,lower,path_upper,target,result,pricing_calls,columns\n";
  trace<<std::setprecision(12)<<clock.elapsed()<<','<<lower<<','<<upper<<",-1,initial,0,0\n";
  // Ascending integer thresholds avoid spending the entire clock establishing
  // fractional feasibility far above the first infeasibility certificate.
  while(lower<upper&&!clock.expired()) {
    const int target=lower;const int state=test(a,target,env,clock,totals,nullptr,refine);++totals.targets;
    if(state<0)lower=target+1;else if(state>0)upper=target;
    trace<<clock.elapsed()<<','<<lower<<','<<upper<<','<<target<<','<<state<<','<<totals.calls<<','<<totals.columns<<'\n';trace.flush();
    if(state>=0)break;
  }
  std::ofstream file(out+".csv");file<<"mode,lower,path_upper,closed,seconds,targets,iterations,pricing_calls,columns,pricing_seconds,master_seconds\n";
  file<<std::setprecision(12)<<(refine?"path_refined,":"path,")<<lower<<','<<upper<<','<<(lower==upper)<<','<<clock.elapsed()<<','<<totals.targets<<','<<totals.iterations<<','<<totals.calls<<','<<totals.columns<<','<<totals.pricing<<','<<totals.master<<'\n';
  std::cout<<"path LB="<<lower<<" seconds="<<clock.elapsed()<<std::endl;
}catch(const GRBException& e){std::cerr<<e.getMessage()<<'\n';return 1;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
