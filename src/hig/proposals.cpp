// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "long_search.hpp"
#include <set>

namespace longrun {
int incremental_lower(const Local& a) {
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  const int g=std::min(a.n(),int(std::floor(a.gamma)));
  const double fraction=a.gamma-std::floor(a.gamma);
  std::multiset<double> top;double sum=0,nominal=0;int processing=0,answer=-INF;
  for(int j:order) {
    top.insert(a.deviation[j]);sum+=a.deviation[j];
    if(int(top.size())>g+1){sum-=*top.begin();top.erase(top.begin());}
    double support=sum;
    if(int(top.size())==g+1)support-=(1-fraction)*(*top.begin());
    nominal+=a.nominal[j];processing+=a.p[j];
    int bins=std::max(1,int(std::ceil((nominal+support-ENERGY_EPS)/a.capacity)));
    answer=std::max(answer,processing+(bins-1)*a.alpha-a.d[j]);
  }
  return answer;
}

std::vector<Schedule> nominal_portfolio(const Local& a) {
  std::vector<Schedule> result;if(!a.n())return result;
  // Same four construction rules as the incumbent. Rotate initial tie order
  // to prevent identical nominal data from producing one repeated partition.
  for(int method=0;method<4;++method) {
    std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
    std::rotate(order.begin(),order.begin()+method%a.n(),order.end());
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){
      if(method==1)return a.d[i]-a.p[i]<a.d[j]-a.p[j];
      if(method==2)return a.d[i]+a.p[i]<a.d[j]+a.p[j];
      if(method==3)return a.nominal[i]+a.deviation[i]>a.nominal[j]+a.deviation[j];
      return a.d[i]<a.d[j];
    });
    std::vector<std::vector<int>> blocks;std::vector<double> energy;
    for(int j:order) {
      int selected=-1;double fill=-1;
      int first=method==0?std::max(0,int(blocks.size())-1):0;
      for(int b=first;b<int(blocks.size());++b)if(energy[b]+a.nominal[j]<=a.capacity+ENERGY_EPS && energy[b]+a.nominal[j]>fill){selected=b;fill=energy[b]+a.nominal[j];}
      if(selected<0){blocks.push_back({j});energy.push_back(a.nominal[j]);}
      else {blocks[selected].push_back(j);energy[selected]+=a.nominal[j];}
    }
    for(auto& b:blocks)std::stable_sort(b.begin(),b.end(),[&](int i,int j){return a.d[i]<a.d[j];});
    auto date=[&](const auto& b){int p=0,l=-INF;for(int j:b){p+=a.p[j];l=std::max(l,p-a.d[j]);}return p-l;};
    std::stable_sort(blocks.begin(),blocks.end(),[&](const auto& b,const auto& c){return date(b)<date(c);});
    std::vector<int> sequence;std::vector<bool> recharge;
    for(const auto& b:blocks)for(int j:b){recharge.push_back(!sequence.empty() && j==b.front());sequence.push_back(j);}
    Local nominal=a;nominal.e=a.nominal;auto s=evaluate(nominal,sequence,recharge);
    if(!valid_schedule(nominal,s,s.lateness))throw std::runtime_error("invalid nominal portfolio schedule");
    auto key=s.blocks;std::sort(key.begin(),key.end());bool duplicate=false;
    for(const auto& old:result){auto previous=old.blocks;std::sort(previous.begin(),previous.end());if(key==previous && old.lateness<=s.lateness)duplicate=true;}
    if(!duplicate)result.push_back(std::move(s));
  }
  return result;
}

namespace {
using Partition=std::vector<std::vector<int>>;
void append_distinct(std::vector<Partition>& pool,const std::vector<Partition>& candidates) {
  auto key=[](Partition p){for(auto& b:p)std::sort(b.begin(),b.end());std::sort(p.begin(),p.end());return p;};
  for(const auto& p:candidates){auto k=key(p);bool seen=false;for(const auto& old:pool)if(key(old)==k){seen=true;break;}if(!seen)pool.push_back(p);}
}
std::vector<Partition> static_candidates(const Local& a) {
  std::vector<Partition> result;std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  for(int method=0;method<4;++method) {
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){
      if(method==1)return a.d[i]-a.p[i]<a.d[j]-a.p[j];
      if(method==2)return a.d[i]+a.p[i]<a.d[j]+a.p[j];
      if(method==3)return a.nominal[i]+a.deviation[i]>a.nominal[j]+a.deviation[j];
      return a.d[i]<a.d[j];
    });
    Partition blocks;
    for(int j:order) {
      int selected=-1;double fill=-1;
      const int first=method==0?std::max(0,int(blocks.size())-1):0;
      for(int b=first;b<int(blocks.size());++b){Mask mask=Mask(1)<<j;for(int h:blocks[b])mask|=Mask(1)<<h;
        double energy=worst_block_energy(a,mask);if(energy<=a.capacity+ENERGY_EPS && energy>fill){fill=energy;selected=b;}}
      if(selected<0)blocks.push_back({j});else blocks[selected].push_back(j);
    }
    result.push_back(std::move(blocks));
  }
  return result;
}
}
std::vector<Schedule> budget_repair_portfolio(const Local& a,bool best_candidates) {
  if(!a.n())return {};
  Local protected_data=a;protected_data.gamma=a.gamma/2;
  std::vector<Partition> bases{static_reference(protected_data)};
  if(best_candidates)append_distinct(bases,static_candidates(protected_data));
  auto schedule=[&](std::vector<std::vector<int>> blocks){
    for(auto& b:blocks)std::stable_sort(b.begin(),b.end(),[&](int i,int j){return a.d[i]<a.d[j];});
    auto date=[&](const auto& b){int p=0,l=-INF;for(int j:b){p+=a.p[j];l=std::max(l,p-a.d[j]);}return p-l;};
    std::stable_sort(blocks.begin(),blocks.end(),[&](const auto& b,const auto& c){return date(b)<date(c);});
    std::vector<int> order;std::vector<bool> recharge;
    for(const auto& b:blocks)for(int j:b){recharge.push_back(!order.empty() && j==b.front());order.push_back(j);}
    Local nominal=a;nominal.e=a.nominal;auto s=evaluate(nominal,order,recharge);
    if(!valid_schedule(nominal,s,s.lateness))throw std::runtime_error("invalid repair portfolio schedule");
    return s;
  };
  std::vector<Schedule> best_result;int best_value=INF;
  for(const auto& base:bases) {
  std::vector<Schedule> result;
  for(int h=0;h<int(base.size());++h) {
    Mask mask=0;for(int j:base[h])mask|=Mask(1)<<j;
    if(worst_block_energy(protected_data,mask)>a.capacity+ENERGY_EPS)throw std::runtime_error("invalid half-budget block");
    if(worst_block_energy(a,mask)<=a.capacity+ENERGY_EPS)continue;
    Local subset;subset.gamma=a.gamma;subset.alpha=a.alpha;subset.capacity=a.capacity;
    for(int j:base[h]){subset.ids.push_back(j);subset.p.push_back(a.p[j]);subset.d.push_back(a.d[j]);subset.e.push_back(a.nominal[j]);subset.nominal.push_back(a.nominal[j]);subset.deviation.push_back(a.deviation[j]);}
    std::vector<Partition> repairs{static_reference(subset)};
    if(best_candidates)append_distinct(repairs,static_candidates(subset));
    Schedule best_repair;best_repair.lateness=INF;
    for(const auto& repair:repairs){
      Partition blocks;
      for(int b=0;b<int(base.size());++b)if(b!=h)blocks.push_back(base[b]);
      for(auto b:repair){for(int& j:b)j=base[h][j];blocks.push_back(std::move(b));}
      auto candidate=schedule(std::move(blocks));if(candidate.lateness<best_repair.lateness)best_repair=std::move(candidate);
    }
    result.push_back(std::move(best_repair));
  }
  if(result.empty())result.push_back(schedule(base));
  int value=-INF;for(const auto& s:result)value=std::max(value,s.lateness);
  if(value<best_value){best_value=value;best_result=std::move(result);}
  }
  return best_result;
}
}
