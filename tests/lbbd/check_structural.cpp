// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#define PATH_BOUND_LIBRARY
#include "path_bound.cpp"
Local subset(const Local& a,Mask mask) {
  Local b;b.capacity=a.capacity;b.gamma=a.gamma;b.alpha=a.alpha;
  for(int j=0;j<a.n();++j)if((mask>>j)&1){b.p.push_back(a.p[j]);b.d.push_back(a.d[j]);b.ids.push_back(a.ids[j]);b.e.push_back(a.e[j]);b.nominal.push_back(a.nominal[j]);b.deviation.push_back(a.deviation[j]);}return b;
}
int main(){try {
  std::mt19937 rng(48711);long comparisons=0,walks=0,lift_checks=0;
  for(int c=0;c<30;++c) {
    Instance a;a.n=4;a.m=1;a.alpha=2;a.gamma=0;a.capacity={7};std::vector<double> weights;
    for(int j=0;j<4;++j){a.p.push_back({1+int(rng()%3)});a.d.push_back(2+rng()%7);a.nominal.push_back({double(1+rng()%4)});a.deviation.push_back({0});weights.push_back((1+rng()%9)/10.);}
    for(Mask tracked:{Mask(0),Mask(5),Mask(15)}) {
      double best=0;
      std::function<void(int,int,Mask,double)> enumerate=[&](int t,int energy,Mask used,double value) {
        ++walks;best=std::max(best,value);
        for(int j=0;j<4;++j)if(!((tracked&used)&(Mask(1)<<j))) {
          int e=int(a.nominal[j][0]);
          if(energy+e<=7&&t+a.p[j][0]<=a.d[j])enumerate(t+a.p[j][0],energy+e,used|(Mask(1)<<j),value+weights[j]);
          if(energy>0&&t+2+a.p[j][0]<=a.d[j])enumerate(t+2+a.p[j][0],e,used|(Mask(1)<<j),value+weights[j]);
        }
      };
      enumerate(0,0,0,0);Clock clock;auto result=price(a,0,0,weights,clock,tracked);
      if(!result.complete||std::abs(result.profit-best)>1e-8)throw std::runtime_error("path pricing disagrees with independent walk enumeration");++comparisons;
    }
  }
  for(int c=0;c<60;++c) {
    auto a=generate(5,1,448000+c,1.5,.5,2,.4);std::vector<int> x(5,0);Local full=a.local(0,x);
    const int value=std::max(0,brute_optimum(full));int base=0;
    for(int j=0;j<5;++j)base=std::max(base,full.p[j]-full.d[j]);
    for(Mask mask=0;mask<32;++mask) {
      auto remaining=subset(full,mask);const int actual=std::max(base,remaining.n()?brute_optimum(remaining):0);
      int rhs=value;for(int j=0;j<5;++j)if(!((mask>>j)&1))rhs-=std::min(value-base,full.p[j]+full.alpha);
      if(rhs>actual)throw std::runtime_error("invalid insertion lifting");++lift_checks;
    }
  }
  std::cout<<"PASS "<<comparisons<<" DP comparisons against "<<walks<<" independently enumerated walks and "<<lift_checks<<" insertion-cut subset checks\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
