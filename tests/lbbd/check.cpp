// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "packing.hpp"
#include <filesystem>
#include <iostream>
using namespace lb;
int main(int argc,char** argv){try {
  long checks=0;std::mt19937 rng(31771);
  for(int c=0;c<60;++c) {
    auto a=generate(5,1,318000+c,.5+.5*(c%5),.5,2,.3);
    // Coarse rational energies make the independent grid checks well separated
    // from the deliberately conservative positive-margin convention.
    a.capacity[0]=10;
    for(int j=0;j<a.n;++j){a.nominal[j][0]=(10+rng()%51)/10.;a.deviation[j][0]=(rng()%21)/10.;}
    std::vector<int> jobs={0,1,2,3,4};
    for(int q:{1,2,3,4,5,6,7,8,16,32,59,97}) {
      const auto result=promote(a,0,jobs,q);
      const double budget=std::accumulate(result.xi.begin(),result.xi.end(),0.0);
      if(budget>a.gamma+1e-9)throw std::runtime_error("budget violation");
      double grid=0;
      for(int code=0;code<3125;++code) {
        int v=code;double used=0,value=0;
        for(int j=0;j<5;++j){const double xi=(v%5)/4.;v/=5;used+=xi;value+=transformed_size((a.nominal[j][0]+a.deviation[j][0]*xi)/10,q);}
        if(used<=a.gamma+1e-10)grid=std::max(grid,value);
      }
      if(result.value+1e-8<grid)throw std::runtime_error("DP below independent scenario grid");
      for(int mask=0;mask<32;++mask) {
        double energy=0,transformed=0;
        for(int j=0;j<5;++j)if((mask>>j)&1) {
          double e=a.nominal[j][0]+a.deviation[j][0]*result.xi[j];energy+=e;transformed+=transformed_size(e/10,q);
        }
        if(energy<=10+1e-10&&transformed>1+1e-8)throw std::runtime_error("invalid packing witness weights");
        ++checks;
      }
    }
  }
  if(argc>1) {
    std::filesystem::create_directories(argv[1]);
    for(int n:{4,10}) {
      Instance a;a.n=n;a.m=1;a.alpha=10;a.gamma=n==4?.8:2.5;a.capacity={10};
      a.p.assign(n,{1});a.d.assign(n,0);a.nominal.assign(n,{n==4?4.8:4.9});a.deviation.assign(n,{1});
      std::vector<int> jobs(n);std::iota(jobs.begin(),jobs.end(),0);Clock clock;
      auto b=packing_bound(a,0,jobs,clock);
      const int expected=n==4?3:10;
      if(int(std::ceil(b.value-1e-8))!=expected)throw std::runtime_error("fractional promotion example");
      save_instance(a,std::string(argv[1])+"/promotion"+std::to_string(n)+".agv");
    }
  }
  std::cout<<"PASS "<<checks<<" subset checks, 720 independent grid comparisons, and two strict-strengthening examples\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
