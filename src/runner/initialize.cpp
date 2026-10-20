// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "core.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace agv;
namespace fs=std::filesystem;

int main(int argc,char** argv) {try {
  if(argc!=3)throw std::runtime_error("agv_initialize INSTANCE SOLUTION");
  const auto a=load_instance(argv[1]);
  std::vector<int> order(a.n),ready(a.m,0),count(a.m,0),assignment(a.n,-1);
  std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  int upper=0;
  for(int j:order) {
    int selected=-1,completion=INF;
    for(int k=0;k<a.m;++k) {
      const int candidate=ready[k]+(count[k]?a.alpha:0)+a.p[j][k];
      if(candidate<completion){completion=candidate;selected=k;}
    }
    if(selected<0)throw std::runtime_error("singleton initialization failed");
    assignment[j]=selected;ready[selected]=completion;++count[selected];
    upper=std::max(upper,completion-a.d[j]);
  }
  upper=std::max(0,upper);fs::path path=argv[2];if(!path.parent_path().empty())fs::create_directories(path.parent_path());
  fs::path temporary=path;temporary+=".tmp";{
    std::ofstream out(temporary);if(!out)throw std::runtime_error("cannot write initial solution");
    out<<a.n<<' '<<upper<<" 0\n";for(int k:assignment)out<<k<<' ';
    out<<"\nobjective maximum_tardiness\ncertificate singleton_jobs_with_recharge\n";
  }
  fs::rename(temporary,path);std::cout<<"singleton upper="<<upper<<" lower=0\n";return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
