// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "model_comparison.hpp"
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {try {
  if(argc!=9)throw std::runtime_error(
    "agv_model_compare static-robust|nominal-assignment INSTANCE TIME OUTPUT_PREFIX SEED NOMINAL_SOL ROBUST_SOL ROBUST_LB");
  const std::string mode=argv[1],prefix=argv[4];
  if(mode!="static-robust"&&mode!="nominal-assignment")throw std::runtime_error("unknown model");
  const auto a=agv::load_instance(argv[2]);const double seconds=std::stod(argv[3]);
  const int seed=std::stoi(argv[5]),bound=std::stoi(argv[8]);
  if(!std::isfinite(seconds)||seconds<=0||seed<0||bound<0)throw std::runtime_error("invalid limits, seed, or bound");
  const auto nominal=model_comparison::read_solution(a,argv[6]);
  const auto robust=model_comparison::read_solution(a,argv[7]);
  if(bound!=std::max(robust.lower,nominal.lower))throw std::runtime_error("reference lower bound mismatch");
  agv::Clock clock(seconds);GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.start();
  agv::set_solver_seed(seed);renewed::set_evaluator(&env,1);
  auto s=mode=="static-robust"?
    model_comparison::solve_static(a,env,clock,seed,{nominal.assignment,robust.assignment},bound):
    model_comparison::evaluate_assignment(a,nominal.assignment,env,clock,bound);
  std::ofstream sol(prefix+".sol");if(!sol)throw std::runtime_error("cannot write solution");
  sol<<a.n<<' '<<s.upper<<' '<<s.lower<<'\n';for(int k:s.assignment)sol<<k<<' ';sol<<'\n';
  std::ofstream table(prefix+".csv");if(!table)throw std::runtime_error("cannot write result");
  table<<"model,lower,upper,optimal,absolute_gap,relative_gap,status,seconds,nodes,verification_calls,deterministic_calls,adversary_calls\n";
  table<<std::setprecision(12)<<mode<<','<<s.lower<<','<<s.upper<<','<<(s.lower==s.upper)<<','
    <<s.upper-s.lower<<','<<(s.upper?double(s.upper-s.lower)/s.upper:0)<<','<<s.status<<','
    <<s.seconds<<','<<s.nodes<<','<<s.verification_calls<<','<<s.deterministic_calls<<','<<s.adversary_calls<<'\n';
  if(mode=="static-robust") {
    std::ofstream schedule(prefix+"_schedule.csv");
    schedule<<"agv,position,job,recharge_before,completion,block,worst_block_energy\n";
    for(int k=0;k<a.m;++k) {
      agv::Local local=a.local(k,s.assignment);const auto& policy=s.schedules[k];int t=0,b=0;
      for(int q=0;q<local.n();++q) {
        if(policy.recharge[q]){++b;t+=a.alpha;}const int j=policy.order[q];t+=local.p[j];
        schedule<<std::setprecision(17)<<k<<','<<q<<','<<local.ids[j]<<','<<policy.recharge[q]<<','
          <<t<<','<<b<<','<<worst_block_energy(local,policy.blocks[b])<<'\n';
      }
    }
  }
  std::cout<<mode<<" lower="<<s.lower<<" upper="<<s.upper<<" status="<<s.status<<" seconds="<<s.seconds<<'\n';
}catch(const GRBException& e){std::cerr<<e.getMessage()<<'\n';return 1;}
 catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
