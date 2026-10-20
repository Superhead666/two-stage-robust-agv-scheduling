// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"
#include <clocale>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

using namespace baseline;

namespace {
struct Args {
  std::map<std::string,std::string> values;
  Args(int argc,char** argv) {
    for(int i=1;i<argc;i+=2) {
      if(i+1>=argc||std::string(argv[i]).rfind("--",0)!=0)
        throw std::runtime_error("arguments must be --key value pairs");
      if(!values.emplace(argv[i],argv[i+1]).second)
        throw std::runtime_error("duplicate argument "+std::string(argv[i]));
    }
  }
  void allow(std::initializer_list<std::string> allowed) const {
    const std::set<std::string> keys(allowed);
    for(const auto& item:values)if(!keys.count(item.first))throw std::runtime_error("unknown argument "+item.first);
  }
  std::string get(const std::string& key,const std::string& fallback="") const {
    auto it=values.find(key);return it==values.end()?fallback:it->second;
  }
  std::string require(const std::string& key) const {
    auto it=values.find(key);if(it==values.end()||it->second.empty())throw std::runtime_error("missing argument "+key);return it->second;
  }
  double real(const std::string& key,double fallback) const {
    auto text=get(key,std::to_string(fallback));std::size_t used=0;double value=std::stod(text,&used);
    if(used!=text.size()||!std::isfinite(value))throw std::runtime_error("invalid number for "+key);return value;
  }
  int integer(const std::string& key,int fallback) const {
    double value=real(key,fallback);if(value!=std::floor(value)||std::abs(value)>2000000000.0)
      throw std::runtime_error("invalid integer for "+key);return static_cast<int>(value);
  }
};

std::string fingerprint(const fs::path& path) {
  std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("cannot read "+path.string());
  std::uint64_t hash=14695981039346656037ULL;char c;
  while(in.get(c)){hash^=static_cast<unsigned char>(c);hash*=1099511628211ULL;}
  std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;return out.str();
}
}

int main(int argc,char** argv) {try {
  std::setlocale(LC_ALL,"C");Args args(argc,argv);
  args.allow({"--load","--out","--mode","--formulation","--limit","--seed","--warm-solution","--initial-lower"});
  const fs::path input=args.require("--load"),prefix=args.require("--out");
  const std::string mode=args.get("--mode","basic");
  if(mode!="basic"&&mode!="hig")throw std::runtime_error("mode must be basic or hig");
  Options opt;opt.formulation=args.get("--formulation","block");opt.limit=args.real("--limit",1200);
  const int seed=args.integer("--seed",1);if(seed<0)throw std::runtime_error("seed must be nonnegative");
  opt.seed=static_cast<unsigned>(seed);opt.threads=1;opt.initial_lower=args.integer("--initial-lower",0);
  if(mode=="hig")opt.warm_solution=args.require("--warm-solution");
  else if(args.values.count("--warm-solution")||opt.initial_lower!=0)
    throw std::runtime_error("basic CCG uses the singleton incumbent and lower bound zero");
  if(!prefix.parent_path().empty())fs::create_directories(prefix.parent_path());
  opt.trace=prefix.string()+"_trace.csv";opt.solution=prefix.string()+".sol";
  auto a=load_instance(input.string());
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.set(GRB_IntParam_Threads,1);
  env.set(GRB_IntParam_Seed,seed);env.start();set_solver_seed(seed);set_oracle(&env);reset_oracle_stats();
  Clock wall;auto result=solve(a,env,opt);const double wall_seconds=wall.elapsed();const auto oracle=oracle_stats();
  std::ofstream out(prefix.string()+".csv");if(!out)throw std::runtime_error("cannot write result CSV");
  out<<"instance,method,formulation,seed,limit,status,reason,lower,upper,absolute_gap,relative_gap,initial_upper,seconds,wall_seconds,time_best,iterations,scenarios,master_seconds,build_seconds,verification_seconds,master_nodes,master_variables,master_constraints,verification_calls,adversary_calls,schedules,deterministic_mip_calls,deterministic_mip_seconds,certificate,input_fingerprint,threads\n";
  const int gap=result.upper<INF?result.upper-result.lower:INF;
  out<<std::setprecision(12)<<input.stem().string()<<','<<mode<<"-ccg,"<<opt.formulation<<','<<opt.seed<<','<<opt.limit<<','<<result.status<<','<<result.reason<<','<<result.lower<<',';
  if(result.upper<INF)out<<result.upper<<','<<gap<<','<<double(gap)/std::max(1,result.upper)<<','<<result.initial;
  else out<<",,,";
  out<<','<<result.seconds<<','<<wall_seconds<<','<<result.time_best<<','<<result.iterations<<','<<result.scenarios<<','<<result.master_seconds<<','<<result.build_seconds<<','<<result.verify_seconds<<','<<result.master_nodes<<','<<result.master_variables<<','<<result.master_constraints<<','<<result.verify_calls<<','<<result.adversary_calls<<','<<result.schedules<<','<<oracle.mip_calls<<','<<oracle.mip_seconds<<','<<result.certificate<<','<<fingerprint(input)<<",1\n";
  if(!out)throw std::runtime_error("result CSV write failed");
  std::cout<<mode<<"-ccg ["<<result.lower<<','<<result.upper<<"] "<<result.status<<" seconds="<<result.seconds<<'\n';
  return 0;
} catch(const GRBException& e) {std::cerr<<"Gurobi error "<<e.getErrorCode()<<": "<<e.getMessage()<<'\n';return 2;}
  catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
