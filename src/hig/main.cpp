// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "long_search.hpp"
#include "evaluator.hpp"
#include <clocale>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

using namespace agv;
namespace fs = std::filesystem;

#ifndef AGV_HIG_ROOT_MODE
#define AGV_HIG_ROOT_MODE 3
#endif

namespace {
struct Args {
  std::map<std::string,std::string> values;
  Args(int argc,char** argv) {
    for(int i=1;i<argc;i+=2) {
      if(i+1>=argc || std::string(argv[i]).rfind("--",0)!=0)
        throw std::runtime_error("arguments must be --key value pairs");
      if(!values.emplace(argv[i],argv[i+1]).second)
        throw std::runtime_error("duplicate argument "+std::string(argv[i]));
    }
  }
  void allow(std::initializer_list<std::string> allowed) const {
    const std::set<std::string> keys(allowed);
    for(const auto& item:values) if(!keys.count(item.first))
      throw std::runtime_error("unknown argument "+item.first);
  }
  std::string get(const std::string& key,const std::string& fallback="") const {
    const auto it=values.find(key);return it==values.end()?fallback:it->second;
  }
  std::string require(const std::string& key) const {
    const auto it=values.find(key);if(it==values.end()||it->second.empty())
      throw std::runtime_error("missing argument "+key);return it->second;
  }
  double real(const std::string& key,double fallback) const {
    const auto text=get(key,std::to_string(fallback));std::size_t used=0;
    const double value=std::stod(text,&used);
    if(used!=text.size()||!std::isfinite(value))throw std::runtime_error("invalid number for "+key);
    return value;
  }
  int integer(const std::string& key,int fallback) const {
    const double value=real(key,fallback);
    if(value!=std::floor(value)||std::abs(value)>2000000000.0)
      throw std::runtime_error("invalid integer for "+key);
    return static_cast<int>(value);
  }
};

void atomic_write(const fs::path& path,const std::string& text) {
  if(!path.parent_path().empty())fs::create_directories(path.parent_path());
  fs::path temporary=path;temporary+=".tmp";
  {std::ofstream out(temporary);if(!out)throw std::runtime_error("cannot write "+path.string());out<<text;if(!out)throw std::runtime_error("write failed for "+path.string());}
  fs::rename(temporary,path);
}
}

int main(int argc,char** argv) {try {
  std::setlocale(LC_ALL,"C");
  Args args(argc,argv);
  if(args.values.count("--generate")) {
    args.allow({"--generate","--n","--m","--seed","--gamma","--variation","--capacity","--due",
      "--tardiness-factor","--due-date-range","--processing-max","--nominal-energy-max",
      "--recharge-factor","--recharge-time"});
    if(args.values.count("--recharge-factor")&&args.values.count("--recharge-time"))
      throw std::runtime_error("specify either --recharge-factor or --recharge-time");
    if(args.values.count("--due")&&args.values.count("--tardiness-factor"))
      throw std::runtime_error("specify either --due or --tardiness-factor");
    const int n=args.integer("--n",0),m=args.integer("--m",0),seed=args.integer("--seed",1);
    if(n<1||n>62||m<1||m>32||seed<0)throw std::runtime_error("invalid n, m, or seed");
    auto data=longrun::generate_paper(n,m,static_cast<unsigned>(seed),args.real("--gamma",2.5),
      args.real("--variation",.3),args.real("--capacity",4),args.real("--due",.45),
      args.integer("--processing-max",20),args.integer("--nominal-energy-max",20),
      args.real("--recharge-factor",2),args.integer("--recharge-time",-1),
      args.real("--tardiness-factor",-1),args.real("--due-date-range",.5));
    save_instance(data,args.require("--generate"));return 0;
  }
  args.allow({"--load","--limit","--seed","--variant","--trace","--solution","--out"});
  if(args.get("--variant","clean")!="clean")
    throw std::runtime_error("the submission executable supports only --variant clean");
  const double limit=args.real("--limit",30);const int seed=args.integer("--seed",1);
  if(limit<=0||seed<0)throw std::runtime_error("limit must be positive and seed nonnegative");
  const fs::path trace_path=args.require("--trace"),solution_path=args.require("--solution"),csv_path=args.require("--out");
  if(!trace_path.parent_path().empty())fs::create_directories(trace_path.parent_path());
  auto a=longrun::read_instance(args.require("--load"));
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.set(GRB_IntParam_Threads,1);
  env.set(GRB_IntParam_Seed,seed);env.start();
  agv::set_solver_seed(seed);renewed::set_evaluator(&env,1);renewed::reset_evaluator_stats();

  auto o=longrun::options_for("manuscript");o.limit=limit;o.seed=static_cast<unsigned>(seed);
  o.root_mode=AGV_HIG_ROOT_MODE;

  std::ofstream trace(trace_path);if(!trace)throw std::runtime_error("cannot write "+trace_path.string());
  trace<<"seconds,upper,lower,verifications,unknown,assignment\n";
  int previous=INF;double last=-1;
  o.progress=[&](double t,int upper,int lower,const longrun::Stats& s,const std::vector<int>& x) {
    upper=std::max(0,upper);lower=std::max(0,lower);
    if(upper>previous)throw std::runtime_error("increasing HIG incumbent");
    if(upper<previous||t-last>=1||t>=o.limit) {
      trace<<std::setprecision(12)<<t<<','<<upper<<','<<lower<<','<<s.verifications<<','<<s.unknown<<',';
      for(int k:x)trace<<k<<';';trace<<'\n';trace.flush();last=t;previous=upper;
    }
  };
  longrun::Search search(a,env,o);auto r=search.run();const auto& s=r.stats;
  const int upper=std::max(0,r.upper),assignment_lower=std::max(0,r.assignment_lower);
  const int global_lower=std::max(0,r.lower),static_upper=std::max(0,r.static_upper);
  if(global_lower>upper||assignment_lower>upper||r.assignment.size()!=static_cast<std::size_t>(a.n))
    throw std::runtime_error("invalid HIG certificate");
  std::ostringstream sol;sol<<a.n<<' '<<upper<<' '<<global_lower<<'\n';
  for(int k:r.assignment)sol<<k<<' ';sol<<"\nobjective maximum_tardiness\ncertificate HIG\n";
  atomic_write(solution_path,sol.str());
  std::ostringstream out;
  out<<"instance,variant,seed,limit,upper,assignment_lower,global_lower,static_upper,seconds,time_best,entries,cache_hits,evaluations,verifications,unknown,verification_seconds,adversary_seconds,entry_seconds,lower_seconds,screen_seconds,master_seconds,portfolio_entries,portfolio_schedules,portfolio_seconds,priority_jobs,critical_jobs,priority_scans,chain_seconds,witness_checks,witness_prunes,lower_hits,lower_entries\n";
  out<<std::setprecision(12)<<fs::path(a.id).stem().string()<<",clean,"<<o.seed<<','<<o.limit<<','<<upper<<','<<assignment_lower<<','<<global_lower<<','<<static_upper<<','<<r.seconds<<','<<r.time_best<<','<<s.entries<<','<<s.cache_hits<<','<<s.evaluations<<','<<s.verifications<<','<<s.unknown<<','<<s.verification_seconds<<','<<s.adversary_seconds<<','<<s.entry_seconds<<','<<s.lower_seconds<<','<<s.screen_seconds<<','<<s.mip_seconds<<','<<s.portfolio_entries<<','<<s.portfolio_schedules<<','<<s.portfolio_seconds<<','<<s.prioritized_jobs<<','<<s.critical_jobs<<','<<s.priority_scans<<','<<s.chain_seconds<<','<<s.witness_checks<<','<<s.witness_prunes<<','<<s.lower_hits<<','<<s.lower_entries<<'\n';
  atomic_write(csv_path,out.str());
  std::cout<<"clean "<<fs::path(a.id).stem().string()<<" UB="<<upper<<" fixedLB="<<assignment_lower<<" seconds="<<r.seconds<<'\n';
  return 0;
} catch(const GRBException& e) {std::cerr<<"Gurobi error "<<e.getErrorCode()<<": "<<e.getMessage()<<'\n';return 2;}
  catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
