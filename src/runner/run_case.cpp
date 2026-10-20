// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;

namespace {
struct Args {
  std::map<std::string,std::string> values;
  Args(int argc,char** argv) {for(int i=1;i<argc;i+=2) {
    if(i+1>=argc||std::string(argv[i]).rfind("--",0)!=0)throw std::runtime_error("arguments must be --key value pairs");
    if(!values.emplace(argv[i],argv[i+1]).second)throw std::runtime_error("duplicate argument "+std::string(argv[i]));
  }}
  void allow(std::initializer_list<std::string> allowed) const {std::set<std::string> keys(allowed);for(const auto& x:values)if(!keys.count(x.first))throw std::runtime_error("unknown argument "+x.first);}
  std::string get(const std::string& key,const std::string& fallback="") const {auto it=values.find(key);return it==values.end()?fallback:it->second;}
  std::string require(const std::string& key) const {auto it=values.find(key);if(it==values.end()||it->second.empty())throw std::runtime_error("missing argument "+key);return it->second;}
  double real(const std::string& key,double fallback) const {auto s=get(key,std::to_string(fallback));std::size_t used=0;double x=std::stod(s,&used);if(used!=s.size()||!std::isfinite(x))throw std::runtime_error("invalid number for "+key);return x;}
  int integer(const std::string& key,int fallback) const {double x=real(key,fallback);if(x!=std::floor(x)||std::abs(x)>2000000000.0)throw std::runtime_error("invalid integer for "+key);return static_cast<int>(x);}
};

std::vector<std::string> split(const std::string& line) {std::vector<std::string> fields;std::string field;std::istringstream in(line);while(std::getline(in,field,','))fields.push_back(field);if(!line.empty()&&line.back()==',')fields.emplace_back();return fields;}
std::map<std::string,std::string> read_row(const fs::path& path) {std::ifstream in(path);std::string header,row;if(!std::getline(in,header)||!std::getline(in,row))throw std::runtime_error("cannot read "+path.string());auto names=split(header),values=split(row);if(names.size()!=values.size())throw std::runtime_error("CSV column mismatch in "+path.string());std::map<std::string,std::string> out;for(std::size_t i=0;i<names.size();++i)out[names[i]]=values[i];return out;}
std::string required(const std::map<std::string,std::string>& row,const std::string& key) {auto it=row.find(key);if(it==row.end()||it->second.empty())throw std::runtime_error("missing CSV value "+key);return it->second;}

struct Solution {int n=0,upper=0,lower=0;std::vector<int> assignment;};
Solution read_solution(const fs::path& path,int machines) {std::ifstream in(path);Solution s;if(!(in>>s.n>>s.upper>>s.lower)||s.n<1||s.upper<0||s.lower<0||s.lower>s.upper)throw std::runtime_error("invalid certified solution "+path.string());s.assignment.resize(s.n);for(int& k:s.assignment)if(!(in>>k)||k<0||k>=machines)throw std::runtime_error("invalid assignment in "+path.string());return s;}

int run_process(const std::vector<std::string>& command,const fs::path& log) {
  if(command.empty())throw std::runtime_error("empty child command");const pid_t child=fork();if(child<0)throw std::runtime_error("fork failed");
  if(child==0){int fd=open(log.c_str(),O_WRONLY|O_CREAT|O_TRUNC,0644);if(fd<0)_exit(126);dup2(fd,STDOUT_FILENO);dup2(fd,STDERR_FILENO);close(fd);setenv("LC_ALL","C",1);std::vector<char*> args;for(const auto& x:command)args.push_back(const_cast<char*>(x.c_str()));args.push_back(nullptr);execv(args.front(),args.data());_exit(127);}
  int status=0;while(waitpid(child,&status,0)<0)if(errno!=EINTR)throw std::runtime_error("waitpid failed");
  if(WIFEXITED(status))return WEXITSTATUS(status);if(WIFSIGNALED(status))return 128+WTERMSIG(status);return 125;
}
double elapsed(const Clock::time_point& start){return std::chrono::duration<double>(Clock::now()-start).count();}
void atomic_write(const fs::path& path,const std::string& text){fs::path temporary=path;temporary+=".tmp";{std::ofstream out(temporary);if(!out)throw std::runtime_error("cannot write "+path.string());out<<text;if(!out)throw std::runtime_error("write failed");}fs::rename(temporary,path);}
}

int main(int argc,char** argv) {try {
  Args args(argc,argv);args.allow({"--instance","--output","--method","--total-limit","--hig-limit","--seed","--hig-binary","--lbbd-binary","--ccg-binary","--initialize-binary"});
  const fs::path instance=fs::absolute(args.require("--instance")),output=fs::absolute(args.require("--output"));
  const std::string method=args.require("--method");
  if(method!="lbbd-basic"&&method!="lbbd-hig"&&method!="ccg-basic"&&method!="ccg-hig"&&method!="cp-basic"&&method!="cp-hig")throw std::runtime_error("unknown method "+method);
  const bool use_hig=method.size()>=3&&method.substr(method.size()-3)=="hig";
  const bool use_lbbd=method.rfind("lbbd-",0)==0;
  const bool use_cp=method.rfind("cp-",0)==0;
  const double total_limit=args.real("--total-limit",1800),hig_limit=args.real("--hig-limit",60);
  const int seed=args.integer("--seed",1);if(total_limit<=0||hig_limit<=0||hig_limit>=total_limit||seed<0)throw std::runtime_error("invalid limits or seed");
  if(fs::exists(output/"result.csv"))throw std::runtime_error("result already exists in "+output.string());fs::create_directories(output);
  int n=0,m=0,alpha=0;double gamma=0;{std::ifstream in(instance);if(!(in>>n>>m>>alpha>>gamma)||n<1||m<1)throw std::runtime_error("invalid instance header");}
  const auto start=Clock::now();double preparation_seconds=0;Solution initial;
  const fs::path warm=output/"initial.sol";
  if(use_hig) {
    const fs::path hig=fs::absolute(args.require("--hig-binary"));
    std::vector<std::string> command={hig.string(),"--load",instance.string(),"--limit",std::to_string(hig_limit),"--seed",std::to_string(seed),"--variant","clean","--trace",(output/"hig_trace.csv").string(),"--solution",warm.string(),"--out",(output/"hig.csv").string()};
    if(int code=run_process(command,output/"hig.log"))throw std::runtime_error("HIG exited with code "+std::to_string(code));
    preparation_seconds=elapsed(start);
    initial=read_solution(warm,m);auto row=read_row(output/"hig.csv");
    if(initial.upper!=std::stoi(required(row,"upper"))||initial.lower!=std::stoi(required(row,"global_lower")))throw std::runtime_error("HIG outputs disagree");
  } else if(use_lbbd) {
    const fs::path initializer=fs::absolute(args.require("--initialize-binary"));
    if(int code=run_process({initializer.string(),instance.string(),warm.string()},output/"initialize.log"))throw std::runtime_error("initializer exited with code "+std::to_string(code));
    initial=read_solution(warm,m);preparation_seconds=elapsed(start);
  }
  const double remaining=total_limit-elapsed(start);if(remaining<=0)throw std::runtime_error("preparation exhausted the total time limit");
  int lower=0,upper=0;double solver_seconds=0;std::string status;
  std::map<std::string,std::string> exact_row;
  if(use_lbbd) {
    const fs::path solver=fs::absolute(args.require("--lbbd-binary")),prefix=output/"lbbd";
    std::vector<std::string> command={solver.string(),instance.string(),std::to_string(remaining),prefix.string(),warm.string(),std::to_string(initial.lower),std::to_string(seed)};
    if(int code=run_process(command,output/"lbbd.log"))throw std::runtime_error("LBBD exited with code "+std::to_string(code));
    exact_row=read_row(prefix.string()+".csv");lower=std::stoi(required(exact_row,"lower"));upper=std::stoi(required(exact_row,"upper"));solver_seconds=std::stod(required(exact_row,"seconds"));status=required(exact_row,"reason");
  } else {
    const fs::path solver=fs::absolute(args.require("--ccg-binary")),prefix=output/"ccg";
    std::vector<std::string> command={solver.string(),"--load",instance.string(),"--out",prefix.string(),"--mode",use_hig?"hig":"basic","--formulation",use_cp?"cp":"block","--limit",std::to_string(remaining),"--seed",std::to_string(seed)};
    if(use_hig){command.push_back("--warm-solution");command.push_back(warm.string());command.push_back("--initial-lower");command.push_back(std::to_string(initial.lower));}
    if(int code=run_process(command,output/"ccg.log"))throw std::runtime_error("CCG exited with code "+std::to_string(code));
    exact_row=read_row(prefix.string()+".csv");lower=std::stoi(required(exact_row,"lower"));upper=std::stoi(required(exact_row,"upper"));solver_seconds=std::stod(required(exact_row,"seconds"));status=required(exact_row,"status");
  }
  if(lower>upper)throw std::runtime_error("solver returned contradictory bounds");
  if(use_hig&&(lower<initial.lower||upper>initial.upper))throw std::runtime_error("final interval does not refine HIG certificate");
  const double wall=elapsed(start);const int gap=upper-lower;
  int initial_gap=-1,upper_improvement=-1,lower_improvement=-1,gap_reduction=-1;
  int hig_found_final_incumbent=0,hig_already_certified=0;
  int exact_proved_hig_incumbent=0,exact_no_final_bound_improvement=0;
  std::string progress="not_applicable";
  if(use_hig) {
    initial_gap=initial.upper-initial.lower;
    upper_improvement=initial.upper-upper;
    lower_improvement=lower-initial.lower;
    gap_reduction=initial_gap-gap;
    hig_found_final_incumbent=upper_improvement==0;
    if(initial_gap==0) {hig_already_certified=1;progress="hig_already_certified";}
    else if(upper_improvement>0) progress="incumbent_improved";
    else if(lower==upper) {exact_proved_hig_incumbent=1;progress="hig_incumbent_proved";}
    else if(lower_improvement>0) progress="lower_bound_improved_only";
    else {exact_no_final_bound_improvement=1;progress="no_final_bound_improvement";}
  }
  auto optional=[&](const std::string& key){auto it=exact_row.find(key);return it==exact_row.end()?std::string():it->second;};
  const std::string master_seconds=optional("master_seconds");
  const std::string verification_seconds=use_lbbd?optional("verify_seconds"):optional("verification_seconds");
  const std::string master_nodes=use_lbbd?optional("nodes"):optional("master_nodes");
  const std::string verification_calls=use_lbbd?optional("robust_calls"):optional("verification_calls");
  const std::string deterministic_calls=use_lbbd?optional("det_calls"):optional("deterministic_mip_calls");
  const std::string root_lower=optional("root_lower");
  std::string root_bound_improvement;
  if(use_lbbd&&!root_lower.empty())root_bound_improvement=std::to_string(std::stoi(root_lower)-initial.lower);
  std::string conflict_reduction;
  if(use_lbbd&&!optional("conflict_jobs_before").empty()&&!optional("conflict_jobs_after").empty())
    conflict_reduction=std::to_string(std::stoi(optional("conflict_jobs_before"))-std::stoi(optional("conflict_jobs_after")));
  std::ostringstream csv;csv<<"instance,method,n,m,gamma,alpha,seed,total_limit,hig_limit,initial_lower,initial_upper,preparation_seconds,final_lower,final_upper,optimal,absolute_gap,relative_gap,status,solver_seconds,wall_seconds,initial_gap,upper_improvement,lower_improvement,gap_reduction,hig_found_final_incumbent,hig_already_certified,exact_proved_hig_incumbent,exact_no_final_bound_improvement,exact_progress_class,architecture,lbbd_iterations,root_lower,root_bound_improvement,master_seconds,verification_seconds,master_nodes,candidate_assignments,certified_assignments,conflict_cuts,conflict_jobs_before,conflict_jobs_after,conflict_jobs_removed,unresolved_verifications,verification_calls,deterministic_calls,dff_cuts,dff_seconds,lifted_conflicts,lifted_coefficients,lift_gain,ccg_iterations,ccg_scenarios,ccg_schedules,build_seconds,adversary_calls\n";
  csv<<std::setprecision(12)<<instance.stem().string()<<','<<method<<','<<n<<','<<m<<','<<gamma<<','<<alpha<<','<<seed<<','<<total_limit<<','<<(use_hig?hig_limit:0)<<','<<(use_hig||use_lbbd?initial.lower:0)<<','<<(use_hig||use_lbbd?initial.upper:-1)<<','<<preparation_seconds<<','<<lower<<','<<upper<<','<<(lower==upper)<<','<<gap<<','<<double(gap)/std::max(1,upper)<<','<<status<<','<<solver_seconds<<','<<wall<<',';
  if(use_hig) csv<<initial_gap<<','<<upper_improvement<<','<<lower_improvement<<','<<gap_reduction;
  else csv<<",,,";
  csv<<','<<hig_found_final_incumbent<<','<<hig_already_certified<<','<<exact_proved_hig_incumbent<<','<<exact_no_final_bound_improvement<<','<<progress<<','
    <<optional("architecture")<<','<<(use_lbbd?optional("iterations"):std::string())<<','<<root_lower<<','<<root_bound_improvement<<','
    <<master_seconds<<','<<verification_seconds<<','<<master_nodes<<','<<optional("candidates")<<','<<optional("certified")<<','<<optional("conflicts")<<','
    <<optional("conflict_jobs_before")<<','<<optional("conflict_jobs_after")<<','<<conflict_reduction<<','<<optional("unknown")<<','
    <<verification_calls<<','<<deterministic_calls<<','<<optional("packing_cuts")<<','<<optional("packing_seconds")<<','<<optional("lifted_conflicts")<<','<<optional("lifted_coefficients")<<','<<optional("lift_gain")<<','
    <<(use_lbbd?std::string():optional("iterations"))<<','<<optional("scenarios")<<','<<optional("schedules")<<','<<optional("build_seconds")<<','<<optional("adversary_calls")<<'\n';
  atomic_write(output/"result.csv",csv.str());std::cout<<method<<" ["<<lower<<','<<upper<<"] "<<status<<" wall="<<wall<<'\n';return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
