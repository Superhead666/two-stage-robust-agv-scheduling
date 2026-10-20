// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "robust.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace agv;
struct Args {
  std::unordered_map<std::string,std::string> values;
  Args(int argc,char**argv){for(int q=2;q<argc;q+=2){if(q+1>=argc)throw std::runtime_error("options need values");values[argv[q]]=argv[q+1];}}
  std::string get(const std::string& key,const std::string& fallback="")const{auto it=values.find(key);return it==values.end()?fallback:it->second;}
  int integer(const std::string& key,int fallback)const{return std::stoi(get(key,std::to_string(fallback)));}
  double real(const std::string& key,double fallback)const{return std::stod(get(key,std::to_string(fallback)));}
};
void check(bool b,const std::string& message){if(!b)throw std::runtime_error("VALIDATION FAILED: "+message);}
RobustOptions options(std::string mode){RobustOptions o;
  const std::vector<std::string> accepted={"basic","sparse","enhanced","strong","workload","refine","cuts","reuse","polished","policy","policy-cuts","lp","order","order-threshold","order-logic","lifted","order-lifted","orders","orders-threshold","orders-logic","orders-lifted","orders-tight","orders-tight-threshold","orders-tight-lifted","cumulative","cumulative-orders","primal","threshold","hybrid","logic","no-early","no-share","no-screen","no-bounds"};
  if(std::find(accepted.begin(),accepted.end(),mode)==accepted.end())throw std::runtime_error("unknown algorithm mode: "+mode);
  if(mode=="cumulative")mode="lifted";if(mode=="cumulative-orders")mode="orders-tight-lifted";
  if(mode=="order"||mode=="order-threshold"||mode=="order-logic"||mode=="order-lifted"){mode="reuse";o.fixed_order=true;}
  if(mode=="orders-tight"||mode=="orders-tight-threshold"||mode=="orders-tight-lifted"){mode="orders";o.order_bounds=true;}
  if(mode=="orders"||mode=="orders-threshold"||mode=="orders-logic"||mode=="orders-lifted"){mode="reuse";o.order_generation=true;}
  if(mode=="lifted")mode="reuse";
  if(mode=="lp"){mode="policy";o.branch_coverage=true;}
  if(mode=="primal"||mode=="threshold"||mode=="hybrid"||mode=="logic")mode="policy";
  if(mode=="basic"){o.enhanced=false;o.screening=false;o.shared_blocks=false;o.early_adversary=false;o.det.dense=true;}
  if(mode=="sparse"){o.enhanced=false;o.screening=false;o.shared_blocks=false;o.early_adversary=false;}
  if(mode=="no-early")o.early_adversary=false;
  if(mode=="no-share")o.shared_blocks=false;
  if(mode=="no-screen")o.screening=false;
  if(mode=="no-bounds")o.det.bounds=false;
  if(mode=="strong"||mode=="workload"||mode=="cuts"||mode=="reuse"||mode=="polished"||mode=="policy"||mode=="policy-cuts")o.workload_screen=true;
  if(mode=="strong"||mode=="refine"||mode=="cuts"||mode=="reuse"||mode=="polished"||mode=="policy"||mode=="policy-cuts")o.refine=true;
  if(mode=="cuts"||mode=="polished"||mode=="policy-cuts")o.overload_cuts=true;
  if(mode=="reuse"||mode=="polished"||mode=="policy"||mode=="policy-cuts")o.reuse_schedules=true;
  if(mode=="policy"||mode=="policy-cuts")o.static_policy=true;
  return o;
}
FleetOptions fleet_options(const std::string& mode){FleetOptions opt;opt.enhanced=mode!="basic";opt.robust=options(mode);opt.robust_workload=opt.robust.workload_screen;
  opt.cumulative_workload=mode=="cumulative"||mode=="cumulative-orders";
  opt.insertion_cuts=mode=="lifted"||mode=="order-lifted"||mode=="orders-lifted"||mode=="orders-tight-lifted"||opt.cumulative_workload;
  opt.root_bound=opt.insertion_cuts||mode=="order-threshold"||mode=="order-logic"||mode=="orders-threshold"||mode=="orders-logic"||mode=="orders-tight-threshold";
  opt.primal_search=mode=="primal"||mode=="threshold"||mode=="hybrid"||mode=="logic"||opt.root_bound;
  opt.threshold_master=mode=="threshold"||mode=="hybrid"||mode=="logic"||opt.root_bound;opt.conflict_cuts=mode=="hybrid"||mode=="logic"||mode=="order-logic"||mode=="orders-logic"||opt.insertion_cuts;opt.logic_only=mode=="logic"||mode=="order-logic"||mode=="orders-logic"||opt.insertion_cuts;opt.max_iterations=opt.logic_only?10000:1000;return opt;
}
void emit(const Args& args,const std::string& header,const std::string& row){
  std::cout<<row<<std::endl;std::string path=args.get("--out");if(path.empty())return;
  auto parent=std::filesystem::path(path).parent_path();if(!parent.empty())std::filesystem::create_directories(parent);
  bool fresh=!std::filesystem::exists(path)||std::filesystem::file_size(path)==0;
  std::ofstream f(path,std::ios::app);if(!f)throw std::runtime_error("output file");if(fresh)f<<header<<'\n';f<<row<<'\n';
}
Instance instance(const Args& args,int run){
  if(!args.get("--load").empty())return load_instance(args.get("--load"));
  Instance a=generate(args.integer("--n",12),args.integer("--m",3),args.integer("--seed",1)+run,args.real("--gamma",2),args.real("--variation",.6),args.real("--capacity",3),args.real("--due",1));
  if(!args.get("--save").empty()){std::filesystem::create_directories(args.get("--save"));save_instance(a,args.get("--save")+"/"+a.id+".agv");}return a;
}
void validation(GRBEnv& env,const Args& args){
  Clock total;std::mt19937 rng(8302026);int count=args.integer("--cases",180), comparisons=0;
  for(int c=0;c<count;++c){int n=1+rng()%7;Instance a=generate(n,1,1000+c,std::min(double(n),(c%7)*.5),.6,2,1);Local b=a.local(0,std::vector<int>(n,0));b.alpha=c%9==0?0:a.alpha;
    for(int j=0;j<n;++j){b.e[j]=static_cast<double>(rng()%41)/4;b.capacity=10;b.p[j]=1+rng()%8;b.d[j]=static_cast<int>(rng()%60)-(c%10==0?15:0);}
    int truth=brute_optimum(b);
    auto parts=enumerate_partitions(b);int partbest=INF;
    for(const auto& s:parts)if(valid_schedule(b,s,INF-1))partbest=std::min(partbest,evaluate(b,s.order,s.recharge).lateness);
    check(partbest==truth,"canonical block ordering");
    for(int shift:{-1,0,1})for(int mode=0;mode<4;++mode){DetOptions opt;opt.dense=mode==0;if(mode==2)opt.greedy_recharge=false;if(mode==3){opt.cross_dominance=false;opt.bounds=false;opt.heuristic=false;}
      Clock budget(10);auto r=deterministic(b,truth+shift,opt,budget);check(r.status!=Status::Unknown,"det timeout");check((r.status==Status::Feasible)==(shift>=0),"DP/DFS vs brute");if(r.status==Status::Feasible)check(valid_schedule(b,r.schedule,truth+shift),"witness");++comparisons;}
  }
  std::cout<<"deterministic_checks="<<comparisons<<" canonical_partition_checks="<<count<<" seconds="<<total.elapsed()<<std::endl;
  Local not_edd;not_edd.p={1,1,1,1};not_edd.d={1,2,3,4};not_edd.e={6,6,4,4};not_edd.capacity=10;not_edd.alpha=10;
  check(edd(not_edd).lateness==20 && brute_optimum(not_edd)==11,"global EDD is not sufficient");
  Local trap;trap.p={1,1,1};trap.d={1,2,3};trap.nominal={4,4,4};trap.deviation={4,4,4};trap.e=trap.nominal;trap.capacity=10;trap.gamma=1;trap.alpha=10;
  for(int q=0;q<4;++q){std::vector<double> xi(3,0);if(q)xi[q-1]=1;check(brute_optimum(realize(trap,xi))<=11,"vertex trap vertices");}
  check(brute_optimum(realize(trap,{1.0/3,1.0/3,1.0/3}))==20,"vertex trap interior");
  for(std::string mode:{"basic","enhanced","strong","polished","policy","lp","order","orders","orders-tight"})for(int lambda:{11,19,20}){Clock budget(20);auto r=verify(trap,lambda,env,options(mode),budget);check((r.status==Status::Feasible)==(lambda==20)&&r.status!=Status::Unknown,"continuous adversary trap");}
  Local adaptive=trap;adaptive.d={3,3,3};adaptive.deviation={3,3,3};
  int static_best=INF;for(const auto& schedule:enumerate_partitions(adaptive)){bool safe=true;for(Mask block:schedule.blocks)if(worst_block_energy(adaptive,block)>adaptive.capacity+ENERGY_EPS)safe=false;if(safe)static_best=std::min(static_best,schedule.lateness);}check(static_best==20,"static partition lower bound");
  for(std::string mode:{"basic","enhanced","strong","polished","policy","lp","order","orders","orders-tight"})for(int lambda:{9,10}){Clock budget(20);auto r=verify(adaptive,lambda,env,options(mode),budget);check(r.status==(lambda==10?Status::Feasible:Status::Infeasible),"adaptive recharge advantage");}
  {Clock budget(20);check(verify_fixed_order(adaptive,{0,1,2},10,env,budget).status==Status::Infeasible,"fixed order is not fully adaptive sequencing");}
  Local order_benefit;order_benefit.p={1,1,1,1};order_benefit.d={4,4,4,4};order_benefit.nominal={4,4,4,4};order_benefit.deviation={3,3,3,3};order_benefit.e=order_benefit.nominal;order_benefit.capacity=10;order_benefit.gamma=1;order_benefit.alpha=10;
  for(int lambda:{19,20}){Clock budget(20);check(verify_fixed_order(order_benefit,{0,1,2,3},lambda,env,budget).status==(lambda==20?Status::Feasible:Status::Infeasible),"adaptive charging with fixed order");}
  int robustchecks=0,orderchecks=0;
  for(int c=0;c<args.integer("--robust-cases",30);++c){int n=2+c%5;double gamma=c%7==0?n:(c%5)*.5;Instance a=generate(n,1,5000+c,gamma,c%11==0?0:.7,2,1);Local b=a.local(0,std::vector<int>(n,0));if(c%13==0)b.alpha=0;if(c%17==0)for(int& d:b.d)d+=100;Clock budget(30);int nom=brute_optimum(b),box=brute_optimum(realize(b,std::vector<double>(n,1)));
    std::unordered_map<std::string,std::vector<Schedule>> pools;
    for(int lambda:{box,nom,(nom+box)/2,nom-1,box+1}){auto full=verify_full_partitions(b,lambda,env,budget);check(full.status!=Status::Unknown,"full coverage timeout");
      for(std::string mode:{"basic","enhanced","strong","polished","policy","lp","order","orders","orders-tight"}){auto r=verify(b,lambda,env,options(mode),budget,&pools[mode]);check(r.status==full.status,"generated vs complete coverage: "+mode);++robustchecks;}
      std::vector<int> order(n);std::iota(order.begin(),order.end(),0);for(int repeat=0;repeat<2;++repeat){std::shuffle(order.begin(),order.end(),rng);auto reference=verify_order_enumeration(b,order,lambda,env,budget);for(bool tight:{false,true}){auto r=verify_fixed_order(b,order,lambda,env,budget,tight);check(r.status==reference.status&&r.status!=Status::Unknown,"fixed-order shortest paths vs complete recharge enumeration");++orderchecks;}}}
  }
  std::cout<<"robust_checks="<<robustchecks<<" fixed_order_checks="<<orderchecks<<" vertex_counterexample=passed seconds="<<total.elapsed()<<std::endl;
  int fleetchecks=0;
  for(int c=0;c<args.integer("--fleet-cases",5);++c){int n=4+c%2;double gamma=c%4==0?0:c%4==1?1:c%4==2?1.5:n;Instance a=generate(n,2,9000+c,gamma,c%7==0?0:.6,2,c%3==0?.5:1);if(c%9==0)a.alpha=0;if(c%11==0)for(int& d:a.d)d+=100;Clock budget(60);int truth=enumerate_fleet(a,env,budget);check(truth<INF,"enumeration timeout");
    for(std::string mode:{"basic","enhanced","strong","polished","policy","primal","threshold","hybrid","logic","order-threshold","order-logic","lifted","order-lifted","orders-threshold","orders-lifted","orders-tight-threshold","orders-tight-lifted","cumulative","cumulative-orders"}){auto r=solve_fleet(a,env,fleet_options(mode),budget);check(r.lower==truth&&r.upper==truth,"fleet vs complete enumeration: "+mode);for(int k=0;k<a.m;++k){auto witness=verify_full_partitions(a.local(k,r.assignment),truth,env,budget);check(witness.status==Status::Feasible,"returned assignment certification");}++fleetchecks;}}
  std::cout<<"fleet_checks="<<fleetchecks<<" all_passed=1 total_seconds="<<total.elapsed()<<std::endl;
  int insertionchecks=0;
  for(int c=0;c<50;++c){int n=3+c%3;auto a=generate(n,2,14000+c,1,.6,2,.5);if(c%9==0)a.alpha=0;auto b=realize(a.local(0,std::vector<int>(n,0)),std::vector<double>(n,1.0/n));int value=brute_optimum(b),lower=-INF;for(int j=0;j<n;++j)lower=std::max(lower,*std::min_element(a.p[j].begin(),a.p[j].end())-a.d[j]);if(value<=lower)continue;
    for(Mask mask=0;mask<(Mask(1)<<n);++mask){Local retained;retained.capacity=b.capacity;retained.alpha=b.alpha;int rhs=value;for(int j=0;j<n;++j)if(mask&(Mask(1)<<j)){retained.p.push_back(b.p[j]);retained.d.push_back(b.d[j]);retained.e.push_back(b.e[j]);}else rhs-=std::min(value-lower,b.p[j]+b.alpha);check(std::max(lower,brute_optimum(retained))>=rhs,"insertion inequality over all retained subsets");++insertionchecks;}}
  for(int invalid=0;invalid<6;++invalid){bool rejected=false;try{generate(invalid==0?0:invalid==1?63:4,invalid==2?0:2,1,invalid==3?5:invalid==4?-1:1,invalid==5?-1:.6);}catch(const std::exception&){rejected=true;}check(rejected,"invalid generator input");}
  std::cout<<"insertion_cut_checks="<<insertionchecks<<" invalid_input_checks=6 all_passed=1 total_seconds="<<total.elapsed()<<std::endl;
  std::ostringstream record;record<<comparisons<<','<<count<<','<<robustchecks<<','<<orderchecks<<','<<fleetchecks<<','<<insertionchecks<<','<<total.elapsed()<<",passed";
  emit(args,"deterministic_checks,partition_checks,robust_checks,fixed_order_checks,fleet_checks,insertion_cut_checks,seconds,status",record.str());
}
int main(int argc,char**argv){try{
  // Gurobi's bundled runtime rejects some host locale identifiers. This affects only this process.
  setenv("LC_ALL","C",1);
  if(argc<2){std::cout<<"Commands: validate, dp, local, fleet. Options: --n N --m M --seed S --runs R --mode basic|sparse|enhanced|no-early|no-share|no-screen|no-bounds --limit seconds --gamma G --capacity scale --variation v --due factor --target optimum|nominal|mid|box --out path --save directory --load path\n";return 0;}
  Args args(argc,argv);GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.start();std::string command=argv[1];
  if(command=="validate"){validation(env,args);return 0;}
  std::string mode=args.get("--mode","enhanced");
  for(int run=0;run<args.integer("--runs",1);++run){Instance a=instance(args,run);Clock budget(args.real("--limit",30));std::ostringstream row;row<<std::setprecision(9);
    if(command=="fleet"){
      auto r=solve_fleet(a,env,fleet_options(mode),budget);
      row<<a.id<<','<<mode<<','<<a.n<<','<<a.m<<','<<a.gamma<<','<<r.status<<','<<r.lower<<','<<r.upper<<','<<r.seconds<<','<<r.master_seconds<<','<<r.verify_seconds<<','<<r.iterations<<','<<r.scenarios<<','<<r.verify_calls<<','<<r.adversary_calls<<','<<r.schedules<<','<<r.conflicts;
      emit(args,"instance,mode,n,m,gamma,status,lower,upper,seconds,master_seconds,verify_seconds,iterations,scenarios,verify_calls,adversary_calls,schedules,conflicts",row.str());
    }else{
      std::vector<int> assignment(a.n,0);Local b=a.local(0,assignment); // Local commands intentionally put every job on one AGV.
      if(command=="dp"){
        std::vector<double> xi(a.n,std::min(1.0,a.gamma/a.n));b=realize(b,xi);int lambda=args.integer("--lambda",edd(b).lateness-args.integer("--shift",1));DetOptions opt=options(mode).det;
        auto r=deterministic(b,lambda,opt,budget);
        row<<a.id<<','<<mode<<','<<a.n<<','<<lambda<<','<<status_name(r.status)<<','<<r.seconds<<','<<r.nodes<<','<<r.dominated<<','<<r.pruned;
        emit(args,"instance,mode,n,lambda,status,seconds,nodes,dominated,pruned",row.str());
      }else if(command=="local"){
        RobustOptions opt=options(mode);Clock prep(std::min(10.0,budget.left()));int lo=deterministic_optimum(b,DetOptions{},prep);int hi=edd(realize(b,std::vector<double>(b.n(),std::min(1.0,b.gamma)))).lateness;
        if(lo==INF){int safe_lower=-INF;for(int j=0;j<b.n();++j)safe_lower=std::max(safe_lower,b.p[j]-b.d[j]);row<<a.id<<','<<mode<<','<<a.n<<','<<a.gamma<<','<<args.get("--target","optimum")<<",preparation_limit,"<<safe_lower<<','<<hi<<','<<budget.elapsed()<<",0,0,0,0,0,0,0,nominal_optimization_limit";emit(args,"instance,mode,n,gamma,target,status,lower,upper,seconds,mip_seconds,det_seconds,verify_calls,adversary_calls,schedules,det_calls,det_nodes,reason",row.str());continue;}
        int lower=lo,upper=hi,calls=0,schedules=0,adv=0,det=0;std::uint64_t nodes=0;double mt=0,dt=0;std::string result="optimal_solver_tolerance",reason;
        std::string target=args.get("--target","optimum");
        std::vector<Schedule> pool;
        if(target!="optimum"){int lambda=target=="nominal"?lo:target=="box"?hi:(lo+hi)/2;lo=lambda;hi=lambda+1;}
        while(lo<hi && !budget.expired()){
          int lambda=target=="optimum"?lo+(hi-lo)/2:lo;auto r=verify(b,lambda,env,opt,budget,&pool);++calls;schedules+=r.schedules;adv+=r.adversary_calls;det+=r.det_calls;nodes+=r.det_nodes;mt+=r.mip_seconds;dt+=r.det_seconds;reason=r.reason;
          if(r.status==Status::Unknown){result="unresolved_"+r.reason;break;}
          if(target!="optimum"){result=status_name(r.status);lower=lambda;upper=lambda;lo=hi;break;}
          if(r.status==Status::Feasible)hi=lambda;else lo=lambda+1;lower=lo;upper=hi;
        }
        if(lo<hi && result=="optimal_solver_tolerance")result="time_limit";
        row<<a.id<<','<<mode<<','<<a.n<<','<<a.gamma<<','<<target<<','<<result<<','<<lower<<','<<upper<<','<<budget.elapsed()<<','<<mt<<','<<dt<<','<<calls<<','<<adv<<','<<schedules<<','<<det<<','<<nodes<<','<<reason;
        emit(args,"instance,mode,n,gamma,target,status,lower,upper,seconds,mip_seconds,det_seconds,verify_calls,adversary_calls,schedules,det_calls,det_nodes,reason",row.str());
      }else throw std::runtime_error("unknown command");
    }
  }return 0;
}catch(const GRBException& e){std::cerr<<"Gurobi error "<<e.getErrorCode()<<": "<<e.getMessage()<<'\n';return 2;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
