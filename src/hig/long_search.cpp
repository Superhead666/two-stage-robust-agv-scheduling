// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "long_search.hpp"
#include <fstream>
#include <set>

namespace longrun {
namespace {
double robust_energy(const Local& a, const std::vector<int>& block) {
  double sum = 0, left = a.gamma;
  std::vector<double> deviations;
  for (int j : block) { sum += a.nominal[j]; deviations.push_back(a.deviation[j]); }
  std::sort(deviations.rbegin(), deviations.rend());
  for (double d : deviations) { double f = std::clamp(left, 0.0, 1.0); sum += f*d; left -= f; }
  return sum;
}

double reference_horizon(const Instance& a, unsigned seed) {
  GRBEnv env(true);
  env.set(GRB_IntParam_OutputFlag, 0);
  env.set(GRB_IntParam_Threads, 1);
  env.set(GRB_IntParam_Seed, static_cast<int>(seed));
  env.start();
  GRBModel model(env);
  std::vector<std::vector<GRBVar>> x(a.n, std::vector<GRBVar>(a.m));
  std::vector<GRBVar> used(a.m), recharges(a.m);
  for (int k=0; k<a.m; ++k) {
    used[k]=model.addVar(0,1,0,GRB_BINARY);
    recharges[k]=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);
  }
  for (int j=0; j<a.n; ++j) for (int k=0; k<a.m; ++k)
    x[j][k]=model.addVar(0,1,0,GRB_BINARY);
  GRBVar horizon=model.addVar(0,GRB_INFINITY,1,GRB_CONTINUOUS);
  for (int j=0; j<a.n; ++j) {
    GRBLinExpr assigned=0;
    for (int k=0; k<a.m; ++k) {
      assigned+=x[j][k];
      model.addConstr(x[j][k]<=used[k]);
    }
    model.addConstr(assigned==1);
  }
  for (int k=0; k<a.m; ++k) {
    GRBLinExpr energy=0, workload=0;
    for (int j=0; j<a.n; ++j) {
      energy+=(a.nominal[j][k]/a.capacity[k])*x[j][k];
      workload+=a.p[j][k]*x[j][k];
    }
    model.addConstr(recharges[k]>=energy-used[k]);
    model.addConstr(horizon>=workload+a.alpha*recharges[k]);
  }
  model.set(GRB_IntAttr_ModelSense,GRB_MINIMIZE);
  model.optimize();
  if (model.get(GRB_IntAttr_Status)!=GRB_OPTIMAL)
    throw std::runtime_error("reference-horizon model was not solved to optimality");
  return model.get(GRB_DoubleAttr_ObjVal);
}
int canonical_value(const Local& a, std::vector<std::vector<int>> blocks) {
  struct B { int work, lateness; };
  std::vector<B> values;
  for (auto& block : blocks) {
    std::stable_sort(block.begin(), block.end(), [&](int i, int j) { return a.d[i] < a.d[j]; });
    int work = 0, lateness = -INF;
    for (int j : block) { work += a.p[j]; lateness = std::max(lateness, work-a.d[j]); }
    if (robust_energy(a, block) > a.capacity+ENERGY_EPS) throw std::runtime_error("invalid static block");
    values.push_back({work, lateness});
  }
  std::stable_sort(values.begin(), values.end(), [](B x, B y) {
    return x.work-x.lateness < y.work-y.lateness;
  });
  int time = 0, result = -INF;
  for (auto b : values) { result = std::max(result, time+b.lateness); time += b.work+a.alpha; }
  return result;
}
}

std::vector<int> Search::ranked(std::vector<int> v) const {
  if(opt.tardiness_score) for(int& value:v) value=std::max(0,value);
  std::sort(v.rbegin(),v.rend()); return v;
}

Options options_for(const std::string& mode) {
  Options o; o.mode = mode;
  if (mode == "manuscript") {
    o.manuscript=true;o.recombination=false;o.targeted=true;
    o.fast_screen=true;o.focus=true;o.block_removal=true;o.single_block=true;
    o.recourse_reference=true;o.ejection_chains=true;o.chain_after=true;
    o.plateau_acceptance=true;o.stop_zero=true;
    o.incremental_lower=true;o.memo_lower=true;o.decoder=2;
    o.local_limit=.1;
  }
  else if (mode == "static") { o.adaptive = false; o.recombination = false; }
  else if (mode == "descent") { o.recombination = false; o.iterated = false; }
  else if (mode == "adaptive") o.recombination = false;
  else if (mode == "packing") { o.recombination = false; o.packing = true; }
  else if (mode == "pool-packing") o.packing = true;
  else if (mode == "staged") { o.staged=true; o.targeted=true; o.inheritance=true; o.recombination=false; }
  else if (mode == "staged-pool") { o.staged=true; o.targeted=true; o.inheritance=true; }
  else if (mode == "staged-packing") { o.staged=true; o.targeted=true; o.inheritance=true; o.packing=true; o.recombination=false; }
  else if (mode == "alternating" || mode == "alternating-polished") {
    o.staged=true; o.targeted=true; o.inheritance=true; o.recombination=false;
    o.alternating=true; o.polishing=mode=="alternating-polished";
  }
  else if (mode == "bound-search") { o.staged=true; o.optimistic=true; o.inheritance=true; o.recombination=false; }
  else if (mode == "bound-packing") { o.staged=true; o.optimistic=true; o.inheritance=true; o.packing=true; o.recombination=false; }
  else if (mode == "bound-pool") { o.staged=true; o.optimistic=true; o.inheritance=true; }
  else if (mode == "balanced" || mode == "balanced-regret") {
    o.staged=true; o.optimistic=true; o.inheritance=true; o.recombination=false;
    o.lower_weight=3; o.regret_repair=mode=="balanced-regret";
  }
  else if (mode == "cyclic" || mode == "cyclic-regret" || mode == "random-regret" || mode == "polished") {
    o.staged=true; o.optimistic=true; o.inheritance=true; o.recombination=false;
    o.cyclic=mode!="random-regret"; o.regret_repair=mode!="cyclic";
    o.randomized_regret=o.regret_repair; o.lower_weight=mode=="random-regret"?3:4;
    o.polishing=mode=="polished";
  }
  else if (mode != "pool") throw std::runtime_error("unknown mode");
  return o;
}

Instance read_instance(const std::string& path) {
  std::ifstream f(path); Instance a; a.id = path;
  if (!(f >> a.n >> a.m >> a.alpha >> a.gamma) || a.n < 1 || a.n > 500 || a.m < 1 || a.m > 32)
    throw std::runtime_error("invalid instance header");
  a.capacity.resize(a.m); for (double& b : a.capacity) f >> b;
  a.d.resize(a.n); a.p.assign(a.n, std::vector<int>(a.m));
  a.nominal.assign(a.n, std::vector<double>(a.m)); a.deviation = a.nominal;
  for (int j=0; j<a.n; ++j) {
    f >> a.d[j];
    for (int k=0; k<a.m; ++k) f >> a.p[j][k] >> a.nominal[j][k] >> a.deviation[j][k];
  }
  if (!f || !std::isfinite(a.gamma) || a.gamma < 0 || a.gamma > a.n || a.alpha < 0 || a.alpha > 1000000)
    throw std::runtime_error("invalid instance data");
  for (double b : a.capacity) if (!std::isfinite(b) || b <= 0) throw std::runtime_error("invalid capacity");
  for (int j=0; j<a.n; ++j) for (int k=0; k<a.m; ++k)
    if (a.p[j][k]<=0 || a.p[j][k]>1000000 || std::abs(static_cast<long long>(a.d[j]))>1000000 ||
        !std::isfinite(a.nominal[j][k]) || !std::isfinite(a.deviation[j][k]) ||
        a.nominal[j][k]<0 || a.deviation[j][k]<0 ||
        a.nominal[j][k]+a.deviation[j][k]>a.capacity[k]+ENERGY_EPS)
      throw std::runtime_error("invalid job data or individual battery infeasibility");
  long long horizon=static_cast<long long>(a.n-1)*a.alpha,max_due=0;
  for(int j=0;j<a.n;++j) {
    horizon+=*std::max_element(a.p[j].begin(),a.p[j].end());
    max_due=std::max(max_due,std::abs(static_cast<long long>(a.d[j])));
  }
  if(horizon+max_due>INF/8) throw std::runtime_error("time magnitudes exceed safe integer range");
  return a;
}

Instance generate_paper(int n, int m, unsigned seed, double gamma, double variation,
                        double scale, double due_factor, int processing_max,
                        int nominal_energy_max, double recharge_factor,
                        int recharge_time, double tardiness_factor,
                        double due_date_range) {
  if (n<1 || n>500 || m<1 || m>32 || !std::isfinite(gamma) || !std::isfinite(variation) ||
      !std::isfinite(scale) || !std::isfinite(due_factor) || !std::isfinite(recharge_factor) ||
      !std::isfinite(tardiness_factor) || !std::isfinite(due_date_range) ||
      gamma<0 || gamma>n || variation<0 || scale<=0 || due_factor<0 ||
      processing_max<1 || nominal_energy_max<1 || recharge_factor<0 || recharge_time < -1 ||
      tardiness_factor < -1 || due_date_range < 0)
    throw std::runtime_error("invalid generator parameters");
  Instance a; a.n=n; a.m=m; a.gamma=gamma;
  a.id="paper_n"+std::to_string(n)+"m"+std::to_string(m)+"s"+std::to_string(seed);
  std::mt19937 r(seed); std::uniform_int_distribution<int> processing_uniform(1,processing_max);
  std::uniform_int_distribution<int> energy_uniform(1,nominal_energy_max);
  std::uniform_real_distribution<double> factor(.8,1.2), due(.25,1.1);
  a.p.assign(n,std::vector<int>(m)); a.nominal.assign(n,std::vector<double>(m)); a.deviation=a.nominal;
  a.capacity.assign(m,0); a.d.resize(n); double processing=0;
  for (int j=0;j<n;++j) for(int k=0;k<m;++k) {
    a.p[j][k]=processing_uniform(r); processing+=a.p[j][k]; a.nominal[j][k]=energy_uniform(r);
    a.deviation[j][k]=a.nominal[j][k]*variation*factor(r);
  }
  a.alpha=recharge_time>=0?recharge_time:
    static_cast<int>(std::round(recharge_factor*processing/(n*m)));
  for(int k=0;k<m;++k) {
    double sum=0,mx=0;
    for(int j=0;j<n;++j) { double e=a.nominal[j][k]+a.deviation[j][k]; sum+=e; mx=std::max(mx,e); }
    a.capacity[k]=std::max(mx,scale*sum/n);
    if(!std::isfinite(a.capacity[k])) throw std::runtime_error("generated capacity overflow");
  }
  // At a fixed seed, the base data and due-date quantiles remain common across
  // uncertainty-budget and due-date-factor levels.
  const bool use_tardiness_design=tardiness_factor>=0;
  const double horizon=use_tardiness_design?reference_horizon(a,seed):
    (processing/(n*m)+a.alpha/4.0)*n/m;
  std::uniform_real_distribution<double> unit(0,1);
  for(int j=0;j<n;++j) {
    double date=use_tardiness_design?
      std::round(horizon*(1-tardiness_factor-due_date_range/2+due_date_range*unit(r))):
      std::round(due_factor*horizon*due(r));
    date=std::max(0.0,date);
    if(!std::isfinite(date) || date>1000000) throw std::runtime_error("generated due date out of range");
    a.d[j]=static_cast<int>(date);
  }
  return a;
}

int static_upper(const Local& a, bool strong) {
  if (!a.n()) return -INF;
  std::vector<int> order(a.n()); std::iota(order.begin(),order.end(),0);
  int best=INF;
  for(int method=0;method<(strong?4:1);++method) {
    std::stable_sort(order.begin(),order.end(),[&](int i,int j) {
      double x=a.d[i],y=a.d[j];
      if(method==1) { x-=a.p[i]; y-=a.p[j]; }
      if(method==2) { x+=a.p[i]; y+=a.p[j]; }
      if(method==3) { x=a.nominal[i]+a.deviation[i]; y=a.nominal[j]+a.deviation[j]; return x>y; }
      return x<y;
    });
    std::vector<std::vector<int>> blocks;
    for(int j:order) {
      int selected=-1; double load=-1;
      int first=method==0?std::max(0,static_cast<int>(blocks.size())-1):0;
      for(int b=first;b<static_cast<int>(blocks.size());++b) {
        auto v=blocks[b]; v.push_back(j); double energy=robust_energy(a,v);
        if(energy<=a.capacity+ENERGY_EPS && energy>load) { load=energy; selected=b; }
      }
      if(selected<0) blocks.push_back({j}); else blocks[selected].push_back(j);
    }
    best=std::min(best,canonical_value(a,blocks));
  }
  return best;
}

int priority_upper(const Local& a,int rule) {
  if(!a.n())return -INF;
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){
    if(rule==1)return a.d[i]-a.p[i]<a.d[j]-a.p[j];
    if(rule==2)return a.d[i]+a.p[i]<a.d[j]+a.p[j];
    return a.nominal[i]+a.deviation[i]>a.nominal[j]+a.deviation[j];
  });
  std::vector<std::vector<int>> blocks;
  for(int j:order) {
    int selected=-1;double fill=-1;
    for(int b=0;b<int(blocks.size());++b) {
      auto block=blocks[b];block.push_back(j);double energy=robust_energy(a,block);
      if(energy<=a.capacity+ENERGY_EPS && energy>fill){fill=energy;selected=b;}
    }
    if(selected<0)blocks.push_back({j});else blocks[selected].push_back(j);
  }
  return canonical_value(a,blocks);
}

std::vector<std::vector<int>> static_reference(const Local& a, bool strong) {
  if(!a.n())return {};
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  int value=INF;std::vector<std::vector<int>> best;
  for(int method=0;method<(strong?4:1);++method) {
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){
      if(method==1)return a.d[i]-a.p[i]<a.d[j]-a.p[j];
      if(method==2)return a.d[i]+a.p[i]<a.d[j]+a.p[j];
      if(method==3)return a.nominal[i]+a.deviation[i]>a.nominal[j]+a.deviation[j];
      return a.d[i]<a.d[j];
    });
    std::vector<std::vector<int>> blocks;
    for(int j:order) {
      int selected=-1;double fill=-1;
      int first=method==0?std::max(0,int(blocks.size())-1):0;
      for(int b=first;b<int(blocks.size());++b) {
        auto candidate=blocks[b];candidate.push_back(j);double energy=robust_energy(a,candidate);
        if(energy<=a.capacity+ENERGY_EPS && energy>fill){fill=energy;selected=b;}
      }
      if(selected<0)blocks.push_back({j});else blocks[selected].push_back(j);
    }
    int candidate=canonical_value(a,blocks);
    if(candidate<value){value=candidate;best=blocks;}
  }
  for(auto& block:best)std::stable_sort(block.begin(),block.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  auto effective=[&](const std::vector<int>& block){
    int time=0,late=-INF;for(int j:block){time+=a.p[j];late=std::max(late,time-a.d[j]);}
    return time-late;
  };
  std::stable_sort(best.begin(),best.end(),[&](const auto& x,const auto& y){return effective(x)<effective(y);});
  return best;
}

int replacement_penalty(const Instance& a,int k,int removed,int inserted) {
  double cap=std::min(1.0,a.gamma);
  // Pointwise comparison of the two affine energy requirements over every
  // feasible value of the replaced job's uncertainty coordinate.
  if(a.nominal[inserted][k]>a.nominal[removed][k] ||
      a.nominal[inserted][k]+cap*a.deviation[inserted][k]>
      a.nominal[removed][k]+cap*a.deviation[removed][k])return INF;
  int delay=a.p[inserted][k]-a.p[removed][k];
  return std::max({0,delay,delay+a.d[removed]-a.d[inserted]});
}

namespace {
std::vector<std::vector<int>> insertion_partition(const Local& a) {
  if (!a.n()) return {};
  std::vector<int> order(a.n()); std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  std::vector<std::vector<int>> blocks;
  std::vector<std::pair<int,int>> block_values; // processing sum and within-block lateness
  auto partition_value=[&](std::vector<std::pair<int,int>> values) {
    std::stable_sort(values.begin(),values.end(),[](const auto& x,const auto& y){return x.first-x.second<y.first-y.second;});
    int time=0,result=-INF;
    for(auto b:values) {result=std::max(result,time+b.second);time+=b.first+a.alpha;}
    return result;
  };
  for (int j:order) {
    int chosen=-1, value=INF; double fill=-1;
    // Compare insertion into a feasible existing block with opening a block.
    // Every candidate uses the exact optimal order for that fixed partition.
    for (int b=0;b<=int(blocks.size());++b) {
      auto block=b==int(blocks.size())?std::vector<int>{}:blocks[b];block.push_back(j);
      double energy=robust_energy(a,block);
      if (energy>a.capacity+ENERGY_EPS) continue;
      auto values=block_values;
      if(b==int(blocks.size())) values.push_back({a.p[j],a.p[j]-a.d[j]});
      else {values[b].first+=a.p[j];values[b].second=std::max(values[b].second,values[b].first-a.d[j]);}
      int cost=partition_value(values);
      if (cost<value || (cost==value && energy>fill)) {
        value=cost; fill=energy; chosen=b;
      }
    }
    if (chosen<0) throw std::runtime_error("no feasible singleton in construction");
    if (chosen==int(blocks.size())) {blocks.push_back({j});block_values.push_back({a.p[j],a.p[j]-a.d[j]});}
    else {blocks[chosen].push_back(j);block_values[chosen].first+=a.p[j];block_values[chosen].second=std::max(block_values[chosen].second,block_values[chosen].first-a.d[j]);}
  }
  return blocks;
}
}

int insertion_upper(const Local& a) {
  if (!a.n()) return -INF;
  return canonical_value(a,insertion_partition(a));
}

std::vector<std::vector<int>> insertion_reference(const Local& a) {
  auto blocks=insertion_partition(a);
  for(auto& block:blocks)std::stable_sort(block.begin(),block.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  auto effective=[&](const std::vector<int>& block) {
    int time=0,late=-INF;for(int j:block){time+=a.p[j];late=std::max(late,time-a.d[j]);}
    return time-late;
  };
  std::stable_sort(blocks.begin(),blocks.end(),[&](const auto& x,const auto& y){return effective(x)<effective(y);});
  return blocks;
}

int slack_upper(const Local& a,long* improvements,long* gain) {
  int best=static_upper(a,false);
  const int initial=best;
  if(!a.n() || a.alpha==0) return best;
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  while(true) {
    const int target=best-1;
    std::vector<int> remaining=order;
    std::vector<std::vector<int>> blocks(1);
    int time=0;bool feasible=true;
    while(!remaining.empty()) {
      int j=remaining.front();auto trial=blocks.back();trial.push_back(j);
      if(robust_energy(a,trial)>a.capacity+ENERGY_EPS) {
        // Before recharging, fill otherwise unused battery capacity only when
        // this leaves enough slack for the earliest pending job after recharge.
        int filler=-1;double fill=-1;
        for(int h:remaining) {
          if(time+a.p[h]>a.d[h]+target || time+a.p[h]+a.alpha+a.p[j]>a.d[j]+target)continue;
          auto block=blocks.back();block.push_back(h);double energy=robust_energy(a,block);
          if(energy<=a.capacity+ENERGY_EPS && energy>fill) {fill=energy;filler=h;}
        }
        if(filler>=0) j=filler;
        else {time+=a.alpha;blocks.push_back({});}
      }
      time+=a.p[j];
      if(time>a.d[j]+target) {feasible=false;break;}
      blocks.back().push_back(j);remaining.erase(std::find(remaining.begin(),remaining.end(),j));
    }
    if(!feasible) break;
    int value=canonical_value(a,blocks);
    if(value>=best)throw std::runtime_error("slack construction did not improve target");
    best=value;
  }
  if(best<initial) {if(improvements)++*improvements;if(gain)*gain+=initial-best;}
  return best;
}

int workload_lower(const Local& a, bool packing) {
  int answer=-INF,work=0;
  std::vector<int> order(a.n()),prefix; std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  for(int j:order) {
    prefix.push_back(j); work+=a.p[j];
    int bins=std::max(1,static_cast<int>(std::ceil((robust_energy(a,prefix)-ENERGY_EPS)/a.capacity)));
    if(packing) for(int q=1;q<=3;++q) {
      double threshold=a.capacity/(q+1),sum=0;
      std::vector<double> costs;
      for(int h:prefix) {
        if(a.nominal[h]>threshold+ENERGY_EPS) costs.push_back(0);
        else if(a.deviation[h]>0 && a.nominal[h]+a.deviation[h]>threshold+1e-7) {
          double need=(threshold+1e-7-a.nominal[h])/a.deviation[h];
          if(need<=1) costs.push_back(std::max(0.0,need));
        }
      }
      std::sort(costs.begin(),costs.end()); int count=0;
      for(double c:costs) { if(sum+c>a.gamma) break; sum+=c; ++count; }
      bins=std::max(bins,(count+q-1)/q);
    }
    answer=std::max(answer,work+(bins-1)*a.alpha-a.d[j]);
  }
  return answer;
}

Search::Search(const Instance& data, GRBEnv& environment, Options options)
    :a(data),env(environment),opt(options),budget(options.limit),rng(options.seed),cache(a.m),lower_cache(a.m) {
  for(int j=0;j<a.n;++j)
    global_lower=std::max(global_lower,*std::min_element(a.p[j].begin(),a.p[j].end())-a.d[j]);
}
std::vector<std::vector<int>> Search::groups(const std::vector<int>& x) const {
  std::vector<std::vector<int>> result(a.m);
  for(int j=0;j<a.n;++j) if(x[j]>=0) result[x[j]].push_back(j);
  return result;
}
Entry& Search::entry(int k,const std::vector<int>& jobs) {
  auto found=cache[k].find(jobs);
  if(found!=cache[k].end()) { ++stats.cache_hits; return *found->second; }
  Clock timing;
  auto item=std::make_unique<Entry>(); auto& e=*item;
  e.data.capacity=a.capacity[k]; e.data.gamma=a.gamma; e.data.alpha=a.alpha;
  for(int j:jobs) {
    e.data.ids.push_back(j); e.data.p.push_back(a.p[j][k]); e.data.d.push_back(a.d[j]);
    e.data.e.push_back(a.nominal[j][k]); e.data.nominal.push_back(a.nominal[j][k]);
    e.data.deviation.push_back(a.deviation[j][k]);
  }
  Clock lower_clock;
  auto cached_lower=lower_cache[k].find(jobs);
  if(opt.memo_lower && cached_lower!=lower_cache[k].end()) {
    e.lower=cached_lower->second;lower_cache[k].erase(cached_lower);++stats.lower_hits;
  } else e.lower=opt.incremental_lower?incremental_lower(e.data):workload_lower(e.data,opt.packing);
  stats.lower_seconds+=lower_clock.elapsed();
  if(opt.decoder>=4) e.static_upper=priority_upper(e.data,opt.decoder==4?1:opt.decoder==5?3:2);
  else if(opt.decoder==3)e.static_upper=slack_upper(e.data,&stats.slack_improvements,&stats.slack_gain);
  else if(opt.decoder==2)e.static_upper=insertion_upper(e.data);
  else e.static_upper=static_upper(e.data,opt.decoder==1?false:opt.strong_static);
  e.upper=e.static_upper;
  if(opt.inheritance) for(std::size_t h=0;h<jobs.size();++h) {
    auto subset=jobs; subset.erase(subset.begin()+h); auto it=cache[k].find(subset);
    if(it!=cache[k].end()) e.lower=std::max(e.lower,it->second->lower);
  }
  if(e.lower>e.upper) throw std::runtime_error("local bounds inconsistent");
  auto* pointer=item.get(); cache[k].emplace(jobs,std::move(item)); ++stats.entries;
  stats.entry_seconds+=timing.elapsed();
  return *pointer;
}
void Search::refine(Entry& e,double seconds,int target) {
  if(opt.interval_targets)target=INF;
  if(e.lower>=e.upper || e.data.n()>62 || budget.expired()) return;
  // A nonpositive upper bound certifies zero tardiness, not exact signed lateness.
  // Leave the signed interval intact for bound transfer and inexpensive guidance.
  if(opt.zero_aware_refinement && e.upper<=0) {++stats.zero_skipped;return;}
  if(opt.zero_aware_refinement && target<0) {target=0;++stats.targets_capped;}
  seconds=std::min(seconds,budget.left()); if(seconds<.001) return;
  if(target==INF && e.largest_failed_budget>=seconds*.99) return;
  if(target!=INF && e.failed_target==target && e.failed_target_budget>=seconds*.99) return;
  Clock local(seconds);
  if(opt.budget_portfolio && !e.budget_seeded) {
    Clock pt;auto seeds=budget_repair_portfolio(e.data,opt.budget_best);int certified=-INF;
    for(auto& s:seeds){certified=std::max(certified,s.lateness);e.schedules.push_back(std::move(s));}
    if(certified<e.upper){stats.budget_gain+=e.upper-certified;e.upper=certified;++stats.budget_updates;}
    if(e.lower>e.upper)throw std::runtime_error("budget portfolio contradicts lower bound");
    e.budget_seeded=true;++stats.portfolio_entries;stats.portfolio_schedules+=seeds.size();stats.portfolio_seconds+=pt.elapsed();
    if(e.lower==e.upper || (target!=INF && target>=e.upper))return;
  }
  if(opt.seed_portfolio && !e.portfolio_seeded) {
    Clock pt;
    auto seeds=nominal_portfolio(e.data);
    for(auto& s:seeds) e.schedules.push_back(std::move(s));
    e.portfolio_seeded=true;++stats.portfolio_entries;stats.portfolio_schedules+=seeds.size();
    stats.portfolio_seconds+=pt.elapsed();
  }
  RobustOptions ro; ro.workload_screen=true; ro.refine=true; ro.reuse_schedules=true;
  // The componentwise box is exact only when the budget can activate every
  // assigned deviation. Avoid a size threshold or a conservative box solve
  // in the genuinely budgeted regime.
  ro.screening=e.data.gamma>=e.data.n();
  while(e.lower<e.upper && !local.expired() && !budget.expired()) {
    if(opt.zero_aware_refinement && e.upper<=0) break;
    int bound=target==INF?e.lower+(e.upper-e.lower)/2:std::clamp(target,e.lower,e.upper-1);
    if(opt.zero_aware_refinement && bound<0) {bound=0;++stats.targets_capped;}
    if(opt.witness_memory && !witnesses.empty()) {
      Clock wt; std::vector<std::pair<double,int>> rank;
      for(int s=0;s<int(witnesses.size());++s) { double excess=0;
        for(int j=0;j<e.data.n();++j) excess+=e.data.deviation[j]*witnesses[s][e.data.ids[j]];
        rank.push_back({excess,s}); }
      std::sort(rank.rbegin(),rank.rend()); bool adverse=false;
      for(int h=0;h<std::min(3,int(rank.size())) && local.left()>.001;++h) {
        std::vector<double> xi;for(int id:e.data.ids)xi.push_back(witnesses[rank[h].second][id]);
        Clock det(std::min(.001,local.left())); auto r=deterministic(realize(e.data,xi),bound,DetOptions{},det);++stats.witness_checks;
        if(r.status==Status::Infeasible) {e.lower=bound+1;++stats.witness_prunes;adverse=true;break;}
        if(r.status==Status::Feasible && e.schedules.size()<256) e.schedules.push_back(r.schedule);
      }
      stats.witness_seconds+=wt.elapsed();
      if(adverse) {if(target!=INF)break;else continue;}
    }
    if(local.expired()) break;
    auto result=opt.tight_coverage?tight_verify(e.data,bound,env,ro,local,&e.schedules):verify(e.data,bound,env,ro,local,&e.schedules);
    ++stats.verifications; ++e.attempts; stats.verification_seconds+=result.seconds;
    if(bound<0) {++stats.negative_targets;stats.negative_target_seconds+=result.seconds;}
    if(bound==0) {++stats.zero_targets;stats.zero_target_seconds+=result.seconds;}
    stats.adversary_seconds+=result.mip_seconds;
    if(result.status==Status::Feasible) {
      e.upper=bound; ++stats.adaptive_improvements;
    } else if(result.status==Status::Infeasible) {
      e.lower=bound+1;
      if(opt.witness_memory && result.witness.size()==std::size_t(e.data.n())) {
        std::vector<double> xi(a.n,0); for(int j=0;j<e.data.n();++j)xi[e.data.ids[j]]=std::clamp(result.witness[j],0.0,1.0);
        double sum=std::accumulate(xi.begin(),xi.end(),0.0);if(sum>a.gamma)for(double& z:xi)z*=a.gamma/sum;
        if(std::find(witnesses.begin(),witnesses.end(),xi)==witnesses.end()) {witnesses.push_back(xi);if(witnesses.size()>32)witnesses.erase(witnesses.begin());}
      }
    }
    else {
      ++stats.unknown; e.largest_failed_budget=std::max(e.largest_failed_budget,seconds);
      if(result.reason=="scheduling_limit") ++stats.det_unknown; else ++stats.coverage_unknown;
      if(result.reason=="numerical_boundary") ++stats.numerical_unknown;
      e.failed_target=target; e.failed_target_budget=seconds; break;
    }
    if(e.schedules.size()>256) e.schedules.erase(e.schedules.begin(),e.schedules.end()-256);
    if(target!=INF) break;
  }
  if(e.lower==e.upper) ++stats.exact_sets;
  if(e.lower>e.upper) throw std::runtime_error("adaptive bounds inconsistent");
}
std::vector<int> Search::score(const std::vector<int>& x,bool improve,const std::vector<int>* incumbent) {
  ++stats.evaluations; auto sets=groups(x);
  std::vector<Entry*> entries; std::vector<int> lower,upper;
  for(int k=0;k<a.m;++k) {
    auto& e=entry(k,sets[k]); entries.push_back(&e); lower.push_back(e.lower); upper.push_back(e.upper);
  }
  if(incumbent && ranked(lower)>=*incumbent) { ++stats.screened; return ranked(upper); }
  if(improve && opt.adaptive) {
    int comparison_target=incumbent?incumbent->front():INF;
    if(opt.manuscript && incumbent) {
      // An expensive check must be able to improve the best HIG objective,
      // not merely the evaluation of an intermediate local-search assignment.
      if(!best_score.empty())comparison_target=std::min(comparison_target,best_score.front());
      --comparison_target;
      if(*std::max_element(lower.begin(),lower.end())>comparison_target)return ranked(upper);
    }
    std::sort(entries.begin(),entries.end(),[](Entry* x,Entry* y){return x->upper>y->upper;});
    for(auto* e:entries) {
      if(budget.expired()) break;
      if(!incumbent || e->lower<=comparison_target) {
        if(opt.targeted && incumbent) {
          if(e->upper>comparison_target) refine(*e,opt.local_limit,comparison_target);
        } else refine(*e,opt.local_limit);
        // This excludes a strict incumbent improvement. Earlier certified
        // updates can still make the candidate preferable to the current x.
        if(opt.manuscript && incumbent && e->lower>comparison_target) break;
      }
    }
    upper.clear(); for(int k=0;k<a.m;++k) upper.push_back(entry(k,sets[k]).upper);
  }
  return ranked(upper);
}
std::vector<int> Search::guide(const std::vector<int>& x) {
  if(!opt.optimistic || !opt.adaptive) return score(x);
  auto sets=groups(x); std::vector<int> lower,upper;
  for(int k=0;k<a.m;++k) {
    auto& e=entry(k,sets[k]);
    lower.push_back(e.data.n()?opt.lower_weight*e.lower+(4-opt.lower_weight)*e.upper:-INF);
    upper.push_back(e.upper);
  }
  lower=ranked(lower); upper=ranked(upper);
  lower.insert(lower.end(),upper.begin(),upper.end()); return lower;
}
void Search::remember(const std::vector<int>& x) {
  if(std::find(x.begin(),x.end(),-1)!=x.end()) return;
  auto s=score(x);
  if(!best.empty()) {
    auto refreshed=score(best);
    if(refreshed<best_score) {
      if(refreshed.front()<best_score.front()) time_best=budget.elapsed();
      best_score=refreshed;
    }
  }
  if(best.empty() || s<best_score) {
    if(best.empty() || s.front()<best_score.front()) time_best=budget.elapsed();
    best=x; best_score=s;
  }
  if(opt.progress) opt.progress(budget.elapsed(),best_score.front(),global_lower,stats,best);
}
std::vector<int> Search::construct(int perturbation) {
  std::vector<int> order(a.n),x(a.n,-1); std::iota(order.begin(),order.end(),0);
  if(!opt.manuscript)std::shuffle(order.begin(),order.end(),rng);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){
    if(perturbation==1) return *std::min_element(a.p[i].begin(),a.p[i].end())>*std::min_element(a.p[j].begin(),a.p[j].end());
    return a.d[i]<a.d[j];
  });
  for(int j:order) {
    std::vector<int> best_value; int chosen=0;
    for(int k=0;k<a.m;++k) { x[j]=k; auto v=guide(x); if(best_value.empty()||v<best_value){best_value=v;chosen=k;} }
    x[j]=chosen;
  }
  return x;
}
bool Search::descend(std::vector<int>& x) {
  bool changed=false; int accepted=0;
  while(!search_expired() && budget.elapsed()<opt.phase_end && (opt.manuscript || accepted<200 || (opt.event_warm && !opt.adaptive))) {
    if(opt.ejection_chains && !opt.chain_after) {
      if(chain_descent(x)) {changed=true;++accepted;continue;}
      if(opt.chain_only) break;
    }
    auto incumbent=opt.optimistic?guide(x):score(x,!opt.manuscript); remember(x);
    std::vector<std::pair<int,int>> moves;
    if(opt.relocation_moves)
      for(int j=0;j<a.n;++j) for(int k=0;k<a.m;++k) if(x[j]!=k) moves.push_back({j,-k-1});
    if(opt.exchange_moves)
      for(int i=0;i<a.n;++i) for(int j=i+1;j<a.n;++j) if(x[i]!=x[j]) moves.push_back({i,j});
    std::shuffle(moves.begin(),moves.end(),rng);
    std::vector<bool> critical(a.m,false);
    if(opt.focus) {
      auto sets=groups(x);for(int k=0;k<a.m;++k)critical[k]=entry(k,sets[k]).upper==incumbent.front();
      auto touches=[&](const auto& move){return critical[x[move.first]] || critical[move.second<0?-move.second-1:x[move.second]];};
      std::stable_sort(moves.begin(),moves.end(),[&](const auto& u,const auto& v){return touches(u)>touches(v);});
      if(opt.critical_prefix) {
        std::vector<bool> preferred(a.n,false);
        for(int k=0;k<a.m;++k)if(critical[k] && !sets[k].empty()) {
          auto& e=entry(k,sets[k]);
          std::vector<int> order;std::vector<bool> recharge;
          const Schedule* chosen=nullptr;
          for(const auto& s:e.schedules)if(s.lateness<=e.upper && (!chosen || s.lateness>chosen->lateness))chosen=&s;
          if(chosen){order=chosen->order;recharge=chosen->recharge;}
          else for(const auto& b:static_reference(e.data))for(int j:b){recharge.push_back(order.empty()?false:j==b.front());order.push_back(j);}
          int t=0,late=-INF,last=-1;
          for(int h=0;h<int(order.size());++h){if(recharge[h])t+=a.alpha;t+=e.data.p[order[h]];if(t-e.data.d[order[h]]>late){late=t-e.data.d[order[h]];last=h;}}
          for(int h=0;h<=last;++h)preferred[e.data.ids[order[h]]]=true;
          stats.prioritized_jobs+=last+1;stats.critical_jobs+=order.size();
        }
        ++stats.priority_scans;
        auto priority=[&](const auto& move){return preferred[move.first] || (move.second>=0 && preferred[move.second]);};
        // Priority only. An adaptive optimal schedule need not have this prefix.
        std::stable_sort(moves.begin(),moves.end(),[&](const auto& u,const auto& v){return priority(u)>priority(v);});
      }
    }
    bool found=false;
    for(auto move:moves) {
      if(search_expired() || budget.elapsed()>=opt.phase_end) break;
      auto candidate=x;
      if(move.second<0) candidate[move.first]=-move.second-1;
      else std::swap(candidate[move.first],candidate[move.second]);
      if(opt.fast_screen && cheap_screen(candidate,incumbent)) continue;
      if(opt.replacement_bounds && move.second>=0) {
        auto before=groups(x),after=groups(candidate);
        for(int j:{move.first,move.second}) {
          int replacement=j==move.first?move.second:move.first,k=x[j];
          int upper_shift=replacement_penalty(a,k,j,replacement);
          int lower_shift=replacement_penalty(a,k,replacement,j);
          if(upper_shift==INF && lower_shift==INF)continue;
          auto& parent=entry(k,before[k]);auto& child=entry(k,after[k]);
          if(upper_shift<INF && child.upper>parent.upper+upper_shift) {
            stats.replacement_gain+=child.upper-parent.upper-upper_shift;
            child.upper=parent.upper+upper_shift;++stats.replacement_upper;
          }
          if(lower_shift<INF && child.lower<parent.lower-lower_shift) {
            child.lower=parent.lower-lower_shift;++stats.replacement_lower;
          }
          if(child.lower>child.upper)throw std::runtime_error("replacement bounds inconsistent");
        }
        if(opt.fast_screen && cheap_screen(candidate,incumbent))continue;
      }
      if(opt.inheritance && move.second<0) {
        auto previous=groups(x),next=groups(candidate); int k=x[move.first];
        auto& removed=entry(k,next[k]); removed.upper=std::min(removed.upper,entry(k,previous[k]).upper);
      }
      auto value=opt.optimistic?guide(candidate):score(candidate,false,&incumbent);
      bool touches=!opt.focus || critical[x[move.first]] || critical[move.second<0?-move.second-1:x[move.second]];
      if(!opt.optimistic && opt.adaptive && value>=incumbent && (opt.manuscript || touches)) value=score(candidate,true,&incumbent);
      if(opt.managed && opt.adaptive && !best.empty() && value.front()<=best_score.front()+a.alpha && elite.size()<24) {
        auto sets=groups(candidate); int lo=-INF; bool unresolved=false;
        for(int k=0;k<a.m;++k) { auto& e=entry(k,sets[k]); lo=std::max(lo,e.lower); unresolved|=e.lower<e.upper; }
        if(unresolved && lo<best_score.front() && std::find(elite.begin(),elite.end(),candidate)==elite.end()) elite.push_back(candidate);
      }
      // Cached bounds may have improved while examining another assignment.
      incumbent=guide(x);
      if(value<incumbent) {
        x=candidate; ++stats.moves; ++accepted;
        if(move.second<0)++stats.relocation_moves;else ++stats.exchange_moves;
        changed=found=true; remember(x); break;
      }
    }
    if(!found) {
      if(opt.chain_after && !search_expired() && budget.elapsed()<opt.phase_end) {
        ++stats.chain_stalls;
        if(chain_descent(x)){changed=true;++accepted;continue;}
      }
      break;
    }
  }
  if(opt.optimistic && opt.adaptive) { if(!search_expired()) score(x,true); remember(x); }
  return changed;
}
std::vector<int> Search::shake(const std::vector<int>& x,int strength) {
  if(opt.block_removal)return block_shake(x);
  if(opt.pair_repair) return pair_shake(x);
  if(opt.marginal_destruction || opt.regret_insertion || opt.energy_destruction) return structured_shake(x,strength);
  auto next=x; std::vector<int> jobs(a.n); std::iota(jobs.begin(),jobs.end(),0);
  std::shuffle(jobs.begin(),jobs.end(),rng);
  if(opt.regret_repair && stats.shakes%3) {
    auto sets=groups(x); int critical=0;
    for(int k=1;k<a.m;++k) if(entry(k,sets[k]).upper>entry(critical,sets[critical]).upper) critical=k;
    if(stats.shakes%3==1) std::stable_sort(jobs.begin(),jobs.end(),[&](int i,int j){return (x[i]==critical)>(x[j]==critical);});
    else std::stable_sort(jobs.begin(),jobs.end(),[&](int i,int j) {
      int ci=a.p[i][x[i]]-*std::min_element(a.p[i].begin(),a.p[i].end());
      int cj=a.p[j][x[j]]-*std::min_element(a.p[j].begin(),a.p[j].end()); return ci>cj;
    });
  }
  jobs.resize(std::min(a.n,strength));
  for(int j:jobs) next[j]=-1;
  if(opt.monotone_deletion) inherit_removals(x,next);
  if(opt.deadline_repair)std::stable_sort(jobs.begin(),jobs.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  if(opt.regret_repair) {
    while(!jobs.empty()) {
      int selected=0,vehicle=0; long long regret=-1;
      for(int h=0;h<static_cast<int>(jobs.size());++h) {
        int j=jobs[h]; std::vector<std::pair<std::vector<int>,int>> choices;
        for(int k=0;k<a.m;++k) { next[j]=k; choices.push_back({guide(next),k}); }
        next[j]=-1; std::sort(choices.begin(),choices.end());
        long long value=a.m<2?0:static_cast<long long>(choices[1].first.front())-choices[0].first.front();
        if(value>regret) {
          regret=value; selected=h;
          int rank=opt.randomized_regret && a.m>1 && rng()%8==0?1:0;
          vehicle=choices[rank].second;
        }
      }
      next[jobs[selected]]=vehicle; jobs.erase(jobs.begin()+selected);
    }
    ++stats.shakes; return next;
  }
  for(int j:jobs) {
    std::vector<std::pair<std::vector<int>,int>> choices;
    for(int k=0;k<a.m;++k) { next[j]=k; choices.push_back({guide(next),k}); }
    std::sort(choices.begin(),choices.end());
    int rank=(a.m>1 && rng()%4==0)?1:0; next[j]=choices[rank].second;
  }
  ++stats.shakes; return next;
}

std::vector<int> Search::block_shake(const std::vector<int>& x) {
  Clock timing;++stats.block_calls;
  auto finish=[&](const std::vector<int>& y){
    stats.block_seconds+=timing.elapsed();++stats.shakes;
    for(int j=0;j<a.n;++j)stats.block_changed+=x[j]!=y[j];
    return y;
  };
  if(a.m<2 || search_expired())return finish(x);
  auto sets=groups(x);auto current=score(x);std::vector<int> critical;
  for(int k=0;k<a.m;++k)if(!sets[k].empty() && entry(k,sets[k]).upper==current.front())critical.push_back(k);
  if(critical.empty())return finish(x);
  int seed_job=opt.related_removal?int(rng()%a.n):-1;
  int source=opt.related_removal?x[seed_job]:critical[rng()%critical.size()];
  int partner=-1;
  if(!opt.related_removal && !opt.single_block){partner=int(rng()%(a.m-1));if(partner>=source)++partner;}
  auto reference=[&](Entry& e){
    const Schedule* selected=nullptr;int value=-INF,ties=0;
    if(opt.recourse_reference)for(const auto& s:e.schedules) {
      if(s.order.size()!=std::size_t(e.data.n()) || s.recharge.size()!=s.order.size() || s.lateness>e.upper)continue;
      if(s.lateness>value){selected=&s;value=s.lateness;ties=1;}
      else if(s.lateness==value && rng()%++ties==0)selected=&s;
    }
    if(selected) {
      std::vector<std::vector<int>> blocks;
      for(int h=0;h<int(selected->order.size());++h) {
        if(h==0 || selected->recharge[h])blocks.push_back({});
        blocks.back().push_back(selected->order[h]);
      }
      ++stats.reference_recourse;return blocks;
    }
    ++stats.reference_static;
    if(opt.decoder==2)return opt.manuscript?insertion_reference(e.data):insertion_partition(e.data);
    return static_reference(e.data,opt.decoder!=1);
  };
  auto& origin=entry(source,sets[source]);auto blocks=reference(origin);
  int time=0,late=-INF,chosen=0,pivot=-1;
  for(int b=0;b<int(blocks.size());++b) {
    if(b)time+=a.alpha;
    for(int j:blocks[b]) {
      time+=origin.data.p[j];int candidate=time-origin.data.d[j];
      if(candidate>late){late=candidate;chosen=b;pivot=origin.data.ids[j];}
    }
  }
  if(opt.related_removal) {
    pivot=seed_job;
    for(int b=0;b<int(blocks.size());++b)for(int j:blocks[b])if(origin.data.ids[j]==pivot)chosen=b;
  }
  if(pivot<0)return finish(x);
  auto effective=[&](const Local& local,const std::vector<int>& block){
    int processing=0,deviation=-INF;
    for(int j:block){processing+=local.p[j];deviation=std::max(deviation,processing-local.d[j]);}
    return processing-deviation;
  };
  int date=effective(origin.data,blocks[chosen]);std::vector<int> removed;
  for(int j:blocks[chosen])removed.push_back(origin.data.ids[j]);
  int distance=INF,ties=0;std::vector<int> matched;
  if(!opt.single_block) {
    for(int k=0;k<a.m;++k)if(k!=source && !sets[k].empty() && (opt.related_removal || k==partner)) {
      auto& receiving=entry(k,sets[k]);auto alternatives=reference(receiving);
      for(const auto& block:alternatives) {
        int delta=std::abs(effective(receiving.data,block)-date);bool take=false;
        if(delta<distance){distance=delta;ties=1;take=true;}
        else if(delta==distance && rng()%++ties==0)take=true;
        if(take){matched.clear();for(int j:block)matched.push_back(receiving.data.ids[j]);}
      }
    }
  }
  removed.insert(removed.end(),matched.begin(),matched.end());
  stats.block_removed+=removed.size();auto next=x;
  for(int j:removed)next[j]=-1;
  if(!opt.manuscript)std::shuffle(removed.begin(),removed.end(),rng);
  std::stable_sort(removed.begin(),removed.end(),[&](int i,int j){
    return a.d[i]!=a.d[j]?a.d[i]<a.d[j]:(opt.manuscript && i<j);
  });
  auto choices_for=[&](int j){
    std::vector<std::pair<std::vector<int>,int>> choices;
    for(int k=0;k<a.m;++k)if(j!=pivot || k!=source) {
      next[j]=k;choices.push_back({score(next),k});++stats.repair_candidates;
    }
    next[j]=-1;std::sort(choices.begin(),choices.end());return choices;
  };
  while(!removed.empty()) {
    if(search_expired())return finish(x);
    int selected=removed.front(),vehicle=-1,scarcity=INF;
    std::vector<std::pair<std::vector<int>,int>> selected_choices;
    if(opt.scarcity_repair) {
      auto partial=groups(next);
      for(int j:removed) {
        if(search_expired())return finish(x);
        auto choices=choices_for(j);int available=0;
        for(const auto& choice:choices) {
          int k=choice.second;auto jobs=partial[k];jobs.insert(std::lower_bound(jobs.begin(),jobs.end(),j),j);
          available+=entry(k,jobs).lower<=current.front();++stats.scarcity_queries;
        }
        if(available<scarcity){scarcity=available;selected=j;selected_choices=std::move(choices);}
      }
    } else selected_choices=choices_for(selected);
    if(selected_choices.empty())return finish(x);
    // Randomness only resolves equal comparison vectors, not a fixed chance
    // of deliberately taking the second-best insertion.
    int tied=1;while(tied<int(selected_choices.size()) && selected_choices[tied].first==selected_choices[0].first)++tied;
    vehicle=selected_choices[rng()%tied].second;next[selected]=vehicle;
    removed.erase(std::find(removed.begin(),removed.end(),selected));
  }
  if(next[pivot]==source)throw std::runtime_error("block pivot did not change vehicle");
  return finish(next);
}

void Search::inherit_removals(const std::vector<int>& before,const std::vector<int>& after) {
  auto parent=groups(before),child=groups(after);
  for(int k=0;k<a.m;++k) if(child[k].size()<parent[k].size() &&
      std::includes(parent[k].begin(),parent[k].end(),child[k].begin(),child[k].end())) {
    int bound=entry(k,parent[k]).upper;auto& e=entry(k,child[k]);
    if(e.upper>bound) {++stats.deletion_updates;stats.deletion_gain+=e.upper-bound;e.upper=bound;}
    if(e.lower>e.upper)throw std::runtime_error("deletion inheritance contradiction");
  }
}

std::vector<int> Search::structured_shake(const std::vector<int>& x,int strength) {
  Clock timing;
  auto next=x; std::vector<int> removed;
  int pivot=-1, origin=-1;
  if(opt.marginal_destruction || opt.energy_destruction) {
    auto sets=groups(x); auto current=score(x);
    std::vector<int> critical;
    for(int k=0;k<a.m;++k) if(entry(k,sets[k]).upper==current.front() && !sets[k].empty()) critical.push_back(k);
    if(critical.empty()) return x;
    origin=critical[rng()%critical.size()];
    const int target=current.front()-1;
    auto remaining=sets[origin];
    const Local local=entry(origin,sets[origin]).data;
    std::vector<int> energy_prefix;
    double energy_target=-1;
    if(opt.energy_destruction && a.alpha>0) {
      std::vector<int> order(local.n()),prefix,binding;std::iota(order.begin(),order.end(),0);
      std::stable_sort(order.begin(),order.end(),[&](int i,int j){return local.d[i]<local.d[j];});
      int processing=0,value=-INF;
      for(int j:order) {
        prefix.push_back(j);processing+=local.p[j];
        int bins=std::max(1,int(std::ceil((robust_energy(local,prefix)-ENERGY_EPS)/local.capacity)));
        int bound=processing+a.alpha*(bins-1)-local.d[j];
        if(bound>value) {value=bound;binding=prefix;}
      }
      int bins=std::max(1,int(std::ceil((robust_energy(local,binding)-ENERGY_EPS)/local.capacity)));
      if(bins>1) {energy_prefix=binding;energy_target=(bins-1)*local.capacity;}
    }
    while(!remaining.empty()) {
      if(budget.expired()) return x;
      std::shuffle(remaining.begin(),remaining.end(),rng);
      std::vector<int> best_value; int selected=-1;double best_energy=1e100;int best_penalty=INF;
      for(int j:remaining) {
        if(budget.expired()) return x;
        if(energy_target>=0) {
          int index=int(std::lower_bound(local.ids.begin(),local.ids.end(),j)-local.ids.begin());
          auto it=std::find(energy_prefix.begin(),energy_prefix.end(),index);if(it==energy_prefix.end())continue;
          auto subset=energy_prefix;subset.erase(subset.begin()+(it-energy_prefix.begin()));
          double energy=robust_energy(local,subset);int penalty=INF;
          for(int k=0;k<a.m;++k)if(k!=origin)penalty=std::min(penalty,a.p[j][k]-a.p[j][origin]);
          if(energy<best_energy-ENERGY_EPS || (std::abs(energy-best_energy)<=ENERGY_EPS && penalty<best_penalty)) {
            selected=j;best_energy=energy;best_penalty=penalty;
          }
        } else {
          auto parent=next;next[j]=-1;
          if(opt.monotone_deletion)inherit_removals(parent,next);
          auto value=score(next); next[j]=origin;
          if(selected<0 || value<best_value) {selected=j;best_value=value;}
        }
      }
      if(selected<0) break;
      if(pivot<0) pivot=selected;
      next[selected]=-1; removed.push_back(selected);
      remaining.erase(std::find(remaining.begin(),remaining.end(),selected));
      std::sort(remaining.begin(),remaining.end());
      if(energy_target>=0) {
        int index=int(std::lower_bound(local.ids.begin(),local.ids.end(),selected)-local.ids.begin());
        energy_prefix.erase(std::find(energy_prefix.begin(),energy_prefix.end(),index));
        if(robust_energy(local,energy_prefix)<=energy_target+ENERGY_EPS)break;
      } else if(entry(origin,remaining).upper<=target) break;
    }
  } else {
    removed.resize(a.n);std::iota(removed.begin(),removed.end(),0);
    std::shuffle(removed.begin(),removed.end(),rng);removed.resize(std::min(a.n,strength));
    for(int j:removed) next[j]=-1;
  }
  stats.removed_jobs+=removed.size();
  if(opt.monotone_deletion) inherit_removals(x,next);
  // For targeted removal the pivot must leave its original machine. Other
  // removed jobs may return. This produces a true perturbation without a
  // tuned probability of choosing the second-best insertion.
  auto choices_for=[&](int j) {
    std::vector<std::pair<std::vector<int>,int>> choices;
    for(int k=0;k<a.m;++k) {
      if(a.m>1 && (j==pivot || opt.energy_destruction) && k==origin) continue;
      next[j]=k;choices.push_back({score(next),k});++stats.repair_candidates;
    }
    next[j]=-1;std::sort(choices.begin(),choices.end());return choices;
  };
  if(opt.regret_insertion) {
    while(!removed.empty()) {
      if(budget.expired()) return x;
      int selected=-1,vehicle=-1;std::vector<long long> largest;
      for(int j:removed) {
        if(budget.expired()) return x;
        auto choices=choices_for(j);std::vector<long long> regret(a.m,0);
        if(choices.size()==1) regret[0]=static_cast<long long>(INF)*2;
        else for(int k=0;k<a.m;++k) regret[k]=static_cast<long long>(choices[1].first[k])-choices[0].first[k];
        if(selected<0 || regret>largest || (regret==largest && a.d[j]<a.d[selected])) {
          selected=j;vehicle=choices.front().second;largest=regret;
        }
      }
      next[selected]=vehicle;removed.erase(std::find(removed.begin(),removed.end(),selected));
    }
  } else {
    std::shuffle(removed.begin(),removed.end(),rng);
    for(int j:removed) {
      if(budget.expired()) return x;
      auto choices=choices_for(j);
      int rank=choices.size()>1 && rng()%4==0?1:0;next[j]=choices[rank].second;
    }
  }
  ++stats.shakes;stats.repair_seconds+=timing.elapsed();return next;
}

std::vector<int> Search::pair_shake(const std::vector<int>& x) {
  if(a.m<2 || budget.left()<.001) return x;
  Clock local(std::min(opt.pulse_limit,budget.left()));++stats.pair_calls;
  auto finish=[&](const std::vector<int>& y){stats.pair_seconds+=local.elapsed();++stats.shakes;return y;};
  auto sets=groups(x);auto incumbent=score(x);std::vector<int> critical;
  for(int k=0;k<a.m;++k) if(entry(k,sets[k]).upper==incumbent.front()) critical.push_back(k);
  int first=critical[rng()%critical.size()],second=rng()%(a.m-1);
  if(second>=first) ++second;
  std::vector<int> jobs=sets[first];jobs.insert(jobs.end(),sets[second].begin(),sets[second].end());
  if(jobs.empty()) return finish(x);
  std::sort(jobs.begin(),jobs.end());
  int lo=-INF,hi=0,min_date=INF;
  for(int j:jobs) {lo=std::max(lo,std::min(a.p[j][first],a.p[j][second])-a.d[j]);hi+=std::max(a.p[j][first],a.p[j][second]);min_date=std::min(min_date,a.d[j]);}
  hi+=a.alpha*(int(jobs.size())-1)-min_date;
  GRBModel model(env);configure(model,local.left());
  auto lambda=model.addVar(lo,hi,1,GRB_INTEGER);
  std::vector<GRBVar> z;GRBLinExpr changed=0;
  for(int j:jobs) {
    z.push_back(model.addVar(0,1,0,GRB_BINARY));z.back().set(GRB_DoubleAttr_Start,x[j]==second);
    if(x[j]==second) changed+=1-z.back();else changed+=z.back();
  }
  model.addConstr(changed>=1);
  std::set<int> dates;for(int j:jobs)dates.insert(a.d[j]);
  for(int k:{first,second}) for(int date:dates) {
    GRBLinExpr processing=0,energy=0;
    auto recharges=model.addVar(0,jobs.size()-1,0,GRB_INTEGER);
    auto price=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);energy+=a.gamma*price;
    for(int h=0;h<int(jobs.size());++h) {
      int j=jobs[h];if(a.d[j]>date)continue;
      GRBLinExpr assign=k==second?GRBLinExpr(z[h]):1-z[h];
      auto excess=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);
      model.addConstr(excess>=a.deviation[j][k]*assign-price);
      processing+=a.p[j][k]*assign;energy+=a.nominal[j][k]*assign+excess;
    }
    model.addConstr(energy<=a.capacity[k]*(1+recharges));
    model.addConstr(processing+a.alpha*recharges<=date+lambda);
  }
  if(local.expired() || budget.expired()) return finish(x);
  model.set(GRB_DoubleParam_TimeLimit,std::min(local.left(),budget.left()));model.optimize();
  if(!model.get(GRB_IntAttr_SolCount)) return finish(x);
  auto result=x;
  for(int h=0;h<int(jobs.size());++h) result[jobs[h]]=z[h].get(GRB_DoubleAttr_X)>.5?second:first;
  for(int j:jobs) stats.pair_changed_jobs+=result[j]!=x[j];
  ++stats.pair_candidates;return finish(result);
}

bool Search::chain_descent(std::vector<int>& x) {
  Clock timing; ++stats.chain_calls;
  auto expired=[&](){return search_expired() || budget.elapsed()>=opt.phase_end;};
  auto finish=[&](bool ok){stats.chain_seconds+=timing.elapsed();return ok;};
  auto incumbent=score(x,!opt.manuscript);remember(x);
  if(opt.target_chains && incumbent.front()<=0)return finish(false);
  auto sets=groups(x);std::vector<int> roots;
  for(int k=0;k<a.m;++k) if(entry(k,sets[k]).upper==incumbent.front())
    roots.insert(roots.end(),sets[k].begin(),sets[k].end());
  std::shuffle(roots.begin(),roots.end(),rng);
  for(int root:roots) {
    if(expired()) break;
    auto partial=x;partial[root]=-1;
    if(opt.monotone_deletion)inherit_removals(x,partial);
    int carry=root,source=x[root];std::vector<bool> visited(a.m,false);visited[source]=true;
    int target=incumbent.front()-1;
    if(opt.target_chains) {
      auto jobs=groups(partial);
      if(entry(source,jobs[source]).lower>target){++stats.chain_root_prunes;continue;}
    }
    std::vector<bool> moved(a.n,false);moved[root]=true;
    for(int depth=0;depth<a.m && !expired();++depth) {
      if(opt.target_chains) {
        incumbent=score(x);target=incumbent.front()-1;
        if(target<0)return finish(false);
        auto jobs=groups(partial);bool impossible=false;
        // A visited AGV cannot release another job later in this chain.
        // Closing the chain at its source can only add a job there.
        for(int k=0;k<a.m;++k)if(visited[k] && entry(k,jobs[k]).lower>target)impossible=true;
        if(impossible){++stats.chain_closed_prunes;break;}
      }
      std::vector<std::pair<std::vector<int>,int>> insertions;
      for(int k=0;k<a.m;++k) if(!visited[k] || (depth>0 && k==source)) {
        if(expired()) break;
        partial[carry]=k;insertions.push_back({score(partial),k});partial[carry]=-1;
      }
      std::sort(insertions.begin(),insertions.end());
      if(insertions.empty() || expired()) break;
      int destination=insertions.front().second;
      partial[carry]=destination;++stats.chain_candidates;
      auto value=insertions.front().first;
      if(opt.target_chains && value.front()>target && opt.adaptive) {
        auto jobs=groups(partial);std::vector<Entry*> candidates;bool possible=true;
        for(int k=0;k<a.m;++k) {
          auto& e=entry(k,jobs[k]);possible&=e.lower<=target;candidates.push_back(&e);
        }
        if(possible) {
          std::sort(candidates.begin(),candidates.end(),[](auto p,auto q){return p->upper>q->upper;});
          for(auto* e:candidates) {
            if(expired())break;
            if(e->upper>target)refine(*e,opt.local_limit,target);
            if(e->lower>target)break;
          }
          value=score(partial);
        }
      } else if(!opt.target_chains && value>=incumbent && opt.adaptive && !cheap_screen(partial,incumbent))
          value=score(partial,true,&incumbent);
      incumbent=score(x);
      if(opt.target_chains?std::max(0,value.front())<std::max(0,incumbent.front()):value<incumbent) {
        stats.chain_primary+=std::max(0,value.front())<std::max(0,incumbent.front());
        x=partial;remember(x);++stats.chain_improvements;++stats.moves;
        stats.chain_jobs+=depth+1;stats.chain_compound+=depth>0;return finish(true);
      }
      if(destination==source || expired()) break;
      visited[destination]=true;
      // A recharge threshold may require freeing a receiving machine before
      // the transfer becomes useful. Select the largest marginal relief on
      // that machine, while leaving the inserted job in place.
      int eject=-1;std::vector<int> best_partial;
      std::vector<int> candidates;
      for(int j=0;j<a.n;++j) if(partial[j]==destination && !moved[j]) candidates.push_back(j);
      std::shuffle(candidates.begin(),candidates.end(),rng);
      for(int j:candidates) {
        if(expired()) break;
        auto parent=partial;partial[j]=-1;
        if(opt.monotone_deletion)inherit_removals(parent,partial);
        if(opt.target_chains) {
          auto jobs=groups(partial);
          if(entry(destination,jobs[destination]).lower>target) {
            ++stats.chain_ejection_prunes;partial[j]=destination;continue;
          }
        }
        auto val=score(partial);partial[j]=destination;
        if(eject<0 || val<best_partial) {eject=j;best_partial=val;}
      }
      if(eject<0 || expired()) break;
      partial[eject]=-1;moved[eject]=true;carry=eject;
    }
  }
  return finish(false);
}
std::vector<int> Search::master_initialization(const std::vector<int>& start,int upper) {
  if(opt.root_mode==0 || budget.left()<.01) return start;
  const bool linear=opt.root_mode==2;
  const bool legacy=opt.root_mode==1;
  Clock timing;GRBModel model(env);
  const double allowed=legacy?std::min(opt.root_limit,budget.left()):budget.left();
  configure(model,allowed);
  if(opt.root_mode==3) {
    model.set(GRB_DoubleParam_NodeLimit,1);
    model.set(GRB_DoubleParam_Heuristics,0);
  }
  auto lambda=model.addVar(global_lower,upper,1,linear?GRB_CONTINUOUS:GRB_INTEGER);
  std::vector<std::vector<GRBVar>> x(a.n);
  for(int j=0;j<a.n;++j) {
    GRBLinExpr sum=0;
    for(int k=0;k<a.m;++k) {
      x[j].push_back(model.addVar(0,1,0,linear?GRB_CONTINUOUS:GRB_BINARY));
      if(legacy)x[j].back().set(GRB_DoubleAttr_Start,start[j]==k);
      sum+=x[j].back();
    }
    model.addConstr(sum==1);
  }
  std::set<int> dates(a.d.begin(),a.d.end());
  for(int k=0;k<a.m;++k) for(int d:dates) {
    GRBLinExpr work=0,energy=0;
    auto charges=model.addVar(0,a.n-1,0,linear?GRB_CONTINUOUS:GRB_INTEGER);
    auto price=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS); energy+=a.gamma*price;
    for(int j=0;j<a.n;++j) if(a.d[j]<=d) {
      auto excess=model.addVar(0,GRB_INFINITY,0,GRB_CONTINUOUS);
      model.addConstr(excess>=a.deviation[j][k]*x[j][k]-price);
      energy+=a.nominal[j][k]*x[j][k]+excess; work+=a.p[j][k]*x[j][k];
    }
    model.addConstr(energy<=a.capacity[k]*(1+charges));
    model.addConstr(work+a.alpha*charges<=lambda+d);
  }
  if(budget.left()<.001) return start;
  model.set(GRB_DoubleParam_TimeLimit,std::max(.001,legacy?
    std::min(opt.root_limit-timing.elapsed(),budget.left()):budget.left()));
  model.optimize(); stats.mip_seconds+=timing.elapsed();
  int status=model.get(GRB_IntAttr_Status);
  if(status==GRB_INFEASIBLE) throw std::runtime_error("root contradicts feasible assignment");
  double bound=model.get(GRB_DoubleAttr_ObjBound);
  if(std::isfinite(bound) && bound>-INF && bound<INF) global_lower=std::max(global_lower,static_cast<int>(std::ceil(bound-1e-6)));
  if(!legacy || !model.get(GRB_IntAttr_SolCount)) return start;
  std::vector<int> next(a.n,-1);
  for(int j=0;j<a.n;++j) for(int k=0;k<a.m;++k) if(x[j][k].get(GRB_DoubleAttr_X)>.5) next[j]=k;
  return next;
}

bool Search::accept_reconstruction(const std::vector<int>& current,
                                   const std::vector<int>& candidate) {
  if(opt.plateau_acceptance) return candidate.front()<=current.front();
  if(opt.greedy_acceptance) return candidate<current;
  if(opt.standard_sa_acceptance) {
    const int deterioration=candidate.front()-current.front();
    if(deterioration<=0)return true;
    const double fraction=std::clamp(budget.elapsed()/opt.limit,0.0,1.0);
    const double initial=std::max(1.0,static_cast<double>(a.alpha))/(-std::log(.5));
    const double final=1.0/(-std::log(.01));
    const double temperature=initial*std::pow(final/initial,fraction);
    return std::generate_canonical<double,32>(rng)<std::exp(-deterioration/temperature);
  }
  const double cooling=opt.cooling_acceptance?1-budget.elapsed()/opt.limit:1;
  const double temperature=std::max(1.0,a.alpha*opt.temperature_scale*cooling)*(opt.optimistic?4:1);
  std::size_t component=0;
  if(opt.lex_acceptance)
    while(component+1<current.size() && current[component]==candidate[component])++component;
  return candidate<current ||
    std::generate_canonical<double,32>(rng)<std::exp((current[component]-candidate[component])/temperature);
}

std::vector<int> Search::recombine(const std::vector<int>& start) {
  if(budget.left()<.01) return start;
  Clock pool_clock;
  auto incumbent_sets=groups(start); auto best_sets=groups(best);
  GRBModel model(env); configure(model,std::min(opt.pool_limit,budget.left()));
  int ub=score(start).front(); auto lambda=model.addVar(global_lower,ub,1,GRB_INTEGER);
  std::vector<GRBLinExpr> coverage(a.n,0);
  struct Column { GRBVar variable; int vehicle; std::vector<int> jobs; };
  std::vector<Column> columns;
  for(int k=0;k<a.m;++k) {
    GRBLinExpr choose=0,value=0;
    int count=0;
    std::vector<const std::pair<const std::vector<int>,std::unique_ptr<Entry>>*> ranked;
    for(const auto& item:cache[k]) ranked.push_back(&item);
    if(opt.managed) std::stable_sort(ranked.begin(),ranked.end(),[](auto p,auto q){return p->second->upper<q->second->upper;});
    for(const auto* ptr:ranked) { const auto& item=*ptr;
      bool essential=item.first==incumbent_sets[k] || item.first==best_sets[k];
      if(!essential && (item.second->upper>ub || count>=opt.max_pool/a.m)) continue;
      auto z=model.addVar(0,1,0,GRB_BINARY); z.set(GRB_DoubleAttr_Start,item.first==incumbent_sets[k]);
      choose+=z; value+=std::max(global_lower,item.second->upper)*z;
      for(int j:item.first) coverage[j]+=z;
      columns.push_back({z,k,item.first}); ++count;
    }
    model.addConstr(choose==1); model.addConstr(lambda>=value);
  }
  for(auto expression:coverage) model.addConstr(expression==1);
  if(budget.left()<.001) return start;
  model.set(GRB_DoubleParam_TimeLimit,std::max(.001,std::min(opt.pool_limit-pool_clock.elapsed(),budget.left())));
  Clock timing; model.optimize(); ++stats.pool_calls; stats.mip_seconds+=pool_clock.elapsed();
  if(!model.get(GRB_IntAttr_SolCount)) return start;
  std::vector<int> result(a.n,-1);
  for(const auto& col:columns) if(col.variable.get(GRB_DoubleAttr_X)>.5)
    for(int j:col.jobs) { if(result[j]>=0) throw std::runtime_error("pool overlap"); result[j]=col.vehicle; }
  if(std::find(result.begin(),result.end(),-1)!=result.end()) throw std::runtime_error("pool incomplete");
  if(score(result)<score(start)) ++stats.pool_improvements;
  return result;
}
Result Search::finish(int initial) {
  Result r; r.initial=initial; r.assignment=best; r.lower=global_lower;
  auto sets=groups(best); r.upper=r.assignment_lower=r.static_upper=-INF;
  for(int k=0;k<a.m;++k) {
    auto& e=entry(k,sets[k]); r.local_lower.push_back(e.lower); r.local_upper.push_back(e.upper);
    r.upper=std::max(r.upper,e.upper); r.assignment_lower=std::max(r.assignment_lower,e.lower);
    r.static_upper=std::max(r.static_upper,e.static_upper);
  }
  if(r.lower>r.upper) throw std::runtime_error("global bounds inconsistent");
  if(opt.progress) opt.progress(budget.elapsed(),r.upper,r.lower,stats,best);
  r.seconds=budget.elapsed(); r.time_best=time_best; r.stats=stats; return r;
}
bool Search::search_expired() const {
  if(budget.expired()) return true;
  if(opt.stop_zero && !best.empty() && best_score.front()<=0) return true;
  if(opt.manuscript && !best.empty() && std::max(0,best_score.front())<=std::max(0,global_lower)) return true;
  if(!opt.polishing || best.empty() || budget.left()>(opt.polish_reserve>0?opt.polish_reserve:std::min(2.0,opt.limit*.1))) return false;
  auto sets=groups(best); int lower=-INF,upper=-INF;
  for(int k=0;k<a.m;++k) {
    const auto& e=*cache[k].at(sets[k]); lower=std::max(lower,e.lower); upper=std::max(upper,e.upper);
  }
  return lower<upper;
}
void Search::polish_incumbent() {
  if(opt.stop_zero && !best.empty() && best_score.front()<=0) return;
  auto sets=groups(best);
  while(!budget.expired()) {
    int lower=-INF,upper=-INF; Entry* critical=nullptr;
    for(int k=0;k<a.m;++k) {
      auto& e=entry(k,sets[k]); lower=std::max(lower,e.lower);
      if(e.upper>upper) { upper=e.upper; critical=&e; }
    }
    if(lower==upper || !critical || critical->data.n()>62) break;
    int old_lower=critical->lower,old_upper=critical->upper;
    refine(*critical,budget.left()); remember(best);
    if(old_lower==critical->lower && old_upper==critical->upper) break;
  }
}
Result Search::run_manuscript() {
  // Initialize best before local search so every selective check can target
  // strict improvement of the best HIG bound. No separate refinement phases.
  auto x=construct(0);int initial=score(x).front();remember(x);
  master_initialization(x,initial);remember(x);
  if(!search_expired())descend(x);
  remember(x);
  while(a.m>=2 && !search_expired()) {
    auto next=block_shake(x);
    if(!search_expired())descend(next);
    remember(next);
    // Verification can also tighten cached bounds shared with x and best.
    if(accept_reconstruction(score(x),score(next)))x=std::move(next);
    remember(x);
  }
  return finish(initial);
}
Result Search::run() {
  if(opt.manuscript)return run_manuscript();
  bool adaptive=opt.adaptive;
  if(opt.staged) opt.adaptive=false;
  auto x=construct(0); int initial=score(x).front(); remember(x);
  if(opt.stop_zero && best_score.front()<=0) return finish(initial);
  if(opt.staged) { opt.adaptive=false; opt.phase_end=opt.event_warm?opt.limit:std::max(budget.elapsed(),opt.limit*opt.warm_fraction); }
  auto seed=master_initialization(x,initial);score(seed,opt.adaptive);remember(seed);
  if(opt.stop_zero && best_score.front()<=0) return finish(initial);
  if(opt.cyclic) {
    opt.adaptive=adaptive; score(seed,true); remember(seed); opt.adaptive=false;
  }
  if(score(seed)<score(x)) x=seed;
  descend(x); remember(x);
  if(opt.stop_zero && best_score.front()<=0) return finish(initial);
  if(opt.staged) {
    if(opt.event_warm) opt.phase_end=budget.elapsed();
    int warm=0;
    while(!search_expired() && budget.elapsed()<opt.phase_end && score(best).front()>global_lower) {
      auto old_value=guide(x);auto next=shake(warm%3==0?best:x,2+warm%5);
      if(opt.post_perturbation_descent)descend(next);
      remember(next);auto new_value=guide(next);
      if(accept_reconstruction(old_value,new_value))x=next;
      ++warm;
    }
    if(opt.probe_warm && !opt.event_warm && !search_expired() && score(best).front()>global_lower) {
      int old=score(best).front();opt.adaptive=true;double local_limit=opt.local_limit;opt.local_limit=.1;
      score(best,true);remember(best);opt.local_limit=local_limit;
      if(old-score(best).front()<std::max(2,a.alpha/5)) {
        opt.adaptive=false;opt.phase_end=std::min(10.0,opt.limit*.25);
        while(!search_expired() && budget.elapsed()<opt.phase_end && score(best).front()>global_lower) {
          auto old_value=guide(x);auto next=shake(warm%3==0?best:x,2+warm%5);
          if(opt.post_perturbation_descent)descend(next);
          remember(next);auto new_value=guide(next);
          if(accept_reconstruction(old_value,new_value))x=next;
          ++warm;
        }
      }
    }
    opt.adaptive=adaptive; opt.phase_end=std::min(opt.limit,budget.elapsed()+std::max(.05,opt.limit*.12));
    if(opt.stop_zero && best_score.front()<=0) return finish(initial);
    if(opt.alternating) opt.phase_end=std::min(opt.phase_end,opt.limit*(opt.warm_fraction+.1));
    x=best; score(x,true); remember(x);
    if(opt.cyclic && guide(seed)<guide(x)) x=seed;
    descend(x);
  }
  int iteration=0;
  do {
    while(opt.iterated && !search_expired() && score(best).front()>global_lower) {
      if(opt.managed && budget.elapsed()-last_pulse>=opt.pulse_interval) { pulse(opt.pulse_limit); last_pulse=budget.elapsed(); }
      if(opt.restart && budget.elapsed()-last_restart>=5 && budget.elapsed()-time_best>=2) {
        x=construct(1+stats.restarts%2); ++stats.restarts; last_restart=budget.elapsed();
        score(x,true); remember(x);
      }
      if(opt.cyclic) opt.lower_weight=std::vector<int>{4,3,0}[(iteration/4)%3];
      if(opt.recombination && (opt.recombine_interval>0?budget.elapsed()-last_pool>=opt.recombine_interval:iteration%3==0)) {
        auto mix=recombine(best); last_pool=budget.elapsed(); score(mix,true); remember(mix); if(score(mix)<score(x)) x=mix;
      }
      if(budget.expired()) break;
      if(opt.staged) opt.phase_end=std::min(opt.limit,budget.elapsed()+std::max(.05,opt.managed?std::min(1.0,opt.limit*.08):opt.limit*.08));
      if(opt.alternating) {
        double width=std::max(.05,opt.limit*.1),start=opt.limit*opt.warm_fraction;
        int phase=std::max(0,static_cast<int>((budget.elapsed()-start)/width));
        opt.adaptive=phase%2==0;
        opt.phase_end=std::min(opt.phase_end,start+(phase+1)*width);
      }
      auto next=shake(iteration%3==0?best:x,2+iteration%std::max(2,std::min(8,a.n/4)));
      if(opt.post_perturbation_descent)descend(next);
      remember(next);
      auto old_value=guide(x),new_value=guide(next);
      if(accept_reconstruction(old_value,new_value))x=next;
      ++iteration;
    }
    if(opt.polishing) polish_incumbent();
    // If certification finishes early, reuse the remaining time for search.
    // The certified incumbent is retained throughout any further neighbourhoods.
  } while(opt.polishing && !search_expired() && score(best).front()>global_lower);
  return finish(initial);
}
bool Search::cheap_screen(const std::vector<int>& x,const std::vector<int>& incumbent) {
  Clock screen_clock;
  auto sets=groups(x); std::vector<int> lower;
  for(int k=0;k<a.m;++k) {
    auto it=cache[k].find(sets[k]);
    if(it==cache[k].end()) {
      Clock lc;auto saved=lower_cache[k].find(sets[k]);
      if(opt.memo_lower && saved!=lower_cache[k].end()){lower.push_back(saved->second);++stats.lower_hits;}
      else {
        auto local=a.local(k,x);
        lower.push_back(opt.incremental_lower?incremental_lower(local):workload_lower(local,opt.packing));
        if(opt.memo_lower){lower_cache[k].emplace(sets[k],lower.back());++stats.lower_entries;}
      }
      stats.lower_seconds+=lc.elapsed();
    } else lower.push_back(it->second->lower);
    if(lower.back()>incumbent.front()) { ++stats.fast_screened;stats.screen_seconds+=screen_clock.elapsed(); return true; }
  }
  stats.screen_seconds+=screen_clock.elapsed();
  if(ranked(lower)>=incumbent) { ++stats.fast_screened; return true; }
  return false;
}
void Search::pulse(double seconds) {
  Clock timing(std::min(seconds,budget.left())); ++stats.pulses;
  int before=score(best).front();
  std::vector<std::vector<int>> candidates{best};
  candidates.insert(candidates.end(),elite.begin(),elite.end()); elite.clear();
  for(auto& x:candidates) {
    if(timing.expired() || budget.expired()) break;
    auto sets=groups(x); int lower=-INF;
    std::vector<Entry*> entries;
    for(int k=0;k<a.m;++k) { auto& e=entry(k,sets[k]); lower=std::max(lower,e.lower); entries.push_back(&e); }
    int target=score(best).front()-1; if(lower>target) continue;
    std::sort(entries.begin(),entries.end(),[](auto p,auto q){return p->upper>q->upper;});
    for(auto* e:entries) {
      if(timing.expired() || budget.expired() || e->lower>target) break;
      if(e->upper>target) refine(*e,std::min(timing.left(),std::max(.05,seconds/2)),target);
    }
    remember(x);
  }
  stats.pulse_improvements+=std::max(0,before-score(best).front()); stats.pulse_seconds+=timing.elapsed();
}
Result Search::audit(const std::vector<int>& x) {
  if(x.size()!=static_cast<std::size_t>(a.n)) throw std::runtime_error("assignment length");
  for(int k:x) if(k<0 || k>=a.m) throw std::runtime_error("invalid vehicle");
  best=x; auto sets=groups(x); int initial=score(x).front();
  for(int pass=0;pass<3 && !budget.expired();++pass) {
    std::vector<Entry*> entries;
    for(int k=0;k<a.m;++k) entries.push_back(&entry(k,sets[k]));
    std::sort(entries.begin(),entries.end(),[](Entry* p,Entry* q){return p->upper>q->upper;});
    for(auto* e:entries) if(e->lower<e->upper) refine(*e,budget.left()/std::max(1,a.m-pass));
  }
  return finish(initial);
}
Result Search::probe_chain(std::vector<int> x) {
  if(x.size()!=std::size_t(a.n))throw std::runtime_error("probe assignment length");
  for(int k:x)if(k<0 || k>=a.m)throw std::runtime_error("probe assignment machine");
  int initial=score(x).front();remember(x);chain_descent(x);remember(x);return finish(initial);
}
Result Search::probe_shake(std::vector<int> x) {
  if(x.size()!=std::size_t(a.n))throw std::runtime_error("probe assignment length");
  for(int k:x)if(k<0 || k>=a.m)throw std::runtime_error("probe assignment machine");
  int initial=score(x).front();remember(x);auto changed=shake(x,3);score(changed,true);
  for(int k:changed)if(k<0 || k>=a.m)throw std::runtime_error("invalid reconstructed assignment");
  // Return the actual reconstruction for direct correctness checks, even
  // when diversification deliberately worsens the incumbent.
  best=changed;best_score=score(changed);return finish(initial);
}
}
