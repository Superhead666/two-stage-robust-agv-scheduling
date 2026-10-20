// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "core.hpp"
#include <fstream>
#include <iomanip>

namespace agv {
Local Instance::local(int k, const std::vector<int>& assignment) const {
  Local a; a.capacity=capacity[k]; a.gamma=gamma; a.alpha=alpha;
  for(int j=0;j<n;++j) if(assignment[j]==k) {
    a.ids.push_back(j);a.p.push_back(p[j][k]);a.d.push_back(d[j]);
    a.nominal.push_back(nominal[j][k]);a.deviation.push_back(deviation[j][k]);a.e.push_back(nominal[j][k]);
  } return a;
}
std::string status_name(Status s) { return s==Status::Feasible?"feasible":s==Status::Infeasible?"infeasible":"unknown"; }
Schedule evaluate(const Local& a,const std::vector<int>& order,const std::vector<bool>& cuts) {
  Schedule s; s.order=order;s.recharge=cuts; int t=0; double b=a.capacity;
  Mask seen=0,block=0;
  if(order.size()!=static_cast<std::size_t>(a.n()) || cuts.size()!=order.size()) throw std::runtime_error("schedule length");
  for(std::size_t q=0;q<order.size();++q) {
    int j=order[q]; if(j<0 || j>=a.n() || (seen&(Mask(1)<<j))) throw std::runtime_error("schedule permutation");
    seen|=Mask(1)<<j;
    if(cuts[q]) { if(q==0) throw std::runtime_error("initial recharge"); s.blocks.push_back(block);block=0;t+=a.alpha;b=a.capacity; }
    if(a.e[j]>b+ENERGY_EPS) { s.lateness=INF;return s; }
    b-=a.e[j];t+=a.p[j];block|=Mask(1)<<j;s.lateness=std::max(s.lateness,t-a.d[j]);
  }
  if(block) s.blocks.push_back(block);
  return s;
}
Schedule edd(const Local& a) {
  std::vector<int> order(a.n());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  std::vector<bool> cuts(a.n());double b=a.capacity;
  for(int q=0;q<a.n();++q) {int j=order[q];cuts[q]=a.e[j]>b+ENERGY_EPS;if(cuts[q])b=a.capacity;b-=a.e[j];}
  return evaluate(a,order,cuts);
}
bool valid_schedule(const Local& a,const Schedule& s,int lambda) {return evaluate(a,s.order,s.recharge).lateness<=lambda;}

DetResult dense_dp(const Local& a,int lambda,Clock& budget) {
  Clock timing;DetResult out;int n=a.n(); if(n==0){out.status=Status::Feasible;return out;}
  if(n>22){out.status=Status::Unknown;return out;}
  std::size_t states=std::size_t(1)<<n;
  if(states*static_cast<std::size_t>(n)>100000000){out.status=Status::Unknown;return out;}
  std::vector<double> rem(states*n,-1);std::vector<std::uint8_t> pred(states*n,255);
  std::vector<int> ps(states);rem[0]=a.capacity;
  for(Mask mask=1;mask<states;++mask){int j=__builtin_ctzll(mask);ps[mask]=ps[mask&(mask-1)]+a.p[j];}
  for(Mask mask=0;mask<states;++mask) {
    if((mask&1023)==0 && budget.expired()){out.seconds=timing.elapsed();return out;}
    for(int r=0;r<n;++r) {double b=rem[mask*n+r];if(b<0)continue;++out.nodes;
      int t=ps[mask]+r*a.alpha;
      for(int j=0;j<n;++j)if(!(mask&(Mask(1)<<j))) {
        Mask next=mask|(Mask(1)<<j);
        if(b+ENERGY_EPS>=a.e[j] && t+a.p[j]<=a.d[j]+lambda && b-a.e[j]>rem[next*n+r]+ENERGY_EPS) {
          rem[next*n+r]=std::max(0.0,b-a.e[j]);pred[next*n+r]=static_cast<std::uint8_t>(j);
        }
        if(mask && r+1<n && a.e[j]<=a.capacity+ENERGY_EPS && t+a.alpha+a.p[j]<=a.d[j]+lambda && a.capacity-a.e[j]>rem[next*n+r+1]+ENERGY_EPS) {
          rem[next*n+r+1]=std::max(0.0,a.capacity-a.e[j]);pred[next*n+r+1]=static_cast<std::uint8_t>(j+64);
        }
      }
    }
  }
  int r=0;while(r<n && rem[(states-1)*n+r]<0)++r;
  out.status=r<n?Status::Feasible:Status::Infeasible;
  if(r<n){Mask mask=states-1;std::vector<int> order;std::vector<bool> cuts;
    while(mask){int code=pred[mask*n+r];if(code==255)throw std::runtime_error("DP reconstruction");int j=code%64;bool c=code>=64;order.push_back(j);cuts.push_back(c);mask^=Mask(1)<<j;r-=c;}
    std::reverse(order.begin(),order.end());std::reverse(cuts.begin(),cuts.end());out.schedule=evaluate(a,order,cuts);
    if(out.schedule.lateness>lambda)throw std::runtime_error("invalid dense witness");
  }
  out.seconds=timing.elapsed();return out;
}

DetResult deterministic(const Local& a,int lambda,const DetOptions& opt,Clock& budget) {
  if(opt.dense)return dense_dp(a,lambda,budget);
  Clock timing;DetResult out;int n=a.n();if(n==0){out.status=Status::Feasible;return out;}if(n>62)return out;
  for(int j=0;j<n;++j)if(a.e[j]>a.capacity+ENERGY_EPS || a.p[j]>a.d[j]+lambda){out.status=Status::Infeasible;return out;}
  if(opt.heuristic){Schedule s=edd(a);if(s.lateness<=lambda){out.status=Status::Feasible;out.schedule=s;out.seconds=timing.elapsed();return out;}}
  std::vector<int> priority(n);std::iota(priority.begin(),priority.end(),0);
  std::stable_sort(priority.begin(),priority.end(),[&](int i,int j){return a.d[i]<a.d[j];});
  Mask full=(Mask(1)<<n)-1;
  std::unordered_map<Mask,std::vector<double>> labels;
  std::vector<int> order;std::vector<bool> cuts;bool interrupted=false;
  std::function<bool(Mask,int,int,double)> dfs=[&](Mask mask,int r,int t,double battery){
    if((++out.nodes&1023)==0 && budget.expired()){interrupted=true;return false;}
    if(mask==full){out.schedule=evaluate(a,order,cuts);return true;}
    auto& lab=labels[mask];if(lab.empty())lab.assign(n,-1);
    for(int q=opt.cross_dominance?0:r;q<=r;++q)if(lab[q]>=battery-ENERGY_EPS){++out.dominated;return false;}
    lab[r]=std::max(lab[r],battery);
    if(opt.bounds){
      int work=0;double energy=0;
      for(int j:priority)if(!(mask&(Mask(1)<<j))){
        work+=a.p[j];energy+=a.e[j];
        int charges=static_cast<int>(std::max(0.0,std::ceil((energy-battery-ENERGY_EPS)/a.capacity)));
        if(t+work+charges*a.alpha>a.d[j]+lambda){++out.pruned;return false;}
      }
    }
    for(int j:priority)if(!(mask&(Mask(1)<<j))){
      bool must=a.e[j]>battery+ENERGY_EPS;
      for(int c=must?1:0;c<=1;++c){
        if(c && (!mask || (!must && opt.greedy_recharge)))continue;
        int finish=t+a.p[j]+c*a.alpha;if(finish>a.d[j]+lambda)continue;
        order.push_back(j);cuts.push_back(c);
        if(dfs(mask|(Mask(1)<<j),r+c,finish,(c?a.capacity:battery)-a.e[j]))return true;
        order.pop_back();cuts.pop_back();if(interrupted)return false;
      }
    } return false;
  };
  bool yes=dfs(0,0,0,a.capacity);out.status=yes?Status::Feasible:interrupted?Status::Unknown:Status::Infeasible;
  out.seconds=timing.elapsed();if(yes && !valid_schedule(a,out.schedule,lambda))throw std::runtime_error("invalid sparse witness");return out;
}
int deterministic_optimum(const Local& a,const DetOptions& opt,Clock& budget){
  if(a.n()==0)return -INF;int lo=-INF,hi=edd(a).lateness;
  for(int j=0;j<a.n();++j)lo=std::max(lo,a.p[j]-a.d[j]);
  while(lo<hi){int mid=lo+(hi-lo)/2;auto r=deterministic(a,mid,opt,budget);if(r.status==Status::Unknown)return INF;if(r.status==Status::Feasible)hi=mid;else lo=mid+1;}return lo;
}
int brute_optimum(const Local& a){
  int n=a.n();if(!n)return -INF;if(n>9)throw std::runtime_error("brute size");
  std::vector<int> order(n);std::iota(order.begin(),order.end(),0);int best=INF;
  do{for(Mask bits=0;bits<(Mask(1)<<(n-1));++bits){std::vector<bool> cuts(n);for(int j=1;j<n;++j)cuts[j]=(bits>>(j-1))&1;best=std::min(best,evaluate(a,order,cuts).lateness);}}while(std::next_permutation(order.begin(),order.end()));return best;
}
std::vector<Schedule> enumerate_partitions(const Local& a){
  if(a.n()>10)throw std::runtime_error("partition enumeration size");std::vector<Schedule> out;std::vector<std::vector<int>> blocks;
  std::function<void(int)> rec=[&](int j){if(j==a.n()){
    struct Block{std::vector<int> order;int p=0,l=-INF;};std::vector<Block> bs;
    for(auto v:blocks){std::stable_sort(v.begin(),v.end(),[&](int i,int k){return a.d[i]<a.d[k];});Block b;b.order=v;for(int k:v){b.p+=a.p[k];b.l=std::max(b.l,b.p-a.d[k]);}bs.push_back(b);}
    std::stable_sort(bs.begin(),bs.end(),[](const Block& x,const Block& y){return x.p-x.l<y.p-y.l;});
    Local zero=a;std::fill(zero.e.begin(),zero.e.end(),0);std::vector<int> order;std::vector<bool> cuts;
    for(const auto& b:bs)for(std::size_t q=0;q<b.order.size();++q){cuts.push_back(!order.empty() && q==0);order.push_back(b.order[q]);}
    out.push_back(evaluate(zero,order,cuts));return;}
    std::size_t count=blocks.size();for(std::size_t q=0;q<count;++q){blocks[q].push_back(j);rec(j+1);blocks[q].pop_back();}
    blocks.push_back({j});rec(j+1);blocks.pop_back();};rec(0);return out;
}
double worst_block_energy(const Local& a,Mask block){double total=0;std::vector<double> ds;for(int j=0;j<a.n();++j)if(block&(Mask(1)<<j)){total+=a.nominal[j];ds.push_back(a.deviation[j]);}std::sort(ds.rbegin(),ds.rend());double left=a.gamma;for(double v:ds){double f=std::min(1.0,std::max(0.0,left));total+=f*v;left-=f;}return total;}
Local realize(Local a,const std::vector<double>& xi){for(int j=0;j<a.n();++j)a.e[j]=a.nominal[j]+a.deviation[j]*xi[j];return a;}
Instance generate(int n,int m,unsigned seed,double gamma,double variation,double scale,double due_factor){
  if(n<1||n>62||m<1||m>32||!std::isfinite(gamma)||!std::isfinite(variation)||!std::isfinite(scale)||!std::isfinite(due_factor)||gamma<0||gamma>n||variation<0||scale<=0||due_factor<0)throw std::runtime_error("invalid generator parameters (1 <= n <= 62)");
  std::mt19937 rng(seed);std::uniform_int_distribution<int> proc(1,20),eng(3,15);std::uniform_real_distribution<double> mult(.8,1.2),due(.25,1.1);
  Instance a;a.n=n;a.m=m;a.gamma=gamma;a.alpha=20;a.id="n"+std::to_string(n)+"m"+std::to_string(m)+"s"+std::to_string(seed);
  a.p.assign(n,std::vector<int>(m));a.nominal.assign(n,std::vector<double>(m));a.deviation=a.nominal;a.capacity.assign(m,0);a.d.resize(n);
  for(int j=0;j<n;++j)for(int k=0;k<m;++k){a.p[j][k]=proc(rng);a.nominal[j][k]=eng(rng);a.deviation[j][k]=std::round(10*a.nominal[j][k]*variation*mult(rng))/10;}
  for(int k=0;k<m;++k){double sum=0,mx=0;for(int j=0;j<n;++j){double e=a.nominal[j][k]+a.deviation[j][k];sum+=e;mx=std::max(mx,e);}a.capacity[k]=std::max(mx,std::round(scale*sum/n));}
  double horizon=(10.5+20/scale)*n/m;for(int j=0;j<n;++j)a.d[j]=static_cast<int>(std::round(due_factor*horizon*due(rng)));
  validate_instance(a);return a;
}
void save_instance(const Instance& a,const std::string& path){std::ofstream f(path);if(!f)throw std::runtime_error("cannot write instance");f<<std::setprecision(17)<<a.n<<' '<<a.m<<' '<<a.alpha<<' '<<a.gamma<<'\n';for(double v:a.capacity)f<<v<<' ';f<<'\n';for(int j=0;j<a.n;++j){f<<a.d[j];for(int k=0;k<a.m;++k)f<<' '<<a.p[j][k]<<' '<<a.nominal[j][k]<<' '<<a.deviation[j][k];f<<'\n';}}
Instance load_instance(const std::string& path){std::ifstream f(path);Instance a;a.id=path;if(!(f>>a.n>>a.m>>a.alpha>>a.gamma)||a.n<1||a.n>62||a.m<1||a.m>32)throw std::runtime_error("invalid instance header");a.capacity.resize(a.m);for(auto& v:a.capacity)f>>v;a.d.resize(a.n);a.p.assign(a.n,std::vector<int>(a.m));a.nominal.assign(a.n,std::vector<double>(a.m));a.deviation=a.nominal;for(int j=0;j<a.n;++j){f>>a.d[j];for(int k=0;k<a.m;++k)f>>a.p[j][k]>>a.nominal[j][k]>>a.deviation[j][k];}if(!f)throw std::runtime_error("incomplete instance");validate_instance(a);return a;}
void validate_instance(const Instance& a){
  if(a.n<1||a.n>62||a.m<1||a.m>32||a.alpha<0||a.alpha>1000000||!std::isfinite(a.gamma)||a.gamma<0||a.gamma>a.n)throw std::runtime_error("unsupported instance dimensions or budget");
  for(double cap:a.capacity)if(!std::isfinite(cap)||cap<=0)throw std::runtime_error("invalid battery capacity");
  for(int j=0;j<a.n;++j){if(a.d[j]<-1000000||a.d[j]>1000000)throw std::runtime_error("due date range");
    for(int k=0;k<a.m;++k)if(a.p[j][k]<=0||a.p[j][k]>1000000||!std::isfinite(a.nominal[j][k])||!std::isfinite(a.deviation[j][k])||a.nominal[j][k]<0||a.deviation[j][k]<0||a.nominal[j][k]+a.deviation[j][k]>a.capacity[k]+ENERGY_EPS)throw std::runtime_error("job must fit a full battery even at its upper energy bound");}
}
}
