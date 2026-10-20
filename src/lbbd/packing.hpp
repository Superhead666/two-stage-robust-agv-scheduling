// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "core.hpp"
#include "dff.hpp"
#include <map>
namespace lb {
using namespace agv;
struct Promotion {
  double value=0;
  std::vector<double> xi;
};
// Minimum budget for each attainable integer increase in q*f_q. Actual
// witnesses use a positive energy margin; boundary-ambiguous cuts are omitted.
inline Promotion promote(const Instance& a,int k,const std::vector<int>& jobs,int q) {
  std::vector<std::vector<std::pair<int,double>>> choices;
  int width=0;
  for(int j:jobs) {
    const double nom=a.nominal[j][k],dev=a.deviation[j][k],cap=a.capacity[k];
    const int low=int(std::llround(q*transformed_size(nom/cap,q)));
    const int high=int(std::llround(q*transformed_size((nom+dev)/cap,q)));
    choices.push_back({{0,0}});
    if(dev>0)for(int level=low+1;level<=high;++level) {
      const double threshold=cap*level/(q+1.0),margin=1e-8*std::max(1.0,cap);
      const double need=std::max(0.0,(threshold+margin-nom)/dev);
      if(need<=1 && need<=a.gamma && q*transformed_size((nom+dev*need)/cap,q)>=level-1e-7)
        choices.back().push_back({level-low,need});
    }
    width+=choices.back().back().first;
  }
  std::vector<double> dp(width+1,1e100);dp[0]=0;
  std::vector<std::vector<int>> previous(jobs.size(),std::vector<int>(width+1,-1));
  std::vector<std::vector<int>> chosen(jobs.size(),std::vector<int>(width+1,-1));
  for(std::size_t t=0;t<jobs.size();++t) {
    std::vector<double> next(width+1,1e100);
    for(int gain=0;gain<=width;++gain)if(dp[gain]<=a.gamma+1e-12)
      for(int c=0;c<int(choices[t].size());++c) {
        const auto [delta,cost]=choices[t][c];
        if(gain+delta<=width && dp[gain]+cost<next[gain+delta]) {
          next[gain+delta]=dp[gain]+cost;previous[t][gain+delta]=gain;chosen[t][gain+delta]=c;
        }
      }
    dp=std::move(next);
  }
  int best=0;for(int v=0;v<=width;++v)if(dp[v]<=a.gamma+1e-12)best=v;
  Promotion out;out.xi.assign(a.n,0);int gain=best;
  for(int t=int(jobs.size())-1;t>=0;--t) {
    const int c=chosen[t][gain];if(c<0)throw std::runtime_error("promotion reconstruction");
    out.xi[jobs[t]]=choices[t][c].second;gain=previous[t][gain];
  }
  const double sum=std::accumulate(out.xi.begin(),out.xi.end(),0.0);
  if(sum>a.gamma+1e-9)throw std::runtime_error("invalid promotion budget");
  // Round down any floating-point budget excess, then evaluate the actual witness.
  if(sum>a.gamma)for(double& v:out.xi)v*=a.gamma/sum;
  for(int j:jobs)out.value+=transformed_size((a.nominal[j][k]+a.deviation[j][k]*out.xi[j])/a.capacity[k],q);
  return out;
}
struct PackResult {int q=0;double value=0;std::vector<double> xi;bool complete=true;};
inline PackResult packing_bound(const Instance& a,int k,const std::vector<int>& jobs,Clock& clock,bool scale=false) {
  PackResult best;best.xi.assign(a.n,0);
  int maximum=a.n-1;
  if(scale) {
    int divisor=0;bool integral=true;
    for(int j=0;j<a.n;++j){double v=a.nominal[j][k];if(v>1e9||std::abs(v-std::round(v))>1e-9)integral=false;else divisor=std::gcd(divisor,int(std::llround(v)));}
    // For integral nominal energies, the capacity-to-grid ratio gives the
    // finite resolution required by this DFF family. For general real data,
    // retain the standard q=1,...,n-1 family without an arbitrary cutoff.
    if(integral&&divisor>0) {
      const double resolution=std::floor(a.capacity[k]/divisor);
      if(resolution<std::numeric_limits<int>::max())maximum=std::max(maximum,int(resolution));
    }
  }
  for(int q=1;q<=maximum;++q) {
    if(clock.expired()){best.complete=false;break;}
    double upper=0;
    for(int j:jobs)upper+=transformed_size((a.nominal[j][k]+a.deviation[j][k])/a.capacity[k],q);
    if(upper<=best.value+1e-8)continue;
    auto z=promote(a,k,jobs,q);
    if(z.value>best.value+1e-8){best.q=q;best.value=z.value;best.xi=std::move(z.xi);}
  }
  return best;
}
}
