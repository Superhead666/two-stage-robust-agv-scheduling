// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace agv {
using Mask = std::uint64_t;
constexpr double ENERGY_EPS = 1e-9;
constexpr int INF = 1000000000;
struct Clock {
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  double limit;
  explicit Clock(double seconds = 1e9) : limit(seconds) {}
  double elapsed() const { return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); }
  double left() const { return std::max(0.0, limit-elapsed()); }
  bool expired() const { return elapsed() >= limit; }
};
struct Local {
  std::vector<int> p, d, ids;
  std::vector<double> e, nominal, deviation;
  double capacity = 10, gamma = 1;
  int alpha = 10;
  int n() const { return static_cast<int>(p.size()); }
};
struct Instance {
  int n = 0, m = 0, alpha = 20;
  double gamma = 2;
  std::vector<int> d;
  std::vector<double> capacity;
  std::vector<std::vector<int>> p;
  std::vector<std::vector<double>> nominal, deviation;
  std::string id;
  Local local(int k, const std::vector<int>& assignment) const;
};
struct Schedule {
  std::vector<int> order;
  std::vector<bool> recharge;
  std::vector<Mask> blocks;
  int lateness = -INF;
};
enum class Status { Feasible, Infeasible, Unknown };
struct DetResult {
  Status status = Status::Unknown;
  Schedule schedule;
  std::uint64_t nodes = 0, dominated = 0, pruned = 0;
  double seconds = 0;
};
struct DetOptions {
  bool dense = false;
  bool bounds = true;
  bool greedy_recharge = true;
  bool cross_dominance = true;
  bool heuristic = true;
};
Schedule evaluate(const Local&, const std::vector<int>&, const std::vector<bool>&);
Schedule edd(const Local&);
bool valid_schedule(const Local&, const Schedule&, int lambda);
DetResult deterministic(const Local&, int lambda, const DetOptions&, Clock&);
DetResult dense_dp(const Local&, int lambda, Clock&);
int deterministic_optimum(const Local&, const DetOptions&, Clock&);
int brute_optimum(const Local&);
std::vector<Schedule> enumerate_partitions(const Local&);
double worst_block_energy(const Local&, Mask);
Local realize(Local, const std::vector<double>&);
Instance generate(int n, int m, unsigned seed, double gamma=2, double variation=.6, double capacity_scale=3, double due_factor=1);
void save_instance(const Instance&, const std::string&);
Instance load_instance(const std::string&);
void validate_instance(const Instance&);
std::string status_name(Status);
}
