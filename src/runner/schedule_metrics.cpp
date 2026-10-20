// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace agv;

namespace {
struct Solution {
  int n = 0;
  int upper = 0;
  int lower = 0;
  std::vector<int> assignment;
};

Solution read_solution(const std::string& path, int machines) {
  std::ifstream input(path);
  Solution solution;
  if (!(input >> solution.n >> solution.upper >> solution.lower) ||
      solution.n < 1 || solution.lower < 0 || solution.upper < solution.lower)
    throw std::runtime_error("invalid solution file");
  solution.assignment.resize(solution.n);
  for (int& machine : solution.assignment)
    if (!(input >> machine) || machine < 0 || machine >= machines)
      throw std::runtime_error("invalid assignment in solution file");
  return solution;
}

void write_unknown(const std::string& output, double seconds) {
  std::ofstream file(output);
  if (!file) throw std::runtime_error("cannot write metrics file");
  file << "status,nominal_maximum_tardiness,nominal_maximum_completion_time,"
          "nominal_mean_completion_time,nominal_total_recharges,"
          "nominal_maximum_recharges_per_agv,metrics_seconds\n";
  file << std::setprecision(12) << "unknown,,,,,," << seconds << '\n';
}
}

int main(int argc, char** argv) {
  const auto start = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
  };
  try {
    if (argc != 5)
      throw std::runtime_error(
          "agv_schedule_metrics INSTANCE SOLUTION OUTPUT TIME_LIMIT");
    const Instance instance = load_instance(argv[1]);
    const Solution solution = read_solution(argv[2], instance.m);
    if (solution.n != instance.n) throw std::runtime_error("solution size mismatch");
    std::size_t used = 0;
    const double limit = std::stod(argv[4], &used);
    if (used != std::string(argv[4]).size() || !std::isfinite(limit) || limit <= 0)
      throw std::runtime_error("time limit must be positive and finite");

    Clock budget(limit);
    long long completion_sum = 0;
    int maximum_completion = 0;
    int maximum_tardiness = 0;
    int total_recharges = 0;
    int maximum_recharges = 0;
    DetOptions options;

    for (int k = 0; k < instance.m; ++k) {
      Local local = instance.local(k, solution.assignment);
      if (!local.n()) continue;
      local.e = local.nominal;

      int low = 0;
      int high = std::max(0, edd(local).lateness);
      while (low < high && !budget.expired()) {
        const int middle = low + (high-low)/2;
        const DetResult result = deterministic(local, middle, options, budget);
        if (result.status == Status::Unknown) {
          write_unknown(argv[3], elapsed());
          return 0;
        }
        if (result.status == Status::Feasible) high = middle;
        else low = middle+1;
      }
      if (budget.expired()) {
        write_unknown(argv[3], elapsed());
        return 0;
      }
      const DetResult result = deterministic(local, low, options, budget);
      if (result.status != Status::Feasible) {
        write_unknown(argv[3], elapsed());
        return 0;
      }

      int time = 0;
      int recharges = 0;
      for (std::size_t position = 0; position < result.schedule.order.size(); ++position) {
        if (result.schedule.recharge[position]) {
          time += local.alpha;
          ++recharges;
        }
        time += local.p[result.schedule.order[position]];
        completion_sum += time;
        maximum_completion = std::max(maximum_completion, time);
      }
      maximum_tardiness = std::max(maximum_tardiness, low);
      total_recharges += recharges;
      maximum_recharges = std::max(maximum_recharges, recharges);
    }

    std::ofstream file(argv[3]);
    if (!file) throw std::runtime_error("cannot write metrics file");
    file << "status,nominal_maximum_tardiness,nominal_maximum_completion_time,"
            "nominal_mean_completion_time,nominal_total_recharges,"
            "nominal_maximum_recharges_per_agv,metrics_seconds\n";
    file << std::setprecision(12) << "complete," << maximum_tardiness << ','
         << maximum_completion << ','
         << static_cast<double>(completion_sum)/instance.n << ','
         << total_recharges << ',' << maximum_recharges << ',' << elapsed() << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
