// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "robust.hpp"
#include <filesystem>
#include <map>

namespace baseline {
using namespace agv;
namespace fs = std::filesystem;
inline constexpr const char* BASELINE_POLICY = "ejor_fair_ccg_v1";

struct Options {
  std::string formulation = "block";
  double limit = 1200;
  unsigned seed = 1;
  int threads = 1;
  int master_target_upper = INF;
  int initial_lower = 0;
  fs::path warm_solution;
  std::vector<int> warm_assignment;
  fs::path trace, solution;
};
struct Result {
  std::string status = "unresolved", reason, certificate = "none";
  int lower = 0, upper = INF, initial = INF, iterations = 0, scenarios = 0;
  int verify_calls = 0, adversary_calls = 0, schedules = 0;
  long master_nodes = 0;
  int master_variables = 0, master_constraints = 0;
  double seconds = 0, master_seconds = 0, build_seconds = 0;
  double extraction_seconds = 0, verify_seconds = 0, time_best = 0;
  std::vector<int> assignment;
};
struct OracleStats {
  long calls = 0, mip_calls = 0;
  double seconds = 0, mip_seconds = 0;
};
void set_oracle(GRBEnv*);
void reset_oracle_stats();
OracleStats oracle_stats();
DetResult deterministic_mip(const Local&, int, GRBEnv&, Clock&);
int time_horizon(const Instance&);
RobustOptions baseline_verifier_options();
Instance generate_paper(int n, int m, unsigned seed, double gamma,
                        double variation, double capacity_scale, double due_scale,
                        int processing_max = 20, int energy_max = 20,
                        double recharge_factor = 2, int recharge_time = -1);
std::vector<std::vector<int>> identical_machine_groups(const Instance&);
std::vector<int> master_block_limits(const Instance&, int exclusive_target_upper);

// A finite-scenario relaxation. Fixed assignments are used only by validation.
class ScenarioMaster {
 public:
  struct Answer {
    Status status = Status::Unknown;
    std::string reason;
    std::vector<int> assignment;
    long nodes = 0;
  };
  virtual ~ScenarioMaster() = default;
  virtual bool add_scenario(const std::vector<double>&, Clock&) = 0;
  virtual Answer solve(int target, Clock&) = 0;
  virtual int scenarios() const = 0;
  virtual int variables() = 0;
  virtual int constraints() = 0;
  virtual double extraction_time() const { return 0; }
  double build_seconds = 0;
};
class Master : public ScenarioMaster {
 public:
  Master(const Instance&, GRBEnv&, const Options&, Clock&,
         const std::vector<int>* fixed = nullptr);
  bool add_scenario(const std::vector<double>&, Clock&) override;
  Answer solve(int target, Clock&) override;
  int scenarios() const override { return static_cast<int>(points.size()); }
  int variables() override;
  int constraints() override;
 private:
  const Instance& a;
  Options opt;
  GRBModel model;
  GRBVar tardiness;
  int horizon;
  std::vector<int> block_limits;
  std::vector<std::vector<GRBVar>> x;
  std::vector<std::vector<double>> points;
  std::vector<int> order, rank;
  void add_blocks(const std::vector<double>&, Clock&);
  void add_positions(const std::vector<double>&, Clock&);
};
std::unique_ptr<ScenarioMaster> make_master(const Instance&, GRBEnv&, const Options&,
  Clock&, const std::vector<int>* fixed = nullptr);
Result solve(const Instance&, GRBEnv&, const Options&);
void write_solution(const fs::path&, const Result&);
void validate(GRBEnv&, bool quick);
std::string cp_version();
}
