// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#pragma once
#include "robust.hpp"
#include <map>
#include <memory>
#include <random>

namespace longrun {
using namespace agv;
struct Stats;
struct Options {
  // Observation only; no scores, budgets, or search decisions are modified.
  std::function<void(double,int,int,const Stats&,const std::vector<int>&)> progress;
  std::string mode = "pool";
  // Algorithm 2. The historical experimental modes remain available internally.
  bool manuscript = false;
  unsigned seed = 1;
  double limit = 10;
  double local_limit = .025;
  double root_limit = .5;
  // 0: analytic bound only; 1: legacy timed MILP seed; 2: LP bound;
  // 3: one-node MILP bound; 4: solve the MILP relaxation to optimality.
  // Modes 2--4 do not use solver incumbents.
  int root_mode = 3;
  double pool_limit = .3;
  bool adaptive = true, recombination = true, packing = false;
  bool iterated = true, strong_static = true;
  bool staged = false, targeted = false, inheritance = false;
  bool optimistic = false;
  int lower_weight = 4;
  bool regret_repair = false;
  bool cyclic = false, randomized_regret = false;
  bool polishing = false;
  bool alternating = false;
  double warm_fraction = .25;
  double phase_end = 1e9;
  int max_pool = 2500;
  bool fast_screen = false, managed = false, restart = false;
  bool witness_memory = false;
  bool probe_warm = false;
  bool focus = false;
  bool tight_coverage = false;
  double polish_reserve = 0, pulse_limit = .25, pulse_interval = 3;
  double recombine_interval = 0;
  // Controlled component alternatives. Zero preserves the frozen design.
  int decoder = 0; // 0 portfolio, 1 EDD, 2 insertion, 3 slack filling, 4-6 single original rules
  bool marginal_destruction = false, regret_insertion = false;
  bool ejection_chains = false, chain_only = false;
  bool greedy_acceptance = false, event_warm = false;
  bool relocation_moves = true, exchange_moves = true;
  bool post_perturbation_descent = true;
  bool cooling_acceptance = true;
  bool plateau_acceptance = false;
  bool standard_sa_acceptance = false;
  double temperature_scale = .15;
  bool pair_repair = false;
  bool energy_destruction = false;
  bool stop_zero = false;
  bool monotone_deletion = false;
  bool deadline_repair = false, block_removal = false, single_block = false;
  bool recourse_reference = false, scarcity_repair = false, chain_after = false;
  bool replacement_bounds = false;
  bool related_removal = false;
  bool target_chains = false;
  bool zero_aware_refinement = false, tardiness_score = false;
  bool incremental_lower = false, seed_portfolio = false;
  bool critical_prefix = false, lex_acceptance = false, interval_targets = false;
  bool memo_lower = false;
  bool budget_portfolio = false, budget_best = false;
};
struct Stats {
  long entries = 0, cache_hits = 0, evaluations = 0, screened = 0;
  long verifications = 0, unknown = 0, exact_sets = 0, adaptive_improvements = 0;
  long moves = 0, shakes = 0, pool_calls = 0, pool_improvements = 0;
  long relocation_moves = 0, exchange_moves = 0;
  double verification_seconds = 0, mip_seconds = 0;
  double entry_seconds = 0, pulse_seconds = 0;
  long fast_screened = 0, pulses = 0, pulse_improvements = 0, restarts = 0;
  long det_unknown = 0, coverage_unknown = 0;
  long witness_checks = 0, witness_prunes = 0;
  double witness_seconds = 0;
  double adversary_seconds = 0;
  long numerical_unknown = 0;
  long chain_calls = 0, chain_candidates = 0, chain_improvements = 0;
  long chain_jobs = 0, removed_jobs = 0, repair_candidates = 0;
  double chain_seconds = 0, repair_seconds = 0;
  long pair_calls = 0, pair_candidates = 0, pair_changed_jobs = 0;
  double pair_seconds = 0;
  long slack_improvements = 0, slack_gain = 0;
  long deletion_updates = 0, deletion_gain = 0;
  long block_calls = 0, reference_recourse = 0, reference_static = 0;
  long block_removed = 0, block_changed = 0, scarcity_queries = 0;
  long chain_compound = 0, chain_stalls = 0;
  double block_seconds = 0;
  long replacement_upper = 0, replacement_lower = 0, replacement_gain = 0;
  long chain_primary = 0, chain_root_prunes = 0, chain_closed_prunes = 0, chain_ejection_prunes = 0;
  long negative_targets = 0, zero_targets = 0, zero_skipped = 0, targets_capped = 0;
  double negative_target_seconds = 0, zero_target_seconds = 0;
  long portfolio_entries = 0, portfolio_schedules = 0;
  long prioritized_jobs = 0, critical_jobs = 0, priority_scans = 0;
  double portfolio_seconds = 0, lower_seconds = 0, screen_seconds = 0;
  long lower_hits = 0, lower_entries = 0;
  long budget_updates = 0, budget_gain = 0;
};
struct Result {
  int initial = INF, upper = INF, lower = -INF;
  int assignment_lower = -INF, static_upper = INF;
  double seconds = 0, time_best = 0;
  std::vector<int> assignment, local_lower, local_upper;
  Stats stats;
};
struct Entry {
  Local data;
  int lower = -INF, upper = INF, static_upper = INF;
  int attempts = 0;
  double largest_failed_budget = 0;
  int failed_target = INF;
  double failed_target_budget = 0;
  std::vector<Schedule> schedules;
  bool portfolio_seeded = false;
  bool budget_seeded = false;
};
Instance generate_paper(int n, int m, unsigned seed, double gamma, double variation,
                        double capacity, double due, int processing_max = 20,
                        int nominal_energy_max = 20, double recharge_factor = 2.0,
                        int recharge_time = -1, double tardiness_factor = -1.0,
                        double due_date_range = 0.5);
Instance read_instance(const std::string& path);
RobustResult tight_verify(const Local&,int,GRBEnv&,const RobustOptions&,Clock&,std::vector<Schedule>*);
int static_upper(const Local& a, bool strong = true);
int insertion_upper(const Local& a);
int priority_upper(const Local& a, int rule);
int slack_upper(const Local& a, long* improvements=nullptr, long* gain=nullptr);
std::vector<std::vector<int>> static_reference(const Local& a, bool strong = true);
// Expose the existing production insertion decoder's schedule for certificates.
// This does not change HIG's construction, evaluation, or search decisions.
std::vector<std::vector<int>> insertion_reference(const Local& a);
int replacement_penalty(const Instance& a,int machine,int removed,int inserted);
int workload_lower(const Local& a, bool packing = false);
int incremental_lower(const Local& a);
std::vector<Schedule> nominal_portfolio(const Local& a);
std::vector<Schedule> budget_repair_portfolio(const Local& a, bool best_candidates=false);
Options options_for(const std::string& mode);

class Search {
 public:
  Search(const Instance& data, GRBEnv& env, Options options);
  Result run();
  Result audit(const std::vector<int>& assignment);
  Result probe_chain(std::vector<int> assignment);
  Result probe_shake(std::vector<int> assignment);
  Entry& entry(int vehicle, const std::vector<int>& jobs);
  void refine(Entry& e, double seconds, int target = INF);
 private:
  friend struct SearchValidation;
  const Instance& a;
  GRBEnv& env;
  Options opt;
  Clock budget;
  std::mt19937 rng;
  Stats stats;
  int global_lower = -INF;
  std::vector<std::map<std::vector<int>, std::unique_ptr<Entry>>> cache;
  std::vector<std::map<std::vector<int>,int>> lower_cache;
  std::vector<int> best;
  std::vector<int> best_score;
  double time_best = 0;
  double last_pulse = 0, last_pool = 0, last_restart = 0;
  std::vector<std::vector<int>> elite;
  std::vector<std::vector<double>> witnesses;
  std::vector<std::vector<int>> groups(const std::vector<int>& x) const;
  std::vector<int> ranked(std::vector<int> v) const;
  bool accept_reconstruction(const std::vector<int>& current,
                             const std::vector<int>& candidate);
  std::vector<int> score(const std::vector<int>& x, bool refine_values = false,
                         const std::vector<int>* incumbent = nullptr);
  std::vector<int> guide(const std::vector<int>& x);
  void remember(const std::vector<int>& x);
  std::vector<int> construct(int perturbation);
  bool descend(std::vector<int>& x);
  bool chain_descent(std::vector<int>& x);
  std::vector<int> structured_shake(const std::vector<int>& x, int strength);
  std::vector<int> pair_shake(const std::vector<int>& x);
  std::vector<int> block_shake(const std::vector<int>& x);
  void inherit_removals(const std::vector<int>& before, const std::vector<int>& after);
  std::vector<int> shake(const std::vector<int>& x, int strength);
  std::vector<int> master_initialization(const std::vector<int>& start, int upper);
  std::vector<int> recombine(const std::vector<int>& start);
  Result finish(int initial);
  Result run_manuscript();
  bool search_expired() const;
  void polish_incumbent();
  void pulse(double seconds);
  bool cheap_screen(const std::vector<int>& x, const std::vector<int>& incumbent);
};
}
