// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"
#include <ilcp/cp.h>

namespace baseline {
namespace {
struct CPEnvironment {
  IloEnv env;
  ~CPEnvironment() { env.end(); }
};
struct CPEngine {
  IloCP cp;
  explicit CPEngine(IloEnv env) : cp(env) {}
  ~CPEngine() { cp.end(); }
};
void check_time(Clock& budget) {
  if (budget.expired()) throw std::runtime_error("construction_time_limit");
}

// Integer bin assignments, native packing, and optional variable-length
// intervals. Floating-point energy coefficients are retained, not rounded.
class CPMaster final : public ScenarioMaster {
  CPEnvironment owner;
  IloEnv env;
  IloModel model;
  IloIntVarArray machine;
  IloIntVarArray decisions;
  IloIntVar tardiness;
  std::unique_ptr<CPEngine> engine;
  const Instance& a;
  Options opt;
  int horizon, nvars=0, nconstraints=0;
  double extraction_seconds=0;
  std::vector<int> block_limits;
  std::vector<int> order;
  std::vector<std::vector<double>> points;
 public:
  CPMaster(const Instance& data, const Options& options,
           Clock& budget, const std::vector<int>* fixed)
    : env(owner.env), model(env), machine(env, data.n, 0, data.m-1), decisions(env),
      tardiness(env, 0, time_horizon(data)), a(data), opt(options), horizon(time_horizon(data)),
      block_limits(master_block_limits(data, options.master_target_upper)) {
    Clock build; check_time(budget);
    if(opt.threads!=1)throw std::runtime_error("CP master requires one worker");
    model.add(machine); model.add(tardiness);
    decisions.add(machine);
    if (fixed) {
      for (int j=0;j<a.n;++j) model.add(machine[j]==(*fixed)[j]);
    }
    if (!fixed) for (const auto& group : identical_machine_groups(a))
      for (std::size_t q=1;q<group.size();++q)
        for (int j=0;j<a.n;++j) {
          IloIntExpr predecessor(env);
          for (int i=0;i<j;++i) predecessor+=(machine[i]==group[q-1]);
          model.add((machine[j]==group[q])<=predecessor);
          predecessor.end();
        }
    order.resize(a.n); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a.d[i]<a.d[j];});
    build_seconds+=build.elapsed();
  }
  bool add_scenario(const std::vector<double>& point, Clock& budget) override {
    try {
      if (point.size()!=static_cast<std::size_t>(a.n)) throw std::runtime_error("scenario dimension");
      double sum=0;
      for(double v:point) {
        if(!std::isfinite(v)||v<0||v>1) throw std::runtime_error("scenario coordinate");
        sum+=v;
      }
      if(sum>a.gamma+1e-10) throw std::runtime_error("scenario budget");
      if(std::find(points.begin(),points.end(),point)!=points.end()) return false;
      // CP Optimizer is reused while only the target changes. A new scenario
      // changes the model and requires one fresh extraction of the enlarged model.
      engine.reset();
      Clock build;
      for(int k=0;k<a.m;++k) {
        check_time(budget);
        int total=0;
        IloIntArray sizes(env);
        for(int j=0;j<a.n;++j) {sizes.add(a.p[j][k]);total+=a.p[j][k];}
        const int slots=block_limits[k];
        // Last bucket contains jobs assigned to another AGV; it has no interval.
        IloIntVarArray where(env,a.n,0,slots), load(env,slots+1,0,total);
        decisions.add(where);
        model.add(IloPack(env,load,where,sizes));
        for(int j=0;j<a.n;++j) {
          model.add((where[j]<slots)==(machine[j]==k));
        }
        IloIntervalVarArray blocks(env);
        for(int b=0;b<slots;++b) {
          IloIntervalVar block(env); block.setOptional();
          block.setStartMin(0);block.setEndMax(horizon);block.setSizeMin(1);block.setSizeMax(total);
          blocks.add(block);
          model.add(IloLengthOf(block,0)==load[b]);
        }
        model.add(IloNoOverlap(env,blocks));
        model.add(IloStartOf(blocks[0],0)==0);
        for(int b=0;b<slots;++b) {
          check_time(budget);
          if(b) {
            model.add(IloPresenceOf(env,blocks[b])<=IloPresenceOf(env,blocks[b-1]));
            model.add(IloEndAtStart(env,blocks[b-1],blocks[b],a.alpha));
          }
          IloNumExpr energy(env);
          for(int j=0;j<a.n;++j)
            energy+=(a.nominal[j][k]+a.deviation[j][k]*point[j])/a.capacity[k]*(where[j]==b);
          model.add(energy<=IloPresenceOf(env,blocks[b]));
          energy.end();
          IloIntExpr inside(env);
          for(int j:order) {
            inside+=a.p[j][k]*(where[j]==b);
            model.add(IloIfThen(env,where[j]==b,IloStartOf(blocks[b],0)+inside<=a.d[j]+tardiness));
          }
          inside.end();
        }
      }
      points.push_back(point);build_seconds+=build.elapsed();return true;
    } catch(const IloException& e) {throw std::runtime_error(std::string("CP model: ")+e.getMessage());}
  }
  Answer solve(int target, Clock& budget) override {
    Answer answer;
    if(budget.expired()) {answer.reason="time_limit";return answer;}
    try {
      // Feasibility is monotone in T. Fix the trial value, as in the MILP;
      // otherwise CP may unnecessarily branch on smaller tardiness values.
      tardiness.setBounds(target,target);
      if(!engine) {
        Clock extraction;
        engine=std::make_unique<CPEngine>(env);IloCP& extracted=engine->cp;
        extracted.setParameter(IloCP::Workers,opt.threads);
        if(extracted.getParameter(IloCP::Workers)!=1)throw std::runtime_error("CP worker-count mismatch");
        extracted.setParameter(IloCP::RandomSeed,static_cast<IloInt>(opt.seed));
        extracted.setParameter(IloCP::LogVerbosity,IloCP::Quiet);
        extracted.setParameter(IloCP::TimeMode,IloCP::ElapsedTime);
        extracted.extract(model);
        // Resolve assignment/packing decisions before their derived interval
        // times. Variable and value selection remain CP Optimizer defaults.
        extracted.setSearchPhases(IloSearchPhase(env,decisions));
        nvars=extracted.getInfo(IloCP::NumberOfIntegerVariables)+extracted.getInfo(IloCP::NumberOfIntervalVariables);
        nconstraints=extracted.getInfo(IloCP::NumberOfConstraints);
        extraction_seconds+=extraction.elapsed();
      }
      IloCP& cp=engine->cp;
      if(budget.expired()) {answer.reason="time_limit";return answer;}
      cp.setParameter(IloCP::TimeLimit,budget.left());
      bool found=cp.solve();answer.nodes=cp.getInfo(IloCP::NumberOfBranches);
      if(found) {
        answer.status=Status::Feasible;
        for(int j=0;j<a.n;++j)answer.assignment.push_back(cp.getValue(machine[j]));
      } else if(cp.getStatus()==IloAlgorithm::Infeasible) answer.status=Status::Infeasible;
      else answer.reason=(budget.expired()||cp.getInfo(IloCP::SearchStopCause)==IloCP::SearchStoppedByLimit)?"time_limit":"cp_unresolved";
      return answer;
    } catch(const IloException& e) {throw std::runtime_error(std::string("CP solve: ")+e.getMessage());}
  }
  int scenarios() const override {return static_cast<int>(points.size());}
  int variables() override {return nvars;}
  int constraints() override {return nconstraints;}
  double extraction_time() const override {return extraction_seconds;}
};
}
std::unique_ptr<ScenarioMaster> make_master(const Instance& a,GRBEnv& env,const Options& opt,
  Clock& budget,const std::vector<int>* fixed) {
  if(opt.formulation!="cp")return std::make_unique<Master>(a,env,opt,budget,fixed);
  try {return std::make_unique<CPMaster>(a,opt,budget,fixed);}
  catch(const IloException& e) {throw std::runtime_error(std::string("CP initialization: ")+e.getMessage());}
}
std::string cp_version() {
  CPEnvironment owner;CPEngine engine(owner.env);return engine.cp.getVersion();
}
}
