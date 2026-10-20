// Copyright (c) 2026 Renjie Yu and Daniel Oron
// SPDX-License-Identifier: MIT
#include "baseline.hpp"
#include <iostream>

int main(int argc,char** argv) {try {
  GRBEnv env(true);env.set(GRB_IntParam_OutputFlag,0);env.set(GRB_IntParam_Threads,1);
  env.set(GRB_IntParam_Seed,1);env.start();agv::set_solver_seed(1);baseline::set_oracle(&env);
  baseline::validate(env,argc>1&&std::string(argv[1])=="--quick");return 0;
} catch(const GRBException& e) {std::cerr<<e.getMessage()<<'\n';return 2;}
  catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
