// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef SMC_ABC_H_
#define SMC_ABC_H_

#include <vector>
#include "biodynamo.h"
#include "sim-param.h"

namespace bdm {

// Behavior that divides the agent at each time step
struct Divide : Behavior {
  BDM_BEHAVIOR_HEADER(Divide, Behavior, 1);

  Divide() {}

  void Run(Agent* agent) override { dynamic_cast<Cell*>(agent)->Divide(); }
};

inline int Simulate(int argc, const char** argv) {
  Param::RegisterParamGroup(new SimParam());

  std::vector<real_t> diams = {30, 40, 50};
  real_t testpar = 0.;
  auto set_param = [&](Param* param) {
    param->export_visualization = true;
    param->visualize_agents["Cell"] = {};
    param->remove_output_dir_contents = true;

    auto* sparam = param->Get<SimParam>();
    sparam->diam = testpar;
  };

  // Create num simulations
  std::vector<Simulation*> simulations;

  int num = 3;
  for (size_t i = 0; i < num; i++) {
    testpar = diams[i];
    simulations.push_back(
        new Simulation("test" + std::to_string(num), set_param));
  }
  // Initialize the model for each simulation
  for (auto* sim : simulations) {
    // If we switch between simulations we must call the function Activate();
    // Only one simulation can be active at any given time.
    sim->Activate();

    // Create initial model
    auto* rm = sim->GetResourceManager();

    auto* param = sim->GetParam();
    auto* sparam = param->Get<SimParam>();
    Log::Warning("Cell diam ", sparam->diam);
    Cell* cell = new Cell(sparam->diam);
    cell->AddBehavior(new Divide());
    rm->AddAgent(cell);
  }

  // For each simulation simulate 5 timesteps
  for (auto* sim : simulations) {
    sim->Activate();
    sim->GetScheduler()->Simulate(5);
  }

  // Delete sinulations
  for (auto* sim : simulations) {
    sim->Activate();
    delete sim;
  }

  std::cout << "Simulation completed successfully!\n";
  return 0;
}

}  // namespace bdm

#endif  // SMC_ABC_H_
