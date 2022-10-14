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

#include <chrono>
#include <vector>
#include "behavior.h"
#include "biodynamo.h"
#include "custom_ops.h"
#include "sim-param.h"

namespace bdm {

inline std::vector<real_t> GetParamsFromPrior(int particles) {
  std::vector<real_t> division_rates;
  Simulation dummy_simulation("dummy");
  auto* dummy_random = dummy_simulation.GetRandom();
  const auto time_seed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::system_clock::now().time_since_epoch());
  dummy_random->SetSeed(time_seed.count());
  for (size_t i = 0; i < particles; i++) {
    division_rates.push_back(dummy_random->Uniform(0., 0.3));
  }
  return division_rates;
}

inline std::vector<std::vector<int>> SetupRunSimulations(
    std::vector<real_t> parameters) {
  
  real_t current_division_rate = 0.;

  auto set_param = [&](Param* param) {
    // BDM default pars
    param->export_visualization = false;
    param->visualize_agents["Cell"] = {};
    param->remove_output_dir_contents = true;

    // Custom pars
    auto* sparam = param->Get<SimParam>();
    sparam->division_rate = current_division_rate;
  };

  // Create num simulations
  std::vector<Simulation*> simulations;
  std::vector<std::vector<int>> results_x;
  for (size_t i = 0; i < parameters.size(); i++) {
    current_division_rate =
        parameters[i];  // Set custom pars for this simulation
    simulations.push_back(
        new Simulation("test" + std::to_string(i), set_param));
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
    Log::Warning("Division rate ", sparam->division_rate);
    Cell* cell = new Cell(sparam->diam);
    cell->AddBehavior(new Divide());
    rm->AddAgent(cell);
  }

  for (auto* sim : simulations) {
    sim->Activate();
    auto* param = sim->GetParam();
    auto* sparam = param->Get<SimParam>();
    auto* count_cells = NewOperation("count_cells");
    count_cells->frequency_ = sparam->count_cell_freq;
    sim->GetScheduler()->ScheduleOp(count_cells);
    sim->GetScheduler()->Simulate(sparam->simulation_time);
    results_x.push_back(
        count_cells->GetImplementation<CountCells>()->GetMeasurements());
  }

  // Delete simulations
  for (auto* sim : simulations) {
    sim->Activate();
    delete sim;
  }

  return results_x;
}

inline void MCMCMoveStep() {

}

inline std::vector<std::pair<real_t, real_t>> GetParamAndDistance(
    std::vector<std::vector<int>> simulation_results_x,
    std::vector<real_t> parameters) {
  std::vector<std::pair<real_t, real_t>>
      param_distance;  // Pairs of parameters and distances to be evaluated
  std::array<int, 10> experimental_data_y = {1,  5,  10, 15, 20,
                                             25, 30, 35, 40, 45};
  std::vector<real_t> distances;
  int out_counter = 0;

  for (auto res : simulation_results_x) {
    real_t temp_distance = 0.0;
    int int_counter = 0;
    for (auto values : res) {
      std::cout << experimental_data_y[int_counter] << " " << values << "\t\t";
      temp_distance +=
          (pow(log2(experimental_data_y[int_counter]) - log2(values), 2));
      int_counter++;
    }
    distances.push_back(temp_distance);
    std::cout << distances.back() << std::endl;
    param_distance.push_back(
        std::make_pair(parameters[out_counter], distances.back()));
    out_counter++;
  }

  // Sort parameters based on distance
  std::sort(param_distance.begin(), param_distance.end(),
            [](auto& left, auto& right) { return left.second < right.second; });

  for (auto values : param_distance) {
    std::cout << values.first << "\t\t" << values.second << std::endl;
  }

  return param_distance;
}

inline int Simulate(int argc, const char** argv) {
  Param::RegisterParamGroup(new SimParam());
  
  // Set algorithm parameters
  int initial_number_of_particles = 10;
  const real_t fraction_accepted_thresholds = 0.5;
  const real_t minimum_mcmc_acceptance_rate = 0.01;
  const real_t prob_not_moving_particle = 0.01;
  int step = 1; //t

  std::vector<real_t> division_rates =
      GetParamsFromPrior(initial_number_of_particles);

  std::vector<std::vector<int>> simulation_results_x = SetupRunSimulations(division_rates);

  std::vector<std::pair<real_t, real_t>> params_distances =
      GetParamAndDistance(simulation_results_x, division_rates);

  int new_number_of_particles = (int)std::floor(initial_number_of_particles*fraction_accepted_thresholds);
  step++;
  real_t threhsold_t_1 = params_distances[new_number_of_particles-1].second;  // Next threshold
  real_t threhsold_t = params_distances[initial_number_of_particles-1].second; // Current threshold
  real_t trial_MCMC_iterations = 30.;
  real_t MCMC_iterations = 0;
  const real_t target_tolerance = 0.01;

  while (threhsold_t > target_tolerance)
  {



    trial_MCMC_iterations = std::floor(MCMC_iterations/2.);
    threhsold_t_1 = params_distances[new_number_of_particles-1].second;  // Next threshold
    threhsold_t = params_distances[initial_number_of_particles-1].second; // Current threshold
    step++;
  }
  

  std::cout << "Simulation completed successfully!\n";
  return 0;
}

}  // namespace bdm

#endif  // SMC_ABC_H_
