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

inline std::vector<real_t> GetParamsFromPrior(int particles, Random* rnd) {
  std::vector<real_t> division_rates;
  for (size_t i = 0; i < particles; i++) {
    division_rates.push_back(rnd->Uniform(0., 0.3));
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
    //Log::Warning("Division rate ", sparam->division_rate);
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

  return param_distance;
}

inline std::vector<std::pair<real_t, real_t>> SortParamAndDistance(std::vector<std::pair<real_t, real_t>> param_distance) {

  // Sort parameters based on distance
  std::sort(param_distance.begin(), param_distance.end(),
            [](auto& left, auto& right) { return left.second < right.second; });

  for (auto values : param_distance) {
    std::cout << values.first << "\t\t" << values.second << std::endl;
  }

  return param_distance;
}

inline real_t GetMean(std::vector<real_t> values) {
  real_t mean = 0.;
  for (auto val : values) {
    mean += val;
  }

  return mean / values.size();
}

inline real_t GetCovarianceMatrix(std::vector<real_t> accepted_parameters) {
  real_t covariance_matrix = 0.;

  real_t mean = GetMean(accepted_parameters);

  for (auto measure : accepted_parameters) {
    covariance_matrix += (pow(measure - mean, 2));
  }

  covariance_matrix /= accepted_parameters.size();

  return covariance_matrix;
}

inline std::vector<real_t> SampleRandomParams(
    std::vector<real_t> old_parameters, int discarded, Random* rnd) {
  std::vector<real_t> new_parameters;

  for (size_t i = 0; i < discarded; i++) {
    int rnd_index = (int)std::floor(rnd->Uniform(0, discarded));
    new_parameters.push_back(old_parameters[rnd_index]);
  }

  return new_parameters;
}

inline real_t MCMCMoveStep(real_t covariance,
                         std::vector<real_t> &proposed_parameters,
                         real_t threhsold,
                         std::vector<real_t> &updated_distances,
                         size_t particle, Random* rnd) {
  real_t new_proposed_param = 0.;
  do
  {
    new_proposed_param = rnd->Gaus(proposed_parameters[particle], covariance);
  } while (new_proposed_param < 0. || new_proposed_param > 1.);
  std::vector<real_t> new_proposed_param_vector{new_proposed_param};
  std::vector<std::vector<int>> simulation_results_x =
      SetupRunSimulations(new_proposed_param_vector);
  std::vector<std::pair<real_t, real_t>> params_distances =
      GetParamAndDistance(simulation_results_x, new_proposed_param_vector);
  real_t new_distance =
      params_distances[0].second;  // We ran a single simulation

  real_t distance_below_threshold = 0;

  if (new_distance < threhsold) {
    distance_below_threshold = 1;
  }

  real_t mcmc_acceptance = std::min(1.0, distance_below_threshold);
  if (rnd->Uniform(0., 1.) < mcmc_acceptance)
  {
    proposed_parameters[particle] = new_proposed_param;
    updated_distances[particle] = new_distance;
  }

  return mcmc_acceptance;
}

inline int Simulate(int argc, const char** argv) {
  Param::RegisterParamGroup(new SimParam());

  Simulation simulation_for_random("random");
  auto* random = simulation_for_random.GetRandom();
  const auto time_seed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::system_clock::now().time_since_epoch());
  random->SetSeed(time_seed.count());

  // Set algorithm parameters
  int initial_number_of_particles = 10;
  const int number_of_parameters = 1;
  const real_t fraction_rejected_thresholds = 0.6; //0.5; TODO check
  const real_t minimum_mcmc_acceptance_rate = 0.01;
  const real_t prob_not_moving_particle = 0.01;
  int step = 1;  // t
  real_t trial_MCMC_iterations = 50; //30.;
  real_t MCMC_iterations = 0;
  const real_t target_tolerance = 0.01;

  std::vector<real_t> division_rates =
      GetParamsFromPrior(initial_number_of_particles, random);

  std::vector<std::vector<int>> simulation_results_x =
      SetupRunSimulations(division_rates);

  std::vector<std::pair<real_t, real_t>> params_distances =
      GetParamAndDistance(simulation_results_x, division_rates);

  std::vector<std::pair<real_t, real_t>> sorted_params_distances =
  SortParamAndDistance(params_distances);

  std::vector<real_t> sorted_division_rates, sorted_distances;

  // Update with sorted values
  for (auto pair : sorted_params_distances)
  {
    sorted_division_rates.push_back(pair.first);
    sorted_distances.push_back(pair.second);
  }

  int discarded_particles = (int)std::floor(initial_number_of_particles *
                                            fraction_rejected_thresholds);
  step++;
  real_t threhsold_t_1 =
      sorted_params_distances[discarded_particles - 1].second;  // Next threshold
  real_t threhsold_t = sorted_params_distances[initial_number_of_particles - 1]
                           .second;  // Current threshold

  // Main algorithm loop
  while (threhsold_t > target_tolerance) {
    std::vector<real_t> accepted_parameters = {
    sorted_division_rates.begin(), sorted_division_rates.end() - discarded_particles};
    std::vector<real_t> accepted_distances = {
    sorted_distances.begin(), sorted_distances.end() - discarded_particles};

    real_t cov_matrix = GetCovarianceMatrix(accepted_parameters);
    std::vector<real_t> proposed_parameters =
        SampleRandomParams(accepted_parameters, discarded_particles, random);
    std::vector<real_t> updated_distances(proposed_parameters.size(), std::numeric_limits<float>::max()); // Init vector
    real_t mcmc_trial_acceptance = 0.;

    if ((int)trial_MCMC_iterations <= 0)
    {
     Log::Warning("Trial MCMC iterations = 0!");
     break;
    }

    // First MCMC
    for (size_t j = 0;
         j < discarded_particles; j++)  // Loop particles
    {
      for (size_t k = 0; k < (int)trial_MCMC_iterations; k++)  // Loop MCMC
      {
        mcmc_trial_acceptance+=MCMCMoveStep(cov_matrix, proposed_parameters, threhsold_t,
                     updated_distances, j, random);
      }
      std::cout << "Finished MCMC for particle " << j << std::endl;
    }

    real_t mcmc_acceptance = mcmc_trial_acceptance; // Must be done before the upcoming division

    // Update MCMC trial acceptance rate
    mcmc_trial_acceptance/= (trial_MCMC_iterations*(initial_number_of_particles-discarded_particles));

    // Update MCMC iterations
    MCMC_iterations = (int)std::floor(std::log(prob_not_moving_particle)/(1+std::log(1-mcmc_trial_acceptance)));

    int additional_MCMC_iterations = (int)std::max(0,(int)MCMC_iterations - (int)trial_MCMC_iterations);

    // Second MCMC
    for (size_t j = 0;
         j < discarded_particles; j++)  // Loop particles
    {
      for (size_t k = 0; k < additional_MCMC_iterations; k++)  // Loop MCMC
      {
        mcmc_acceptance+=MCMCMoveStep(cov_matrix, proposed_parameters, threhsold_t,
                     updated_distances, j, random);
      }
    }

    // Update MCMC acceptance rate
    mcmc_acceptance/= (MCMC_iterations*(initial_number_of_particles-discarded_particles));

    // Merge and sort vectors
    sorted_division_rates.clear();
    sorted_distances.clear();
    sorted_division_rates.insert(sorted_division_rates.end(), accepted_parameters.begin(), accepted_parameters.end());
    sorted_division_rates.insert(sorted_division_rates.end(), proposed_parameters.begin(), proposed_parameters.end());
    sorted_distances.insert(sorted_distances.end(), accepted_distances.begin(), accepted_distances.end());
    sorted_distances.insert(sorted_distances.end(), updated_distances.begin(), updated_distances.end());

    // Update with sorted values
    for (int i = 0; i<sorted_params_distances.size(); i++)
    {
      sorted_params_distances[i].first = sorted_division_rates[i];
      sorted_params_distances[i].second = sorted_distances[i];
    }
    
    // Sort
    sorted_params_distances = SortParamAndDistance(sorted_params_distances);

    // Split
      for (int i = 0; i<sorted_params_distances.size(); i++)
      {
      sorted_division_rates[i] = sorted_params_distances[i].first;
      sorted_distances[i] = sorted_params_distances[i].second;
    }
    
    trial_MCMC_iterations = std::floor(MCMC_iterations / 2.);

    threhsold_t_1 =
        params_distances[discarded_particles - 1].second;  // Next threshold
    threhsold_t = params_distances[initial_number_of_particles - 1]
                      .second;  // Current threshold

    std::cout << "FINISHED STEP " << step << std::endl;
    step++;
    // Exit if minimum acceptance reached
    if (mcmc_acceptance < minimum_mcmc_acceptance_rate)
    {
      Log::Warning("MCMC acceptance < minimum MCMC acceptance rate!");
      break;
    }
  }

  std::cout << "Final params and distances" << std::endl;
  for (size_t i = 0; i < sorted_division_rates.size(); i++)
  {
     std::cout << sorted_division_rates[i] << "\t" << sorted_distances[i] << std::endl;
  }
  
  std::cout << "Simulation completed successfully!\n";
  return 0;
}

}  // namespace bdm

#endif  // SMC_ABC_H_
