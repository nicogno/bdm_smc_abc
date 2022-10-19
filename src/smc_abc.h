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

inline std::vector<int> SetupRunSimulation(std::vector<real_t> parameters) {
  auto set_param = [&](Param* param) {
    // BDM default pars
    param->export_visualization = false;
    param->visualize_agents["Cell"] = {};
    param->remove_output_dir_contents = true;

    // Custom pars
    auto* sparam = param->Get<SimParam>();
    sparam->division_rate = parameters[0];
  };

  Simulation sim("test", set_param);

  std::vector<int> results_x;

  // Create initial model
  auto* rm = sim->GetResourceManager();
  auto* param = sim->GetParam();
  auto* sparam = param->Get<SimParam>();
  // Log::Warning("Division rate ", sparam->division_rate);
  Cell* cell = new Cell(sparam->diam);
  cell->AddBehavior(new Divide());
  cell->AddBehavior(new Apoptosis());
  rm->AddAgent(cell);

  auto* param = sim->GetParam();
  auto* sparam = param->Get<SimParam>();
  auto* count_cells = NewOperation("count_cells");
  count_cells->frequency_ = sparam->count_cell_freq;
  sim->GetScheduler()->ScheduleOp(count_cells);
  sim->GetScheduler()->Simulate(sparam->simulation_time);
  results_x.push_back(
      count_cells->GetImplementation<CountCells>()->GetMeasurements());

  return results_x;
}

inline real_t GetDistance(std::vector<int> simulation_results_x) {
  std::array<int, 10> experimental_data_y = {1,  5,  10, 15, 20,
                                             25, 30, 35, 40, 45};
  real_t distance = 0.0;

  for (auto res : simulation_results_x) {
    distance += (pow(log2(experimental_data_y[int_counter]) - log2(res), 2));
    int_counter++;
  }
  return distance;
}

inline void SortParamAndDistance(std::vector<std::vector<real_t>>& params,
                                 std::vector<real_t>& distances) {
  std::vector<std::pair<std::vector<real_t>, real_t>> params_distances;

  // Merge
  for (size_t i = 0; i < params.size(); i++) {
    params_distances.push_back(std::make_pair(params[i], distances[i]));
  }

  // Sort parameters based on distance
  std::sort(param_distance.begin(), param_distance.end(),
            [](auto& left, auto& right) { return left.second < right.second; });

  // Split
  for (size_t i = 0; i < params_distances.size(); i++) {
    params[i] = params_distances[i].first;
    distances[i] = params_distances[i].second;
  }
}

inline real_t GetMean(std::vector<real_t> values) {
  real_t mean = 0.;
  for (auto val : values) {
    mean += val;
  }

  return mean / values.size();
}

inline real_t GetVariance(std::vector<real_t> values_1, std::vector<real_t> values_2, real_t mean_1, real_t mean_2) {
  real_t variance = 0.;
  for (size_t i = 0; i < count; i++)
  {
    variance += ((values_1[i]-mean_1)*(values_2[i]-mean_2));
  }

  return variance / values_1.size();
}

inline std::vector<std::vector<real_t>> GetCovarianceMatrix(std::vector<std::vector<real_t>> parameters, int number_of_pars) {
  
  /////////////////////////////////////////
  //////      ( s(x. x) s(x, y) )
  ////// C =  (                 )
  //////      ( s(y, x) s(y, y) )
  /////////////////////////////////////////
  /// In parameters, each outer element is a parameter and for each parameter we have N measurements (inner elements)
  std::vector<std::vector<real_t>> covariance_matrix(number_of_pars, std::vector<real_t>(number_of_pars, 0.0));
  std::vector<real_t> means;

  for (auto par : parameters)
  {
    means.push_back(GetMean(par));
  }

  for (size_t i = 0; i < parameters.size(); i++)
  {
    for (size_t j = 0; j < parameters.size(); j++) 
    {
      covariance_matrix[i][j] = GetVariance(parameters[i], parameters[j], means[i], means[j]);
    }
  }

  return covariance_matrix;
}

inline void SampleRandomParams(std::vector<std::vector<real_t>> old_parameters,
                              std::vector<std::vector<real_t>>& new_parameters,
                               int discarded, Random* rnd,
                               std::vector<real_t> old_distances,
                               std::vector<real_t>& new_distances) {
  for (size_t i = 0; i < discarded; i++) {
    int rnd_index = (int)std::floor(rnd->Uniform(0, old_parameters.size()));
    new_parameters[i] = old_parameters[rnd_index];
    new_distances[i] = old_distances[rnd_index];
  }
}

inline real_t MCMCMoveStep(std::vector<std::vector<real_t>> covariance,
                          std::vector<std::vector<real_t>>& proposed_parameters,
                           real_t threhsold,
                           std::vector<real_t>& updated_distances, int particle, Random* rnd, int nb_params) {
  std::vector<real_t> new_proposed_params(nb_params, 0.);
  for (size_t i = 0; i < nb_params; i++) {
    real_t temp_val = 0.0;
    do {
      temp_val = rnd->Gaus(proposed_parameters[particle], covariance);
    } while (temp_val < 0. || temp_val > 1.);
    new_proposed_params[0] = temp_val;
    temp_val = 0.0;
    do {
      temp_val = rnd->Gaus(proposed_parameters[particle], covariance);
    } while (temp_val < 0. || temp_val > 1.);
    new_proposed_params[1] = temp_val;
  }

  std::vector<int> simulation_results_x =
      SetupRunSimulation(new_proposed_params);
  real_t> new_distance =
      GetDistance(simulation_results_x);

  real_t distance_below_threshold = 0;

  if (new_distance < threhsold) {
    distance_below_threshold = 1;
  }

  real_t mcmc_acceptance = std::min(1.0, distance_below_threshold);
  if (rnd->Uniform(0., 1.) < mcmc_acceptance) {
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
  int initial_number_of_particles = 1000;  // 10
  const int number_of_parameters = 2;
  const real_t fraction_rejected_thresholds = 0.5;
  const real_t minimum_mcmc_acceptance_rate = 0.01;
  const real_t prob_not_moving_particle = 0.01;
  int step = 1;                       // t
  real_t trial_MCMC_iterations = 30;  // 30.;
  real_t MCMC_iterations = 0;
  const real_t target_tolerance = 0.01;

  std::vector<real_t> division_rates =
      GetParamsFromPrior(initial_number_of_particles, random);

  std::vector<real_t> apoptosis_rates =
      GetParamsFromPrior(initial_number_of_particles, random);

  std::vector<std::vector<int>> simulation_results_x;
  std::vector<std::vector<real_t>> parameters;
  std::vector<real_t> distances;

  for (size_t i = 0; i < initial_number_of_particles; i++) {
    parameters.push_back({division_rates[i], apoptosis_rates[i]});
    simulation_results_x.push_back(SetupRunSimulation(parameters.back()));
    distances.push_back(GetDistance(simulation_results_x.back()));
  }

  SortParamAndDistance(parameters, distances);

  int discarded_particles = (int)std::floor(initial_number_of_particles *
                                            fraction_rejected_thresholds);
  step++;
  real_t threhsold_t_1 = distances[initial_number_of_particles -
                                                 discarded_particles - 1];  // Next threshold
  real_t threhsold_t = distances[initial_number_of_particles - 1];  // Current threshold

  // Main algorithm loop
  while (threhsold_t > target_tolerance) {
    std::cout << "Last threshold " << threhsold_t << ", accepted threshold "
              << threhsold_t_1 << std::endl;
    std::vector<std::vector<real_t>> accepted_parameters = {
        parameters.begin(),
    parameters.end() - discarded_particles};
    std::vector<real_t> accepted_distances = {
    distances.begin(), distances.end() - discarded_particles};

    std::vector<std::vector<real_t>> cov_matrix = GetCovarianceMatrix(accepted_parameters, number_of_parameters);
    std::vector<real_t> updated_distances(discarded_particles,
                                          0.0);  // Init vector
    std::vector<std::vector<real_t>> proposed_parameters(discarded_particles,
                                            std::vector(number_of_parameters, 0.0));  // Init vector
    SampleRandomParams(accepted_parameters, proposed_parameters,
                       discarded_particles, random, accepted_distances,
                       updated_distances);  // Fill vectors
    real_t mcmc_trial_acceptance = 0.;

    if ((int)trial_MCMC_iterations <= 0) {
      Log::Warning("Trial MCMC iterations = 0!");
      break;
    }

    // First MCMC
    for (size_t j = 0; j < discarded_particles; j++)  // Loop particles
    {
      for (size_t k = 0; k < (int)trial_MCMC_iterations; k++)  // Loop MCMC
      {
        mcmc_trial_acceptance +=
            MCMCMoveStep(cov_matrix, proposed_parameters, threhsold_t_1,
                         updated_distances, j, random, number_of_parameters);
      }
      std::cout << "Finished MCMC for particle " << j << std::endl;
    }

    real_t mcmc_acceptance =
        mcmc_trial_acceptance;  // Must be done before the upcoming division

    // Update MCMC trial acceptance rate
    mcmc_trial_acceptance /= (trial_MCMC_iterations * discarded_particles);

    // Update MCMC iterations
    MCMC_iterations = (int)std::ceil(std::log(prob_not_moving_particle) /
                                     (1 + std::log(1 - mcmc_trial_acceptance)));

    int additional_MCMC_iterations =
        (int)std::max(0, (int)MCMC_iterations - (int)trial_MCMC_iterations);

    std::cout << "Additional MCMC steps " << additional_MCMC_iterations
              << std::endl;

    // Second MCMC
    for (size_t j = 0; j < discarded_particles; j++)  // Loop particles
    {
      for (size_t k = 0; k < additional_MCMC_iterations; k++)  // Loop MCMC
      {
        mcmc_acceptance +=
            MCMCMoveStep(cov_matrix, proposed_parameters, threhsold_t_1,
                         updated_distances, j, random, number_of_parameters);
      }
    }

    // Update MCMC acceptance rate
    mcmc_acceptance /= (MCMC_iterations * discarded_particles);

    // Merge and sort vectors
    sorted_division_rates.clear();
    sorted_distances.clear();
    sorted_division_rates.insert(sorted_division_rates.end(),
                                 accepted_parameters.begin(),
                                 accepted_parameters.end());
    sorted_division_rates.insert(sorted_division_rates.end(),
                                 proposed_parameters.begin(),
                                 proposed_parameters.end());
    sorted_distances.insert(sorted_distances.end(), accepted_distances.begin(),
                            accepted_distances.end());
    sorted_distances.insert(sorted_distances.end(), updated_distances.begin(),
                            updated_distances.end());

    // Update with sorted values
    for (int i = 0; i < sorted_params_distances.size(); i++) {
      sorted_params_distances[i].first = sorted_division_rates[i];
      sorted_params_distances[i].second = sorted_distances[i];
    }

    // Sort
    sorted_params_distances = SortParamAndDistance(sorted_params_distances);

    // Split
    for (int i = 0; i < sorted_params_distances.size(); i++) {
      sorted_division_rates[i] = sorted_params_distances[i].first;
      sorted_distances[i] = sorted_params_distances[i].second;
    }

    trial_MCMC_iterations = std::ceil(MCMC_iterations / 2.);

    threhsold_t_1 =
        sorted_distances[initial_number_of_particles - discarded_particles -
                         1];  // Next threshold
    threhsold_t =
        sorted_distances[initial_number_of_particles - 1];  // Current threshold

    std::cout << "FINISHED STEP " << step << std::endl;
    step++;
    // Exit if minimum acceptance reached
    if (mcmc_acceptance < minimum_mcmc_acceptance_rate) {
      Log::Warning("MCMC acceptance < minimum MCMC acceptance rate!");
      break;
    }
  }

  std::cout << "Final params and distances" << std::endl;
  for (size_t i = 0; i < sorted_division_rates.size(); i++) {
    std::cout << sorted_division_rates[i] << "\t" << sorted_distances[i]
              << std::endl;
  }

  std::cout << "Simulation completed successfully!\n";
  return 0;
}

}  // namespace bdm

#endif  // SMC_ABC_H_
