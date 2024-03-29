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

#include <armadillo>
#include <chrono>
#include <vector>
#include "behavior.h"
#include "biodynamo.h"
#include "custom_ops.h"
#include "sim-param.h"
#include <TH1F.h>

namespace bdm {

typedef std::vector<real_t> stdvec;
typedef std::vector<std::vector<real_t>> stdvecvec;

inline stdvecvec mat_to_std_vec_vec(arma::mat& M) {
  stdvecvec SVV(M.n_cols);
  for (size_t i = 0; i < M.n_cols; ++i) {
    SVV[i] = arma::conv_to<stdvec>::from(M.col(i));
  };
  return SVV;
}

inline stdvec mat_to_std_vec(arma::mat& M) {
  stdvec SV(1);
  if (M.n_cols == 1) {
    SV = arma::conv_to<stdvec>::from(M.col(0));
  } else {
    Log::Fatal("Tried to convert matrix with n_col > 1 to vector");
  }
  return SV;
}

inline arma::mat std_vec_vec_to_mat(stdvecvec& SVV) {
  // mat(n_rows, n_cols)
  arma::mat M(SVV[0].size(), SVV.size(), arma::fill::zeros);
  for (size_t i = 0; i < SVV.size(); ++i) {
    M.col(i) = arma::conv_to<arma::vec>::from(SVV[i]);
  };
  return M;
}

inline arma::mat std_vec_to_vec(stdvec& SV) {
  arma::vec V(SV.size(), arma::fill::zeros);
  V = arma::conv_to<arma::vec>::from(SV);

  return V;
}

inline stdvec vec_to_std_vec(arma::vec& V) {
  stdvec SV(V.size(), 0.0);
  SV = arma::conv_to<stdvec>::from(V);

  return SV;
}

inline std::vector<real_t> SetupRunSimulation(stdvec parameters) {
  auto set_param = [&](Param* param) {
    // BDM default pars
    param->export_visualization = false;
    param->visualize_agents["Cell"] = {};
    param->remove_output_dir_contents = true;

    // Custom pars
    auto* sparam = param->Get<SimParam>();
    sparam->division_rate = parameters[0];
    sparam->apoptosis_rate = parameters[1];
  };

  Simulation sim("test", set_param);

  // Substances
  // const real_t sub_time_res = (double)sparam->steps_per_day;
  // const int subs_resolution = sparam->sub_res;
  // struct FlatConc {
  //   double initial_conc_;
  //   FlatConc(double initial_conc) { initial_conc_ = initial_conc; }
  //   double operator()(double x, double y, double z) { return 0.0; }
  // };
  // ModelInitializer::DefineSubstance(
  //         kTGFb, "TGFb",
  //         sparam->TGFb_ac_D / sub_time_res, sparam->TGFb_ac_d / sub_time_res,
  //         subs_resolution);
  // ModelInitializer::InitializeSubstance(sparam->TGFb_ac_flat_conc,
  //                                       FlatConc());

  std::vector<real_t> results_x;

  // Create initial model
  auto* rm = sim.GetResourceManager();
  auto* param = sim.GetParam();
  auto* sparam = param->Get<SimParam>();
  // Log::Warning("Division rate ", sparam->division_rate);
  Cell* cell = new Cell(sparam->diam);
  cell->AddBehavior(new Divide());
  cell->AddBehavior(new Apoptosis());
  // cell->AddBehavior(new Secretion());
  rm->AddAgent(cell);

  auto* count_cells = NewOperation("count_cells");
  count_cells->frequency_ = sparam->count_cell_freq;
  sim.GetScheduler()->ScheduleOp(count_cells);
  // auto* get_conc = NewOperation("measure_conc");
  // get_conc->frequency_ = sparam->count_cell_freq;
  // sim.GetScheduler()->ScheduleOp(get_conc);
  sim.GetScheduler()->Simulate(sparam->simulation_time);
  results_x = count_cells->GetImplementation<CountCells>()->GetMeasurements();

  return results_x;
}

inline real_t GetDistance(std::vector<real_t> simulation_results_x) {
  std::array<int, 10> experimental_data_y = {1,  5,  10, 15, 20,
                                             25, 30, 35, 40, 45};
  // std::array<real_t, 10> experimental_data_y = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  // for (size_t i = 0; i < experimental_data_y.size(); i++)
  // {
  //   experimental_data_y[i]=experimental_data_y[i]*2.51*pow(10,-12);
  // }
  
  real_t distance = 0.0;
  size_t int_counter = 0;

  for (auto res : simulation_results_x) {
    distance += (pow(log2(experimental_data_y[int_counter]) - log2(res), 2));
    int_counter++;
  }
  return distance;
}

inline void SortParamAndDistance(stdvecvec& params, stdvec& distances) {
  std::vector<std::pair<stdvec, real_t>> params_distances;

  // Merge
  for (size_t i = 0; i < params.size(); i++) {
    params_distances.push_back(std::make_pair(params[i], distances[i]));
  }

  // Sort parameters based on distance
  std::sort(params_distances.begin(), params_distances.end(),
            [](auto& left, auto& right) { return left.second < right.second; });

  // Split
  for (size_t i = 0; i < params_distances.size(); i++) {
    params[i] = params_distances[i].first;
    distances[i] = params_distances[i].second;
  }
}

inline void SampleRandomParams(stdvecvec old_parameters,
                               stdvecvec& new_parameters, int discarded,
                               Random* rnd, stdvec old_distances,
                               stdvec& new_distances) {
  for (size_t i = 0; i < discarded; i++) {
    int rnd_index = (int)std::floor(rnd->Uniform(0, old_parameters.size()));
    new_parameters[i] = old_parameters[rnd_index];
    new_distances[i] = old_distances[rnd_index];
  }
}

// Define custom priors and get params
inline stdvec DrawParamsFromPrior(Random* rnd, int nb_params) {
  stdvec parameters(nb_params, 0.0);

  //parameters[0] = rnd->Uniform(0., 0.3);  // Division rate
  parameters[0] = rnd->Uniform(0., 1.0);  // Division rate
  parameters[1] = rnd->Uniform(0., 1.0);  // Apoptosis rate

  return parameters;
}

// Define custom prior and get param probability
inline real_t GetPriorProbability(stdvec& parameters) {
  real_t joint_probability = 1.0;
  stdvec prob_vector(parameters.size(), 0.0);

  // Prior probability
  if (parameters[0] >= 0.0 && parameters[0] <= 1.0) {
    prob_vector[0] = 1.0 / (1.0 - 0.0);
  } else {
    prob_vector[0] = 0.0;
  }
  if (parameters[1] >= 0.0 && parameters[1] <= 1.0) {
    prob_vector[1] = 1.0 / (1.0 - 0.0);
  } else {
    prob_vector[1] = 0.0;
  }

  // Compute joint probability
  for (auto prob : prob_vector) {
    joint_probability *= prob;
  }

  return joint_probability;
}

inline real_t MCMCMoveStep(arma::vec new_proposed_params_vec,
                           stdvecvec& proposed_parameters, real_t threhsold,
                           stdvec& updated_distances, int particle, Random* rnd,
                           int nb_params) {

  stdvec new_proposed_params = vec_to_std_vec(new_proposed_params_vec);

  std::vector<real_t> simulation_results_x =
      SetupRunSimulation(new_proposed_params);

  real_t new_distance = GetDistance(simulation_results_x);

  real_t distance_below_threshold = 0.0;

  if (new_distance < threhsold) {
    distance_below_threshold = 1.0;
  }

  real_t prior_ratio = GetPriorProbability(new_proposed_params) /
                       GetPriorProbability(proposed_parameters[particle]);

  real_t mcmc_acceptance =
      std::min(1.0, distance_below_threshold * prior_ratio);
  if (rnd->Uniform(0., 1.) < mcmc_acceptance) {
    proposed_parameters[particle] = new_proposed_params;
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
  int initial_number_of_particles = 300;  // 1000
  const int number_of_parameters = 2;
  const real_t fraction_rejected_thresholds = 0.5;
  const real_t minimum_mcmc_acceptance_rate = 0.01;
  const real_t prob_not_moving_particle = 0.01;
  int step = 1;                       // t
  real_t trial_MCMC_iterations = 30;  // 30.;
  real_t MCMC_iterations = 0;
  const real_t target_tolerance = 0.01;

  std::vector<std::vector<int>> simulation_results_x;
  stdvecvec parameters;
  stdvec distances;

  for (size_t i = 0; i < initial_number_of_particles; i++) {
    parameters.push_back(DrawParamsFromPrior(random, number_of_parameters));
    simulation_results_x.push_back(SetupRunSimulation(parameters.back()));
    distances.push_back(GetDistance(simulation_results_x.back()));
  }

  SortParamAndDistance(parameters, distances);

  int discarded_particles = (int)std::floor(initial_number_of_particles *
                                            fraction_rejected_thresholds);
  step++;
  real_t threhsold_t_1 = distances[initial_number_of_particles -
                                   discarded_particles - 1];  // Next threshold
  real_t threhsold_t =
      distances[initial_number_of_particles - 1];  // Current threshold

  std::cout << "initial pars \n";
  for (auto val : parameters) {
    for (auto par : val) {
      std::cout << par << "\t";
    }
    std::cout << "\n";
  }

  // Main algorithm loop
  while (threhsold_t > target_tolerance) {
    std::cout << "Last threshold " << threhsold_t << ", accepted threshold "
              << threhsold_t_1 << std::endl;
    stdvecvec accepted_parameters = {parameters.begin(),
                                     parameters.end() - discarded_particles};
    stdvec accepted_distances = {distances.begin(),
                                 distances.end() - discarded_particles};
    arma::mat accepted_parameters_mat = std_vec_vec_to_mat(accepted_parameters);
    inplace_trans(
        accepted_parameters_mat);  // We must traponse the matrix since the cov
                                   // function needs measurements in the rows
                                   // and vars in the columns
    arma::mat cov_matrix = cov(accepted_parameters_mat);

    stdvec updated_distances(discarded_particles,
                             0.0);  // Init vector
    stdvecvec proposed_parameters(
        discarded_particles, stdvec(number_of_parameters, 0.0));  // Init vector
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
      arma::vec proposed_parameters_vec =
          std_vec_to_vec(proposed_parameters[j]);
      arma::mat new_proposed_params_mat =
          arma::mvnrnd(proposed_parameters_vec, cov_matrix,
      (int)trial_MCMC_iterations);  // Generate proposed parameters
      for (size_t k = 0; k < (int)trial_MCMC_iterations; k++)  // Loop MCMC
      {
        mcmc_trial_acceptance +=
            MCMCMoveStep(new_proposed_params_mat.col(k), proposed_parameters, threhsold_t_1,
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
      arma::vec proposed_parameters_vec =
          std_vec_to_vec(proposed_parameters[j]);
      arma::mat new_proposed_params_mat =
          arma::mvnrnd(proposed_parameters_vec, cov_matrix,
      (int)additional_MCMC_iterations);  // Generate proposed parameters
      for (size_t k = 0; k < (int)additional_MCMC_iterations; k++)  // Loop MCMC
      {
        mcmc_acceptance +=
            MCMCMoveStep(new_proposed_params_mat.col(k), proposed_parameters, threhsold_t_1,
                         updated_distances, j, random, number_of_parameters);
      }
      if (additional_MCMC_iterations > 0) {
        std::cout << "Finished additional MCMC for particle " << j << std::endl;
      }
    }

    // Update MCMC acceptance rate
    mcmc_acceptance /= (MCMC_iterations * discarded_particles);

    // Merge vectors
    parameters.clear();
    distances.clear();
    parameters.insert(parameters.end(), accepted_parameters.begin(),
                      accepted_parameters.end());
    parameters.insert(parameters.end(), proposed_parameters.begin(),
                      proposed_parameters.end());
    distances.insert(distances.end(), accepted_distances.begin(),
                     accepted_distances.end());
    distances.insert(distances.end(), updated_distances.begin(),
                     updated_distances.end());

    // Sort
    SortParamAndDistance(parameters, distances);

    trial_MCMC_iterations = std::ceil(MCMC_iterations / 2.);

    threhsold_t_1 = distances[initial_number_of_particles -
                              discarded_particles - 1];  // Next threshold
    threhsold_t =
        distances[initial_number_of_particles - 1];  // Current threshold

    std::cout << "FINISHED STEP " << step << std::endl;
    step++;
    // Exit if minimum acceptance reached
    if (mcmc_acceptance < minimum_mcmc_acceptance_rate) {
      Log::Warning("MCMC acceptance < minimum MCMC acceptance rate!");
      break;
    }
  }

  std::cout << "Final params and distances" << std::endl;
  TCanvas *c1 = new TCanvas("c1","Candle Decay",800,600);
  c1->Divide(number_of_parameters,1);
  TH1F *h[number_of_parameters];
  h[0] = new TH1F("Division", "Division rate", 20, 0, 1);
  h[1] = new TH1F("Apoptosis", "Apoptosis rate", 20, 0, 1);
  for (size_t particle = 0; particle < parameters.size(); particle++) {
    std::cout << "Parameters \t";
    for (size_t param = 0; param < parameters[particle].size(); param++) {
      std::cout << parameters[particle][param] << "\t";
      h[param]->Fill(parameters[particle][param]);
    }
    std::cout << "Distances \t" << distances[particle] << std::endl;
  }
  
  for (size_t subc = 0; subc < number_of_parameters; subc++)
  {
    c1->cd(subc+1);
    h[subc]->Draw();
  }
  
  c1->SaveAs("/Users/marcodurante/smc_abc/histo.pdf");
  c1->SaveAs("/Users/marcodurante/smc_abc/histo.root");

  std::cout << "Simulation completed successfully!\n";
  return 0;
}

}  // namespace bdm

#endif  // SMC_ABC_H_
