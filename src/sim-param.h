#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include <TStyle.h>
#include "core/param/param_group.h"
#include "biodynamo.h"

namespace bdm {

/// This class defines parameters that are specific to this simulation.
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);
  real_t diam = 10.;
  int measurements = 10;
  int simulation_time = 1000;
  int count_cell_freq = 100;
  real_t division_rate = 0.0;
  real_t apoptosis_rate = 0.0;
  real_t TGFb_ac_D = pow(10, 8) * 0.0432;  // um^2/day
  real_t TGFb_ac_d = 333;  // day^(-1) //
  //real_t TGFb_ac_flat_conc = 2.51 * pow(10, -12);
  real_t TGFb_ac_flat_conc = 0.0;
  int sub_res = 8;
  uint64_t hours_per_day = 24;
  uint64_t seconds_per_step = 200;
  uint64_t steps_per_day = hours_per_day * 60 * 60 / seconds_per_step; // Divided by the simulation_time_step (200 s)
};

}  // namespace bdm

#endif  // SIM_PARAM_H_
