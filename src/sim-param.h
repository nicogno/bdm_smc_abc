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
};

}  // namespace bdm

#endif  // SIM_PARAM_H_
