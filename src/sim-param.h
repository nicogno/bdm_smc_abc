#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include <TStyle.h>
#include "core/param/param_group.h"

namespace bdm {

/// This class defines parameters that are specific to this simulation.
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);
  real_t diam = 3.14;
};

}  // namespace bdm

#endif  // SIM_PARAM_H_
