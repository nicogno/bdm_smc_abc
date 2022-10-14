#ifndef CUSTOM_OPS_H_
#define CUSTOM_OPS_H_

#include "biodynamo.h"
#include "core/environment/uniform_grid_environment.h"

namespace bdm {

struct CountCells : public StandaloneOperationImpl {
  // This macro will generate the boilerplate code. It must be included.
  BDM_OP_HEADER(CountCells);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();

    total_cells_.push_back(rm->GetNumAgents());
  }

  std::vector<int> GetMeasurements() { return total_cells_; }

  std::vector<int> total_cells_;
};

}  // namespace bdm

#endif
