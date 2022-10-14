#ifndef BEHAVIOR_H_
#define BEHAVIOR_H_

#include "core/behavior/behavior.h"
#include "sim-param.h"

namespace bdm {

// Behavior that divides the agent at each time step
struct Divide : Behavior {
  BDM_BEHAVIOR_HEADER(Divide, Behavior, 1);

  Divide() {}

  void Run(Agent* agent) override {
    auto* random = Simulation::GetActive()->GetRandom();
    real_t ran = random->Uniform(0, 1) * 1.0;
    auto* param = Simulation::GetActive()->GetParam();
    auto* sparam = param->Get<SimParam>();
    if (ran < sparam->division_rate) {
      dynamic_cast<Cell*>(agent)->Divide();
    }
  }
};

}  // namespace bdm

#endif  // BEHAVIOR_H_
