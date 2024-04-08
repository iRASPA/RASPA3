module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#endif

export module equation_of_states;

#ifndef USE_LEGACY_HEADERS
import <vector>;
#endif

import component;
import simulationbox;

export struct EquationOfState
{
  enum class Type : int
  {
    PengRobinson = 0,
    PengRobinsonGasem = 1,
    SoaveRedlichKwong = 2
  };
  enum class MultiComponentMixingRules : int
  {
    VanDerWaals = 0
  };
  enum class FluidState : int
  {
    Unknown = 0,
    SuperCriticalFluid = 1,
    Vapor = 2,
    Liquid = 3,
    VaporLiquid = 4
  };

  EquationOfState::FluidState fluidState{ EquationOfState::FluidState::Unknown };
  EquationOfState::Type equationOfState{ EquationOfState::Type::PengRobinson };
  EquationOfState::MultiComponentMixingRules
    multiComponentMixingRules{ EquationOfState::MultiComponentMixingRules::VanDerWaals };

  EquationOfState() = default;

  EquationOfState(EquationOfState::Type type,
                  EquationOfState::MultiComponentMixingRules multiComponentMixingRules,
                  double temperature,
                  double pressure,
                  const SimulationBox &simulationBox,
                  double HeliumVoidFraction,
                  std::vector<Component> &components); 

  bool operator==(EquationOfState const&) const = default;

  void computeComponentFluidProperties(EquationOfState::Type equationOfState,                           
                                       EquationOfState::MultiComponentMixingRules multiComponentMixingRules,         
                                       double temperature,
                                       double pressure,
                                       const SimulationBox &simulationBox,
                                       double HeliumVoidFraction,
                                       std::vector<Component> &components);
};
