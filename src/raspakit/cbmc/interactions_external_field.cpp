module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <thread>
#include <type_traits>
#include <vector>
#endif

module cbmc_interactions_external_field;

#ifdef USE_STD_IMPORT
import std;
#endif

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeExternalFieldEnergy(
    bool hasExternalField, [[maybe_unused]] const ForceField &forceField,
    [[maybe_unused]] const SimulationBox &simulationBox, [[maybe_unused]] double cutOffVDW,
    [[maybe_unused]] double cutOffCoulomb, [[maybe_unused]] std::span<Atom> atoms) noexcept
{
  RunningEnergy energySum;

  const double overlapCriteria = forceField.energyOverlapCriteria;

  if (hasExternalField)
  {
    if (atoms.empty()) return energySum;

    for (std::span<Atom>::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
      [[maybe_unused]] std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
      [[maybe_unused]] std::size_t compA = static_cast<std::size_t>(it1->componentId);
      [[maybe_unused]] std::size_t typeA = static_cast<std::size_t>(it1->type);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
      [[maybe_unused]] bool isFractional = static_cast<bool>(it1->isFractional);
      [[maybe_unused]] double scalingVDWA = it1->scalingVDW;
      [[maybe_unused]] double scaleCoulombA = it1->scalingCoulomb;
      [[maybe_unused]] double chargeA = it1->charge;
      [[maybe_unused]] double3 posA = it1->position;
      [[maybe_unused]] double3 s = (simulationBox.inverseCell * posA).fract();

      // Fill in the energy based on the atom properties and the fractional position 's'
      Potentials::EnergyFactor energyFactor = Potentials::EnergyFactor(0.0, 0.0);
      if (energyFactor.energy > overlapCriteria) return std::nullopt;

      energySum.externalFieldVDW += energyFactor.energy;
      energySum.dudlambdaVDW += energyFactor.dUdlambda;
    }
  }

  return energySum;
}
