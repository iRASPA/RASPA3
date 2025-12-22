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
import double4;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
import interpolation_energy_grid;

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeExternalFieldEnergy(
    bool hasExternalField, [[maybe_unused]] const ForceField &forceField,
    [[maybe_unused]] const SimulationBox &simulationBox, 
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    [[maybe_unused]] double cutOffVDW,
    [[maybe_unused]] double cutOffCoulomb, [[maybe_unused]] std::span<Atom> atoms,
    std::make_signed_t<std::size_t> skip) noexcept
{
  RunningEnergy energySum{};

  //const double overlapCriteria = forceField.energyOverlapCriteria;

  if (hasExternalField)
  {
    if (atoms.empty()) return energySum;

    int index = 0;
    for (std::span<Atom>::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
      if (index != skip)
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

        if (externalFieldInterpolationGrid.has_value())
        {
          energyFactor.energy = scalingVDWA * externalFieldInterpolationGrid->interpolate(posA);
        }
        else
        {
          switch (forceField.potentialEnergySurfaceType)
          {
            case ForceField::PotentialEnergySurfaceType::None:
              break;
            case ForceField::PotentialEnergySurfaceType::GridFile:
              energyFactor.energy = scalingVDWA * externalFieldInterpolationGrid->interpolate(posA);
              break;
            case ForceField::PotentialEnergySurfaceType::MullerBrown:
            {
              double4 A{-200.0, -100.0, -170.0, -15.0};
              double4 a{-1.0, -1.0, -6.5, -0.7};
              double4 x0{1.0, 0.0, -0.5, -1.0};
              double4 b{0.0, 0.0, 11.0, 0.6};
              double4 y0{0.0, 0.5, 1.5, 1.0};
              double4 c{-10.0, -10.0, -6.5, 0.7};
              for (std::size_t i = 0; i < 4; ++i)
              {
                energyFactor.energy += A[i] * std::exp(a[i] * (posA.x - x0[i]) * (posA.x - x0[i]) +
                                                       b[i] * (posA.x - x0[i]) * (posA.y - y0[i]) +
                                                       c[i] * (posA.y - y0[i]) * (posA.y - y0[i]));
              }
            }
            break;
            case ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction:
              energyFactor.energy = posA.x * posA.y * posA.z;
              break;
            default:
              break;
          }

        }
        energySum.externalFieldVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      ++index;
    }
  }

  return energySum;
}
