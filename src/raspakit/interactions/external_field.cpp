module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>
#include <vector>
#endif

module interactions_external_field;

#ifdef USE_STD_IMPORT
import std;
#endif

import int3;
import double3;
import double4;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import gradient_factor;
import running_energy;
import component;
import forcefield;
import tricubic_derivatives_external_field;
import triquintic_derivatives_external_field;
import interpolation_energy_grid;

void Interactions::computeExternalFieldEnergy(bool hasExternalField, [[maybe_unused]] const ForceField &forceField,
        [[maybe_unused]] const SimulationBox &simulationBox,
        [[maybe_unused]] std::span<const Atom> moleculeAtoms,
        [[maybe_unused]] RunningEnergy &energyStatus,
        [[maybe_unused]] const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid) noexcept
{
  if (hasExternalField)
  {
    if (moleculeAtoms.empty()) return;

    for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
    {
      [[maybe_unused]] std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
      [[maybe_unused]] std::size_t compA = static_cast<std::size_t>(it1->componentId);
      [[maybe_unused]] std::size_t typeA = static_cast<std::size_t>(it1->type);
      [[maybe_unused]] bool isFractional = static_cast<bool>(it1->isFractional);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
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
      energyStatus.externalFieldVDW += energyFactor.energy;
      energyStatus.dudlambdaVDW += energyFactor.dUdlambda;
    }
  }
}

void Interactions::computeExternalFieldTailEnergy(bool hasExternalField, [[maybe_unused]] const ForceField &forceField,
                                                  [[maybe_unused]] const SimulationBox &simulationBox,
                                                  [[maybe_unused]] std::span<const Atom> moleculeAtoms,
                                                  [[maybe_unused]] RunningEnergy &energyStatus) noexcept
{
  if (hasExternalField)
  {
  }
}


[[nodiscard]] std::optional<RunningEnergy> Interactions::computeExternalFieldEnergyDifference(
    bool hasExternalField, [[maybe_unused]] const ForceField &forceField,
    [[maybe_unused]] const SimulationBox &simulationBox, 
    [[maybe_unused]] const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    [[maybe_unused]] std::span<const Atom> newatoms, [[maybe_unused]] std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum;

  const double overlapCriteria = forceField.energyOverlapCriteria;

  if (hasExternalField)
  {
    for (std::span<const Atom>::iterator it1 = newatoms.begin(); it1 != newatoms.end(); ++it1)
    {
      [[maybe_unused]] double3 posA = it1->position;
      [[maybe_unused]] std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
      [[maybe_unused]] std::size_t compA = static_cast<std::size_t>(it1->componentId);
      [[maybe_unused]] std::size_t typeA = static_cast<std::size_t>(it1->type);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
      [[maybe_unused]] double scalingVDWA = it1->scalingVDW;
      [[maybe_unused]] double scaleCoulombA = it1->scalingCoulomb;
      [[maybe_unused]] double chargeA = it1->charge;

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

    for (std::span<const Atom>::iterator it1 = oldatoms.begin(); it1 != oldatoms.end(); ++it1)
    {
      [[maybe_unused]] std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
      [[maybe_unused]] std::size_t compA = static_cast<std::size_t>(it1->componentId);
      [[maybe_unused]] std::size_t typeA = static_cast<std::size_t>(it1->type);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
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
            energyFactor.energy = externalFieldInterpolationGrid->interpolate(posA);
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
      energySum.externalFieldVDW -= energyFactor.energy;
      energySum.dudlambdaVDW -= energyFactor.dUdlambda;
    }
  }

  return energySum;
}
