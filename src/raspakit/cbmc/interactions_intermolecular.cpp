module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <thread>
#include <type_traits>
#include <vector>
#endif

module cbmc_interactions_intermolecular;

#ifndef USE_LEGACY_HEADERS
import <numbers>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <optional>;
import <thread>;
import <future>;
import <type_traits>;
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
import force_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum;

  bool useCharge = forceField.useCharge;

  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    double3 posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (int index = 0; const Atom &atom : atoms)
    {
      if (index != skip)
      {
        double3 posB = atom.position;
        size_t compB = static_cast<size_t>(atom.componentId);
        size_t molB = static_cast<size_t>(atom.moleculeId);
        size_t typeB = static_cast<size_t>(atom.type);
        bool groupIdB = static_cast<bool>(atom.groupId);
        double scalingVDWB = atom.scalingVDW;
        double scalingCoulombB = atom.scalingCoulomb;
        double chargeB = atom.charge;

        if (!(compA == compB && molA == molB))
        {
          dr = posA - posB;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          rr = double3::dot(dr, dr);

          if (rr < cutOffVDWSquared)
          {
            EnergyFactor energyFactor =
                potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
            if (energyFactor.energy > overlapCriteria) return std::nullopt;

            energySum.moleculeMoleculeVDW += energyFactor.energy;
            energySum.dudlambdaVDW += energyFactor.dUdlambda;
          }
          if (useCharge && rr < cutOffChargeSquared)
          {
            double r = std::sqrt(rr);
            EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                               scalingCoulombB, r, chargeA, chargeB);

            energySum.moleculeMoleculeCharge += energyFactor.energy;
            energySum.dudlambdaCharge += energyFactor.dUdlambda;
          }
        }
      }
      ++index;
    }
  }

  return energySum;
}
