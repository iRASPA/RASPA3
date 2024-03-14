module;

module interactions_external_field;

import <complex>;
import <span>;
import <numbers>;
import <cmath>;
import <vector>;
import <iostream>;
import <algorithm>;
import <type_traits>;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import force_factor;
import running_energy;
import component;
import forcefield;

void Interactions::computeExternalFieldEnergy(bool hasExternalField, [[maybe_unused]] const ForceField &forceField, 
                                              [[maybe_unused]] const SimulationBox &simulationBox,
                                              [[maybe_unused]] std::span<const Atom> moleculeAtoms, 
                                              [[maybe_unused]] RunningEnergy &energyStatus) noexcept
{
  if(hasExternalField)
  {
    if (moleculeAtoms.empty()) return;

    for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
    {
      [[maybe_unused]] size_t molA = static_cast<size_t>(it1->moleculeId);
      [[maybe_unused]] size_t compA = static_cast<size_t>(it1->componentId);
      [[maybe_unused]] size_t typeA = static_cast<size_t>(it1->type);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
      [[maybe_unused]] double scalingVDWA = it1->scalingVDW;
      [[maybe_unused]] double scaleCoulombA = it1->scalingCoulomb;
      [[maybe_unused]] double chargeA = it1->charge;
      [[maybe_unused]] double3 posA = it1->position;
      [[maybe_unused]] double3 s = (simulationBox.inverseUnitCell * posA).fract();


      // Fill in the energy based on the atom properties and the fractional position 's'
      EnergyFactor energyFactor = EnergyFactor(0.0, 0.0);

      energyStatus.externalFieldVDW += energyFactor.energy;
      energyStatus.dudlambdaVDW += energyFactor.dUdlambda;
    }
  }
}

void Interactions::computeExternalFieldTailEnergy(bool hasExternalField,[[maybe_unused]] const ForceField &forceField, 
                                                  [[maybe_unused]] const SimulationBox &simulationBox,
                                                  [[maybe_unused]] std::span<const Atom> moleculeAtoms, 
                                                  [[maybe_unused]] RunningEnergy &energyStatus) noexcept
{
  if(hasExternalField)
  {
  }
}

[[nodiscard]] std::optional<RunningEnergy>
Interactions::computeExternalFieldEnergyDifference(bool hasExternalField, [[maybe_unused]] const ForceField &forceField, 
                                                   [[maybe_unused]] const SimulationBox &simulationBox,
                                                   [[maybe_unused]] std::span<const Atom> newatoms, 
                                                   [[maybe_unused]] std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum;

  const double overlapCriteria = forceField.overlapCriteria;

  if(hasExternalField)
  {
    for (std::span<const Atom>::iterator it1 = newatoms.begin(); it1 != newatoms.end(); ++it1)
    {
      [[maybe_unused]]  double3 posA = it1->position;
      [[maybe_unused]] size_t molA = static_cast<size_t>(it1->moleculeId);
      [[maybe_unused]] size_t compA = static_cast<size_t>(it1->componentId);
      [[maybe_unused]] size_t typeA = static_cast<size_t>(it1->type);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
      [[maybe_unused]] double scalingVDWA = it1->scalingVDW;
      [[maybe_unused]] double scaleCoulombA = it1->scalingCoulomb;
      [[maybe_unused]] double chargeA = it1->charge;

      // FILL IN
      EnergyFactor energyFactor = EnergyFactor(0.0, 0.0);
      if (energyFactor.energy > overlapCriteria) return std::nullopt;

      energySum.externalFieldVDW += energyFactor.energy;
      energySum.dudlambdaVDW += energyFactor.dUdlambda;
    }

    for (std::span<const Atom>::iterator it1 = oldatoms.begin(); it1 != oldatoms.end(); ++it1)
    {
      [[maybe_unused]] size_t molA = static_cast<size_t>(it1->moleculeId);
      [[maybe_unused]] size_t compA = static_cast<size_t>(it1->componentId);
      [[maybe_unused]] size_t typeA = static_cast<size_t>(it1->type);
      [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
      [[maybe_unused]] double scalingVDWA = it1->scalingVDW;
      [[maybe_unused]] double scaleCoulombA = it1->scalingCoulomb;
      [[maybe_unused]] double chargeA = it1->charge;
      [[maybe_unused]] double3 posA = it1->position;
      [[maybe_unused]] double3 s = (simulationBox.inverseUnitCell * posA).fract();


      // Fill in the energy based on the atom properties and the fractional position 's'
      EnergyFactor energyFactor = EnergyFactor(0.0, 0.0);

      energySum.externalFieldVDW -= energyFactor.energy;
      energySum.dudlambdaVDW -= energyFactor.dUdlambda;
    }
  }

  return energySum;
}
