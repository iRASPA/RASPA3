module;

module system;

import energy_status;
import potential_energy_vdw;
import potential_energy_coulomb;
import potential_gradient_vdw;
import potential_gradient_coulomb;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
import threading;
import energy_factor;
import force_factor;

import <numbers>;
import <optional>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <thread>;
import <future>;
import <deque>;
import <semaphore>;
import <atomic>;

// system_interactions_framework_molecule.cpp

// Used to compute the total framework-molecule energy
//

void System::computeFrameworkMoleculeEnergy(const SimulationBox &box, std::span<const Atom> frameworkAtomPositions, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtomPositions.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtomPositions.begin(); it1 != frameworkAtomPositions.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scaleCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    for (std::span<const Atom>::iterator it2 = moleculeAtomPositions.begin(); it2 != moleculeAtomPositions.end(); ++it2)
    {
      posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;
      double scaleCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = box.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffVDWSquared)
      {
        EnergyFactor energyFactor = potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energyStatus.frameworkMoleculeVDW += energyFactor.energy;
        energyStatus.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scaleCoulombA, scaleCoulombB, r, chargeA, chargeB);

        energyStatus.frameworkMoleculeCharge += energyFactor.energy;
        energyStatus.dudlambdaCharge += energyFactor.dUdlambda;
      }
    }
  }
}

ForceFactor System::computeFrameworkMoleculeGradient() noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  std::span<Atom> moleculeAtoms = spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

  ForceFactor energy{ 0.0, 0.0, 0.0 };

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    for (std::span<Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffVDWSquared)
      {
        ForceFactor forceFactor = potentialVDWGradient(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energy += forceFactor;

        const double3 factor = forceFactor.forceFactor * dr;

        it2->gradient -= factor;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        ForceFactor energyFactor = potentialCoulombGradient(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energy += energyFactor;

        const double3 factor = energyFactor.forceFactor * dr;

        it2->gradient -= factor;
      }
    }
  }
  return energy;
}

// Used in Translation and Rotation
//

[[nodiscard]] std::optional<RunningEnergy> System::computeFrameworkMoleculeEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool noCharges = forceField.noCharges;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (const Atom& atom : newatoms)
    {
      double3 posB = atom.position;
      size_t typeB = static_cast<size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      double scalingVDWB = atom.scalingVDW;
      double scalingCoulombB = atom.scalingCoulomb;
      double chargeB = atom.charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffVDWSquared)
      {
        EnergyFactor energyFactor = potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW +=  energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;
      }
    }

    for (const Atom& atom : oldatoms)
    {
      double3 posB = atom.position;
      size_t typeB = static_cast<size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      double scalingVDWB = atom.scalingVDW;
      double scalingCoulombB = atom.scalingCoulomb;
      double chargeB = atom.charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffVDWSquared)
      {
        EnergyFactor energyFactor = potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.dudlambdaVDW -= energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.dudlambdaCharge -= energyFactor.dUdlambda;
      }
    }
  }

  return energySum;
}

void System::computeFrameworkMoleculeTailEnergy(std::span<const Atom> frameworkAtomPositions, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept
{
  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::span<const Atom>::iterator it1 = frameworkAtomPositions.begin(); it1 != frameworkAtomPositions.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (std::span<const Atom>::iterator it2 = moleculeAtomPositions.begin(); it2 != moleculeAtomPositions.end(); ++it2)
    {
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energyStatus.tail += scalingVDWA * scalingVDWB * temp;
      energyStatus.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0)
                                 + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }
}

[[nodiscard]] RunningEnergy System::computeFrameworkMoleculeTailEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (const Atom& atom : newatoms)
    {
      size_t typeB = static_cast<size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      double scalingVDWB = atom.scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0)
                              + (groupIdB ? scalingVDWA * temp : 0.0);
    }

    for (const Atom& atom : oldatoms)
    {
      size_t typeB = static_cast<size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      double scalingVDWB = atom.scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail -= scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW -= (groupIdA ? scalingVDWB * temp : 0.0)
                              + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  return energySum;
}

