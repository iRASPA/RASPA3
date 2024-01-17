module;

module system;

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

import <numbers>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <optional>;
import <thread>;
import <future>;


// used in volume moves for computing the state at a new box and new, scaled atom positions
void System::computeInterMolecularEnergy(const SimulationBox &box, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtomPositions.empty()) return;

  for (std::span<const Atom>::iterator it1 = moleculeAtomPositions.begin(); it1 != moleculeAtomPositions.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtomPositions.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        size_t typeB = static_cast<size_t>(it2->type);
        bool groupIdB = static_cast<bool>(it2->groupId);
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = box.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffVDWSquared)
        {
          EnergyFactor energyFactor = potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energyStatus.moleculeMoleculeVDW += energyFactor.energy;
          energyStatus.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (!noCharges && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energyStatus.moleculeMoleculeCharge += energyFactor.energy;
          energyStatus.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }
  }
}

ForceFactor System::computeInterMolecularGradient() noexcept
{
  double3 dr, posA, posB;
  double rr;

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  ForceFactor energy{ 0.0, 0.0, 0.0 };

  std::span<Atom> moleculeAtoms = spanOfMoleculeAtoms();
  if (moleculeAtoms.empty()) return energy;

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scaleCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    
    for (std::span<Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        size_t typeB = static_cast<size_t>(it2->type);
        bool groupIdB = static_cast<bool>(it2->groupId);
        double scalingVDWB = it2->scalingVDW;
        double scaleCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffVDWSquared)
        {
          ForceFactor forceFactor = potentialVDWGradient(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energy += forceFactor;

          const double3 f = forceFactor.forceFactor * dr;

          it1->gradient += f;
          it2->gradient -= f;
        }
        if (!noCharges && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          ForceFactor forceFactor = potentialCoulombGradient(forceField, groupIdA, groupIdB, scaleCoulombA, scaleCoulombB, r, chargeA, chargeB);

          energy += forceFactor;

          const double3 f = forceFactor.forceFactor * dr;

          it1->gradient += f;
          it2->gradient -= f;
        }
      }
    }
  }

  return energy;
}

// used in mc_moves_translation.cpp, mc_moves_rotation.cpp, mc_moves_random_translation.cpp, mc_moves_random_rotation.cpp
//         mc_moves_swap_cfcmc.cpp, mc_moves_swap_cfcmc_cbmc.cpp, mc_moves_gibbs_swap_cfcmc.cpp
[[nodiscard]] std::optional<RunningEnergy> System::computeInterMolecularEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool noCharges = forceField.noCharges;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    size_t molA = static_cast<size_t>(it1->moleculeId);
    double3 posA = it1->position;
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (const Atom& atom : newatoms)
    {
      size_t compB = static_cast<size_t>(atom.componentId);
      size_t molB = static_cast<size_t>(atom.moleculeId);

      if (!(compA == compB && molA == molB))
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

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (!noCharges && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }

    for (const Atom& atom : oldatoms)
    {
      size_t compB = static_cast<size_t>(atom.componentId);
      size_t molB = static_cast<size_t>(atom.moleculeId);

      if (!(compA == compB && molA == molB))
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

          energySum.moleculeMoleculeVDW -= energyFactor.energy;
          energySum.dudlambdaVDW -= energyFactor.dUdlambda;
        }
        if (!noCharges && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          EnergyFactor energyFactor =  potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge -= energyFactor.energy;
          energySum.dudlambdaCharge -= energyFactor.dUdlambda;
        }
      }
    }
  }

  return energySum;
}

void System::computeInterMolecularTailEnergy(std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept
{
  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::span<const Atom>::iterator it1 = moleculeAtomPositions.begin(); it1 != moleculeAtomPositions.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    [[maybe_unused]] bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtomPositions.end(); ++it2)
    {
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energyStatus.tail += scalingVDWA * scalingVDWB * temp;
      energyStatus.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0) 
                                 + (groupIdB ? scalingVDWA * temp : 0.0);
    }

    double temp = preFactor * forceField(typeA, typeA).tailCorrectionEnergy;
    energyStatus.tail += scalingVDWA * scalingVDWA * temp;
    energyStatus.dudlambdaVDW += (groupIdA ? 2.0 * scalingVDWA * temp : 0.0);
  }
}

[[nodiscard]] RunningEnergy System::computeInterMolecularTailEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (const Atom& atom : newatoms)
    {
      size_t compB = static_cast<size_t>(atom.componentId);
      size_t molB = static_cast<size_t>(atom.moleculeId);

      if (!(compA == compB && molA == molB))
      {
        size_t typeB = static_cast<size_t>(atom.type);
        bool groupIdB = static_cast<bool>(atom.groupId);
        double scalingVDWB = atom.scalingVDW;

        double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
        energySum.tail += scalingVDWA * scalingVDWB * temp;
        energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0)
                                 + (groupIdB ? scalingVDWA * temp : 0.0);
      }
    }

    for (const Atom& atom : oldatoms)
    {
      size_t compB = static_cast<size_t>(atom.componentId);
      size_t molB = static_cast<size_t>(atom.moleculeId);

      if (!(compA == compB && molA == molB))
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
  }

  for (const Atom& atomA : newatoms)
  {
    size_t typeA = static_cast<size_t>(atomA.type);
    bool groupIdA = static_cast<bool>(atomA.groupId);
    double scalingVDWA = atomA.scalingVDW;

    for (const Atom& atomB : newatoms)
    {
      size_t typeB = static_cast<size_t>(atomB.type);
      bool groupIdB = static_cast<bool>(atomB.groupId);
      double scalingVDWB = atomB.scalingVDW;

      double temp = preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0)
                              + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  for (const Atom& atomA : oldatoms)
  {
    size_t typeA = static_cast<size_t>(atomA.type);
    bool groupIdA = static_cast<bool>(atomA.groupId);
    double scalingVDWA = atomA.scalingVDW;

    for (const Atom& atomB : oldatoms)
    {
      size_t typeB = static_cast<size_t>(atomB.type);
      bool groupIdB = static_cast<bool>(atomB.groupId);
      double scalingVDWB = atomB.scalingVDW;

      double temp = preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail -= scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW -= (groupIdA ? scalingVDWB * temp : 0.0)
                              + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  return energySum;
}
