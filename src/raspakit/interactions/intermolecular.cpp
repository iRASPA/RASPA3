module;

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
#include <vector>
#endif

module interactions_intermolecular;

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
#endif

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import potential_electrostatics;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import component;
import units;
import threadpool;

// used in volume moves for computing the state at a new box and new, scaled atom positions
RunningEnergy Interactions::computeInterMolecularEnergy(const ForceField &forceField, const SimulationBox &box,
                                                        std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  RunningEnergy energySum{};
  bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (forceField.omitInterInteractions) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
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

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }
  }

  return energySum;
}

RunningEnergy Interactions::computeInterMolecularTailEnergy(const ForceField &forceField,
                                                            const SimulationBox &simulationBox,
                                                            std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    double temp_self = preFactor * forceField(typeA, typeA).tailCorrectionEnergy;
    energySum.tail += scalingVDWA * scalingVDWA * temp_self;
    energySum.dudlambdaVDW += (groupIdA ? scalingVDWA * temp_self : 0.0) + (groupIdA ? scalingVDWA * temp_self : 0.0);

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  return energySum;
}

// used in mc_moves_translation.cpp, mc_moves_rotation.cpp,
//         mc_moves_random_translation.cpp, mc_moves_random_rotation.cpp
//         mc_moves_swap_cfcmc.cpp, mc_moves_swap_cfcmc_cbmc.cpp, mc_moves_gibbs_swap_cfcmc.cpp
[[nodiscard]] std::optional<RunningEnergy> Interactions::computeInterMolecularEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

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

    for (const Atom &atom : newatoms)
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

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }

    for (const Atom &atom : oldatoms)
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

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW -= energyFactor.energy;
          energySum.dudlambdaVDW -= energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge -= energyFactor.energy;
          energySum.dudlambdaCharge -= energyFactor.dUdlambda;
        }
      }
    }
  }

  return std::optional{energySum};
}

[[nodiscard]] RunningEnergy Interactions::computeInterMolecularTailEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (const Atom &atom : newatoms)
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
        energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
      }
    }

    for (const Atom &atom : oldatoms)
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
        energySum.dudlambdaVDW -= (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
      }
    }
  }

  for (const Atom &atomA : newatoms)
  {
    size_t typeA = static_cast<size_t>(atomA.type);
    bool groupIdA = static_cast<bool>(atomA.groupId);
    double scalingVDWA = atomA.scalingVDW;

    for (const Atom &atomB : newatoms)
    {
      size_t typeB = static_cast<size_t>(atomB.type);
      bool groupIdB = static_cast<bool>(atomB.groupId);
      double scalingVDWB = atomB.scalingVDW;

      double temp = preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  for (const Atom &atomA : oldatoms)
  {
    size_t typeA = static_cast<size_t>(atomA.type);
    bool groupIdA = static_cast<bool>(atomA.groupId);
    double scalingVDWA = atomA.scalingVDW;

    for (const Atom &atomB : oldatoms)
    {
      size_t typeB = static_cast<size_t>(atomB.type);
      bool groupIdB = static_cast<bool>(atomB.groupId);
      double scalingVDWB = atomB.scalingVDW;

      double temp = preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail -= scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW -= (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  return energySum;
}

RunningEnergy Interactions::computeInterMolecularGradient(const ForceField &forceField,
                                                          const SimulationBox &simulationBox,
                                                          std::span<Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB;
  double rr;
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (forceField.omitInterInteractions) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
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
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::GradientFactor gradientFactor = Potentials::potentialVDWGradient(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW += gradientFactor.energy;
          energySum.dudlambdaVDW += gradientFactor.dUdlambda;

          const double3 f = gradientFactor.gradientFactor * dr;

          it1->gradient += f;
          it2->gradient -= f;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::GradientFactor gradientFactor = Potentials::potentialCoulombGradient(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += gradientFactor.energy;
          energySum.dudlambdaCharge += gradientFactor.dUdlambda;

          const double3 f = gradientFactor.gradientFactor * dr;

          it1->gradient += f;
          it2->gradient -= f;
        }
      }
    }
  }

  return energySum;
}

std::pair<EnergyStatus, double3x3> Interactions::computeInterMolecularEnergyStrainDerivative(
    const ForceField &forceField, const std::vector<Component> &components, const SimulationBox &simulationBox,
    std::span<Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB;
  double rr;

  bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  EnergyStatus energy(1, 1, components.size());
  double3x3 strainDerivativeTensor{};

  if (forceField.omitInterInteractions) return {energy, strainDerivativeTensor};
  if (moleculeAtoms.empty()) return {energy, strainDerivativeTensor};

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    energy.componentEnergy(compA, compA).VanDerWaalsTailCorrection += Potentials::EnergyFactor(
        preFactor * scalingVDWA * scalingVDWA * forceField(typeA, typeA).tailCorrectionEnergy, 0.0);

    for (std::span<Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);
      size_t typeB = static_cast<size_t>(it2->type);
      double scalingVDWB = it2->scalingVDW;

      Potentials::EnergyFactor temp(
          preFactor * scalingVDWA * scalingVDWB * forceField(typeA, typeB).tailCorrectionEnergy, 0.0);
      energy.componentEnergy(compA, compB).VanDerWaalsTailCorrection += 2.0 * temp;

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        bool groupIdB = static_cast<bool>(it2->groupId);
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::GradientFactor gradientFactor = Potentials::potentialVDWGradient(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energy.componentEnergy(compA, compB).VanDerWaals +=
              0.5 * Potentials::EnergyFactor(gradientFactor.energy, 0.0);
          energy.componentEnergy(compB, compA).VanDerWaals +=
              0.5 * Potentials::EnergyFactor(gradientFactor.energy, 0.0);

          const double3 g = gradientFactor.gradientFactor * dr;

          it1->gradient += g;
          it2->gradient -= g;

          strainDerivativeTensor.ax += g.x * dr.x;
          strainDerivativeTensor.bx += g.y * dr.x;
          strainDerivativeTensor.cx += g.z * dr.x;

          strainDerivativeTensor.ay += g.x * dr.y;
          strainDerivativeTensor.by += g.y * dr.y;
          strainDerivativeTensor.cy += g.z * dr.y;

          strainDerivativeTensor.az += g.x * dr.z;
          strainDerivativeTensor.bz += g.y * dr.z;
          strainDerivativeTensor.cz += g.z * dr.z;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);

          Potentials::GradientFactor energyFactor = Potentials::potentialCoulombGradient(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energy.componentEnergy(compA, compB).CoulombicReal +=
              0.5 * Potentials::EnergyFactor(energyFactor.energy, 0.0);
          energy.componentEnergy(compB, compA).CoulombicReal +=
              0.5 * Potentials::EnergyFactor(energyFactor.energy, 0.0);

          const double3 g = energyFactor.gradientFactor * dr;

          it1->gradient += g;
          it2->gradient -= g;

          strainDerivativeTensor.ax += g.x * dr.x;
          strainDerivativeTensor.bx += g.y * dr.x;
          strainDerivativeTensor.cx += g.z * dr.x;

          strainDerivativeTensor.ay += g.x * dr.y;
          strainDerivativeTensor.by += g.y * dr.y;
          strainDerivativeTensor.cy += g.z * dr.y;

          strainDerivativeTensor.az += g.x * dr.z;
          strainDerivativeTensor.bz += g.y * dr.z;
          strainDerivativeTensor.cz += g.z * dr.z;
        }
      }
    }
  }

  return {energy, strainDerivativeTensor};
}

void Interactions::computeInterMolecularElectricPotential(const ForceField &forceField, const SimulationBox &box,
                                                          std::span<double> electricPotentialMolecules,
                                                          std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool useCharge = forceField.useCharge;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!useCharge) return;
  if (forceField.omitInterInteractions) return;
  if (forceField.omitInterPolarization) return;
  if (moleculeAtoms.empty()) return;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = box.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);

          double potential = Potentials::potentialElectrostatics(forceField, 1.0, r, 1.0);

          size_t indexA = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it1));
          electricPotentialMolecules[indexA] += scalingCoulombB * chargeB * potential;

          size_t indexB = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it2));
          electricPotentialMolecules[indexB] += scalingCoulombA * chargeA * potential;
        }
      }
    }
  }
}
RunningEnergy Interactions::computeInterMolecularElectricField(const ForceField &forceField, const SimulationBox &box,
                                                               std::span<double3> electricFieldMolecules,
                                                               std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!useCharge) return energySum;
  if (forceField.omitInterInteractions) return energySum;
  if (forceField.omitInterPolarization) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
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

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);
          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;

          Potentials::GradientFactor gradient =
              Potentials::potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);

          Potentials::GradientFactor gradientFactorA = scalingCoulombB * chargeB * gradient;
          size_t indexA = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it1));
          electricFieldMolecules[indexA] -= gradientFactorA.gradientFactor * dr;

          Potentials::GradientFactor gradientFactorB = scalingCoulombA * chargeA * gradient;
          size_t indexB = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it2));
          electricFieldMolecules[indexB] += gradientFactorB.gradientFactor * dr;
        }
      }
    }
  }

  return energySum;
}

std::optional<RunningEnergy> Interactions::computeInterMolecularElectricFieldDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<double3> electricFieldMolecules,
    std::span<double3> electricFieldMolecule, std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;
  if (forceField.omitInterPolarization) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    size_t indexA = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it1));
    size_t molA = static_cast<size_t>(it1->moleculeId);
    double3 posA = it1->position;
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      size_t indexB = static_cast<size_t>(std::distance(newatoms.begin(), it2));
      size_t compB = static_cast<size_t>(it2->componentId);
      size_t molB = static_cast<size_t>(it2->moleculeId);

      if (!(compA == compB && molA == molB))
      {
        double3 posB = it2->position;
        size_t typeB = static_cast<size_t>(it2->type);
        bool groupIdB = static_cast<bool>(it2->groupId);
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;

          Potentials::GradientFactor gradient =
              Potentials::potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);

          Potentials::GradientFactor gradientFactorA = scalingCoulombB * chargeB * gradient;
          electricFieldMolecules[indexA] -= gradientFactorA.gradientFactor * dr;

          Potentials::GradientFactor gradientFactorB = scalingCoulombA * chargeA * gradient;
          electricFieldMolecule[indexB] += gradientFactorB.gradientFactor * dr;
        }
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      size_t indexB = static_cast<size_t>(std::distance(oldatoms.begin(), it2));
      size_t compB = static_cast<size_t>(it2->componentId);
      size_t molB = static_cast<size_t>(it2->moleculeId);

      if (!(compA == compB && molA == molB))
      {
        double3 posB = it2->position;
        size_t typeB = static_cast<size_t>(it2->type);
        bool groupIdB = static_cast<bool>(it2->groupId);
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW -= energyFactor.energy;
          energySum.dudlambdaVDW -= energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge -= energyFactor.energy;
          energySum.dudlambdaCharge -= energyFactor.dUdlambda;

          Potentials::GradientFactor gradient =
              Potentials::potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
          ;

          Potentials::GradientFactor gradientFactorA = scalingCoulombB * chargeB * gradient;
          electricFieldMolecules[indexA] += gradientFactorA.gradientFactor * dr;

          Potentials::GradientFactor gradientFactorB = scalingCoulombA * chargeA * gradient;
          electricFieldMolecule[indexB] -= gradientFactorB.gradientFactor * dr;
        }
      }
    }
  }

  return std::optional{energySum};
}
