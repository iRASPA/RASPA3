module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <deque>
#include <future>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <semaphore>
#include <span>
#include <thread>
#include <utility>
#include <vector>
#endif

module interactions_framework_molecule;

#ifndef USE_LEGACY_HEADERS
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
import <utility>;
import <limits>;
#endif

import double3;
import double4;
import double3x3;
import double3x3x3;
import energy_status;
import potential_energy_vdw;
import potential_energy_coulomb;
import potential_gradient_vdw;
import potential_gradient_coulomb;
import potential_hessian_vdw;
import potential_hessian_coulomb;
import potential_tricubic_derivative_lj;
import potential_tricubic_derivative_real_ewald;
import potential_triquintic_derivative_lj;
import potential_electrostatics;
import simulationbox;
import framework;
import forcefield;
import atom;
import energy_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
// import threading;
import energy_factor;
import gradient_factor;
import hessian_factor;
import tricubic_derivative_factor;
import triquintic_derivative_factor;
import framework;
import component;
import interpolation_energy_grid;

RunningEnergy Interactions::computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!framework.has_value()) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
  {
    posB = it2->position;
    size_t typeB = static_cast<size_t>(it2->type);
    bool groupIdB = static_cast<bool>(it2->groupId);
    double scalingVDWB = it2->scalingVDW;
    double scaleCoulombB = it2->scalingCoulomb;
    double chargeB = it2->charge;

    if (interpolationGrids[typeB].has_value() && (groupIdB == 0))
    {
      energySum.frameworkMoleculeVDW += interpolationGrids[typeB]->interpolate(posB);
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge += chargeB * interpolationGrids.back()->interpolate(posB);
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        posA = it1->position;
        size_t typeA = static_cast<size_t>(it1->type);
        bool groupIdA = static_cast<bool>(it1->groupId);
        double scalingVDWA = it1->scalingVDW;
        double scaleCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scaleCoulombA, scaleCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }
  }
  return energySum;
}

RunningEnergy Interactions::computeFrameworkMoleculeTailEnergy(const ForceField &forceField,
                                                               const SimulationBox &simulationBox,
                                                               std::span<const Atom> frameworkAtoms,
                                                               std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
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

// Used in Translation and Rotation
//

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (!framework.has_value()) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  std::vector<uint8_t> forceExplicit(std::max(newatoms.size(), oldatoms.size()));
  for (std::size_t i = 0; i < newatoms.size(); ++i)
  {
    forceExplicit[i] = newatoms[i].groupId;
  }
  for (std::size_t i = 0; i < oldatoms.size(); ++i)
  {
    forceExplicit[i] = forceExplicit[i] || oldatoms[i].groupId;
  }

  for (size_t i = 0; i < newatoms.size(); ++i)
  {
    const Atom &atom = newatoms[i];
    double3 posB = atom.position;
    size_t typeB = static_cast<size_t>(atom.type);
    bool groupIdB = static_cast<bool>(atom.groupId);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !forceExplicit[i])
    {
      double energy = interpolationGrids[typeB]->interpolate(posB);
      if (energy > overlapCriteria)
      {
        return std::nullopt;
      }
      energySum.frameworkMoleculeVDW += energy;
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge += chargeB * interpolationGrids.back()->interpolate(posB);
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posA = it1->position;
        size_t typeA = static_cast<size_t>(it1->type);
        bool groupIdA = static_cast<bool>(it1->groupId);
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }
  }

  for (size_t i = 0; i < oldatoms.size(); ++i)
  {
    const Atom &atom = oldatoms[i];
    double3 posB = atom.position;
    size_t typeB = static_cast<size_t>(atom.type);
    bool groupIdB = static_cast<bool>(atom.groupId);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !forceExplicit[i])
    {
      energySum.frameworkMoleculeVDW -= interpolationGrids[typeB]->interpolate(posB);
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge -= chargeB * interpolationGrids.back()->interpolate(posB);
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posA = it1->position;
        size_t typeA = static_cast<size_t>(it1->type);
        bool groupIdA = static_cast<bool>(it1->groupId);
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW -= energyFactor.energy;
          energySum.dudlambdaVDW -= energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge -= energyFactor.energy;
          energySum.dudlambdaCharge -= energyFactor.dUdlambda;
        }
      }
    }
  }

  return std::optional{energySum};
}

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifferenceOriginal(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (!framework.has_value()) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (auto &atom : newatoms)
  {
    double3 posB = atom.position;
    size_t typeB = static_cast<size_t>(atom.type);
    bool groupIdB = static_cast<bool>(atom.groupId);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !groupIdB)
    {
      double energy = interpolationGrids[typeB]->interpolate(posB);
      if (energy > overlapCriteria)
      {
        return std::nullopt;
      }
      energySum.frameworkMoleculeVDW += energy;
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge += chargeB * interpolationGrids.back()->interpolate(posB);
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posA = it1->position;
        size_t typeA = static_cast<size_t>(it1->type);
        bool groupIdA = static_cast<bool>(it1->groupId);
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }
  }

  for (auto &atom : oldatoms)
  {
    double3 posB = atom.position;
    size_t typeB = static_cast<size_t>(atom.type);
    bool groupIdB = static_cast<bool>(atom.groupId);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !groupIdB)
    {
      energySum.frameworkMoleculeVDW -= interpolationGrids[typeB]->interpolate(posB);
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge -= chargeB * interpolationGrids.back()->interpolate(posB);
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posA = it1->position;
        size_t typeA = static_cast<size_t>(it1->type);
        bool groupIdA = static_cast<bool>(it1->groupId);
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW -= energyFactor.energy;
          energySum.dudlambdaVDW -= energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge -= energyFactor.energy;
          energySum.dudlambdaCharge -= energyFactor.dUdlambda;
        }
      }
    }
  }

  return std::optional{energySum};
}

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifferenceInterpolationExplicit(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> atoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!framework.has_value()) return energySum;
  if (atoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it2 = atoms.begin(); it2 != atoms.end(); ++it2)
  {
    posB = it2->position;
    size_t typeB = static_cast<size_t>(it2->type);
    bool groupIdB = static_cast<bool>(it2->groupId);
    double scalingVDWB = it2->scalingVDW;
    double scaleCoulombB = it2->scalingCoulomb;
    double chargeB = it2->charge;

    if (interpolationGrids[typeB].has_value() && !groupIdB)
    {
      energySum.frameworkMoleculeVDW -= interpolationGrids[typeB]->interpolate(posB);
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge -= chargeB * interpolationGrids.back()->interpolate(posB);
      }

      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        posA = it1->position;
        size_t typeA = static_cast<size_t>(it1->type);
        bool groupIdA = static_cast<bool>(it1->groupId);
        double scalingVDWA = it1->scalingVDW;
        double scaleCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.dudlambdaVDW += energyFactor.dUdlambda;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
              forceField, groupIdA, groupIdB, scaleCoulombA, scaleCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.dudlambdaCharge += energyFactor.dUdlambda;
        }
      }
    }
  }
  return energySum;
}

[[nodiscard]] RunningEnergy Interactions::computeFrameworkMoleculeTailEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;

    for (const Atom &atom : newatoms)
    {
      size_t typeB = static_cast<size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      double scalingVDWB = atom.scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
    }

    for (const Atom &atom : oldatoms)
    {
      size_t typeB = static_cast<size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      double scalingVDWB = atom.scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail -= scalingVDWA * scalingVDWB * temp;
      energySum.dudlambdaVDW -= (groupIdA ? scalingVDWB * temp : 0.0) + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }

  return energySum;
}

RunningEnergy Interactions::computeFrameworkMoleculeGradient(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
    std::span<Atom> moleculeAtoms,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids) noexcept
{
  RunningEnergy energySum{};

  double3 dr, posA, posB;
  double rr;

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtoms.empty()) return energySum;

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    if (interpolationGrids[typeA].has_value() && (groupIdA == 0))
    {
      auto [energy_vdw, gradient_vdw] = interpolationGrids[typeA]->interpolateGradient(posA);
      energySum.frameworkMoleculeVDW += energy_vdw;
      it1->gradient += gradient_vdw;
      if (useCharge)
      {
        auto [energy_real_ewald, gradient_real_ewald] = interpolationGrids.back()->interpolateGradient(posA);
        energySum.frameworkMoleculeCharge += chargeA * energy_real_ewald;
        it1->gradient += chargeA * gradient_real_ewald;
      }
    }
    else
    {
      for (std::span<Atom>::iterator it2 = frameworkAtoms.begin(); it2 != frameworkAtoms.end(); ++it2)
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

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::GradientFactor gradientFactor = Potentials::potentialVDWGradient(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW += gradientFactor.energy;
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

          energySum.frameworkMoleculeCharge += gradientFactor.energy;
          energySum.dudlambdaCharge += gradientFactor.dUdlambda;

          const double3 g = gradientFactor.gradientFactor * dr;

          it1->gradient += g;
          it2->gradient -= g;
        }
      }
    }
  }
  return energySum;
}

[[nodiscard]] std::pair<EnergyStatus, double3x3> Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
    const ForceField &forceField, const std::optional<Framework> &framework,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::vector<Component> &components, const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
    std::span<Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB;
  double rr;

  double3x3 strainDerivativeTensor;
  EnergyStatus energy(1, framework.has_value() ? 1 : 0, components.size());

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  if (moleculeAtoms.empty()) return std::make_pair(energy, strainDerivativeTensor);

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    if (interpolationGrids[typeA].has_value() && (groupIdA == 0))
    {
      auto [energy_vdw, gradient_vdw] = interpolationGrids[typeA]->interpolateGradient(posA);
      energy.frameworkComponentEnergy(compA, 0).VanDerWaals += Potentials::EnergyFactor(energy_vdw, 0.0);
      const double3 f = gradient_vdw;

      it1->gradient += f;

      strainDerivativeTensor.ax += f.x * posA.x;
      strainDerivativeTensor.bx += f.y * posA.x;
      strainDerivativeTensor.cx += f.z * posA.x;

      strainDerivativeTensor.ay += f.x * posA.y;
      strainDerivativeTensor.by += f.y * posA.y;
      strainDerivativeTensor.cy += f.z * posA.y;

      strainDerivativeTensor.az += f.x * posA.z;
      strainDerivativeTensor.bz += f.y * posA.z;
      strainDerivativeTensor.cz += f.z * posA.z;

      if (useCharge)
      {
        auto [energy_real_ewald, gradient_real_ewald] = interpolationGrids.back()->interpolateGradient(posA);
        energy.frameworkComponentEnergy(compA, 0).CoulombicReal +=
            Potentials::EnergyFactor(chargeA * energy_real_ewald, 0.0);
        const double3 g = chargeA * gradient_real_ewald;

        it1->gradient += g;

        strainDerivativeTensor.ax += g.x * posA.x;
        strainDerivativeTensor.bx += g.y * posA.x;
        strainDerivativeTensor.cx += g.z * posA.x;

        strainDerivativeTensor.ay += g.x * posA.y;
        strainDerivativeTensor.by += g.y * posA.y;
        strainDerivativeTensor.cy += g.z * posA.y;

        strainDerivativeTensor.az += g.x * posA.z;
        strainDerivativeTensor.bz += g.y * posA.z;
        strainDerivativeTensor.cz += g.z * posA.z;
      }
    }
    else
    {
      for (std::span<Atom>::iterator it2 = frameworkAtoms.begin(); it2 != frameworkAtoms.end(); ++it2)
      {
        posB = it2->position;
        size_t compB = static_cast<size_t>(it2->componentId);
        size_t typeB = static_cast<size_t>(it2->type);
        bool groupIdB = static_cast<bool>(it2->groupId);
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        Potentials::EnergyFactor temp(
            preFactor * scalingVDWA * scalingVDWB * forceField(typeB, typeA).tailCorrectionEnergy, 0.0);
        energy.frameworkComponentEnergy(compB, compA).VanDerWaalsTailCorrection += 2.0 * temp;

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::GradientFactor gradientFactor = Potentials::potentialVDWGradient(
              forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energy.frameworkComponentEnergy(compB, compA).VanDerWaals +=
              Potentials::EnergyFactor(gradientFactor.energy, 0.0);

          const double3 f = gradientFactor.gradientFactor * dr;

          it1->gradient += f;
          it2->gradient -= f;

          strainDerivativeTensor.ax += f.x * dr.x;
          strainDerivativeTensor.bx += f.y * dr.x;
          strainDerivativeTensor.cx += f.z * dr.x;

          strainDerivativeTensor.ay += f.x * dr.y;
          strainDerivativeTensor.by += f.y * dr.y;
          strainDerivativeTensor.cy += f.z * dr.y;

          strainDerivativeTensor.az += f.x * dr.z;
          strainDerivativeTensor.bz += f.y * dr.z;
          strainDerivativeTensor.cz += f.z * dr.z;
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);

          Potentials::GradientFactor gradientFactor = Potentials::potentialCoulombGradient(
              forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energy.frameworkComponentEnergy(compB, compA).CoulombicReal +=
              Potentials::EnergyFactor(gradientFactor.energy, 0.0);

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
      }
    }
  }

  return std::make_pair(energy, strainDerivativeTensor);
}

void Interactions::computeFrameworkMoleculeElectricPotential(const ForceField &forceField,
                                                             const SimulationBox &simulationBox,
                                                             std::span<double> electricPotentialMolecules,
                                                             std::span<const Atom> frameworkAtoms,
                                                             std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool useCharge = forceField.useCharge;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!useCharge) return;
  if (moleculeAtoms.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        size_t index = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it2));
        electricPotentialMolecules[index] +=
            2.0 * Potentials::potentialElectrostatics(forceField, scalingCoulombA, r, chargeA);
      }
    }
  }
}

RunningEnergy Interactions::computeFrameworkMoleculeElectricField(const ForceField &forceField,
                                                                  const SimulationBox &simulationBox,
                                                                  std::span<double3> electricFieldMolecules,
                                                                  std::span<const Atom> frameworkAtoms,
                                                                  std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!useCharge) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
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

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
            forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA *
            Potentials::potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
        size_t index = static_cast<size_t>(std::distance(moleculeAtoms.begin(), it2));
        electricFieldMolecules[index] += 2.0 * gradientFactor.gradientFactor * dr;
      }
    }
  }

  return energySum;
}

std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeElectricFieldDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<double3> electricFieldMolecule, std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      size_t indexB = static_cast<size_t>(std::distance(newatoms.begin(), it2));
      double3 posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
            forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA *
            Potentials::potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMolecule[indexB] += 2.0 * gradientFactor.gradientFactor * dr;
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      size_t indexB = static_cast<size_t>(std::distance(oldatoms.begin(), it2));
      double3 posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.dudlambdaVDW -= energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
            forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.dudlambdaCharge -= energyFactor.dUdlambda;

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA *
            Potentials::potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMolecule[indexB] -= 2.0 * gradientFactor.gradientFactor * dr;
      }
    }
  }

  return std::optional{energySum};
}

std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionVDW(const ForceField &forceField,
                                                                                   const SimulationBox &simulationBox,
                                                                                   double3 posA, size_t typeA,
                                                                                   std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    size_t typeB = static_cast<size_t>(it1->type);
    bool groupIdB = static_cast<bool>(it1->groupId);
    double scalingB = it1->scalingVDW;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < cutOffFrameworkVDWSquared)
    {
      Potentials::HessianFactor v =
          Potentials::potentialVDWHessian(forceField, 0, groupIdB, 1.0, scalingB, rr, typeA, typeB);

      energy += v.energy;

      first_derivative.x += dr.x * v.firstDerivativeFactor;
      first_derivative.y += dr.y * v.firstDerivativeFactor;
      first_derivative.z += dr.z * v.firstDerivativeFactor;

      // add contribution to the second derivatives (Hessian matrix)
      second_derivative.ax += v.secondDerivativeFactor * dr.x * dr.x + v.firstDerivativeFactor;
      second_derivative.ay += v.secondDerivativeFactor * dr.x * dr.y;
      second_derivative.az += v.secondDerivativeFactor * dr.x * dr.z;

      second_derivative.bx += v.secondDerivativeFactor * dr.y * dr.x;
      second_derivative.by += v.secondDerivativeFactor * dr.y * dr.y + v.firstDerivativeFactor;
      second_derivative.bz += v.secondDerivativeFactor * dr.y * dr.z;

      second_derivative.cx += v.secondDerivativeFactor * dr.z * dr.x;
      second_derivative.cy += v.secondDerivativeFactor * dr.z * dr.y;
      second_derivative.cz += v.secondDerivativeFactor * dr.z * dr.z + v.firstDerivativeFactor;
    }
  }

  return {energy, first_derivative, second_derivative};
}

std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionCoulomb(
    const ForceField &forceField, const SimulationBox &simulationBox, double3 posA, double chargeA,
    std::span<const Atom> frameworkAtoms)
{
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    bool groupIdB = static_cast<bool>(it1->groupId);
    double scalingB = it1->scalingVDW;
    double chargeB = it1->charge;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < cutOffChargeSquared)
    {
      double r = std::sqrt(rr);
      Potentials::HessianFactor v =
          Potentials::potentialCoulombHessian(forceField, 0, groupIdB, 1.0, scalingB, rr, r, chargeA, chargeB);

      energy += v.energy;

      first_derivative.x += dr.x * v.firstDerivativeFactor;
      first_derivative.y += dr.y * v.firstDerivativeFactor;
      first_derivative.z += dr.z * v.firstDerivativeFactor;

      // add contribution to the second derivatives (Hessian matrix)
      second_derivative.ax += v.secondDerivativeFactor * dr.x * dr.x + v.firstDerivativeFactor;
      second_derivative.ay += v.secondDerivativeFactor * dr.x * dr.y;
      second_derivative.az += v.secondDerivativeFactor * dr.x * dr.z;

      second_derivative.bx += v.secondDerivativeFactor * dr.y * dr.x;
      second_derivative.by += v.secondDerivativeFactor * dr.y * dr.y + v.firstDerivativeFactor;
      second_derivative.bz += v.secondDerivativeFactor * dr.y * dr.z;

      second_derivative.cx += v.secondDerivativeFactor * dr.z * dr.x;
      second_derivative.cy += v.secondDerivativeFactor * dr.z * dr.y;
      second_derivative.cz += v.secondDerivativeFactor * dr.z * dr.z + v.firstDerivativeFactor;
    }
  }

  return {energy, first_derivative, second_derivative};
}
