module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <atomic>
#include <cmath>
#include <deque>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <semaphore>
#include <span>
#include <thread>
#include <utility>
#include <vector>
#include <limits>
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
import potential_third_derivative_vdw;
import potential_third_derivative_coulomb;
import potential_electrostatics;
import simulationbox;
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
import third_derivative_factor;
import framework;
import component;

RunningEnergy Interactions::computeFrameworkMoleculeEnergy(const ForceField &forceField,
                                                           const SimulationBox &simulationBox,
                                                           std::span<const Atom> frameworkAtoms,
                                                           std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scaleCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
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

      if (rr < cutOffFrameworkVDWSquared)
      {
        EnergyFactor energyFactor =
            potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor =
            potentialCoulombEnergy(forceField, groupIdA, groupIdB, scaleCoulombA, scaleCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;
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
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
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

    for (const Atom &atom : newatoms)
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

      if (rr < cutOffFrameworkVDWSquared)
      {
        EnergyFactor energyFactor =
            potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                           scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;
      }
    }

    for (const Atom &atom : oldatoms)
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

      if (rr < cutOffFrameworkVDWSquared)
      {
        EnergyFactor energyFactor =
            potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.dudlambdaVDW -= energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                           scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.dudlambdaCharge -= energyFactor.dUdlambda;
      }
    }
  }

  return std::optional{energySum};
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

RunningEnergy Interactions::computeFrameworkMoleculeGradient(const ForceField &forceField,
                                                             const SimulationBox &simulationBox,
                                                             std::span<Atom> frameworkAtoms,
                                                             std::span<Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  double3 dr, posA, posB;
  double rr;

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtoms.empty()) return energySum;

  for (std::span<Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
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

      if (rr < cutOffFrameworkVDWSquared)
      {
        GradientFactor gradientFactor =
            potentialVDWGradient(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW += gradientFactor.energy;
        energySum.dudlambdaVDW += gradientFactor.dUdlambda;

        const double3 f = gradientFactor.gradientFactor * dr;

        it1->gradient += f;
        it2->gradient -= f;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        GradientFactor gradientFactor = potentialCoulombGradient(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                              scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += gradientFactor.energy;
        energySum.dudlambdaCharge += gradientFactor.dUdlambda;

        const double3 f = gradientFactor.gradientFactor * dr;

        it1->gradient += f;
        it2->gradient -= f;
      }
    }
  }
  return energySum;
}

[[nodiscard]] std::pair<EnergyStatus, double3x3> Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
    const ForceField &forceField, const std::vector<Framework> &frameworkComponents,
    const std::vector<Component> &components, const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
    std::span<Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  double3x3 strainDerivativeTensor;
  EnergyStatus energy(1, frameworkComponents.size(), components.size());

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  if (moleculeAtoms.empty()) return std::make_pair(energy, strainDerivativeTensor);

  for (std::span<Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    bool groupIdA = static_cast<bool>(it1->groupId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;
    for (std::span<Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      size_t compB = static_cast<size_t>(it2->componentId);

      posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);
      bool groupIdB = static_cast<bool>(it2->groupId);
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      EnergyFactor temp(preFactor * scalingVDWA * scalingVDWB * forceField(typeA, typeB).tailCorrectionEnergy, 0.0);
      energy.frameworkComponentEnergy(compA, compB).VanDerWaalsTailCorrection += 2.0 * temp;

      if (rr < cutOffFrameworkVDWSquared)
      {
        GradientFactor gradientFactor =
            potentialVDWGradient(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energy.frameworkComponentEnergy(compA, compB).VanDerWaals += EnergyFactor(gradientFactor.energy, 0.0);

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

        GradientFactor gradientFactor = potentialCoulombGradient(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                           scalingCoulombB, r, chargeA, chargeB);

        energy.frameworkComponentEnergy(compA, compB).CoulombicReal += EnergyFactor(gradientFactor.energy, 0.0);

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
        electricPotentialMolecules[index] += 2.0 * potentialElectrostatics(forceField, scalingCoulombA, r, chargeA);
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
        EnergyFactor energyFactor =
            potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                           scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;

        GradientFactor gradientFactor =
            scalingCoulombA * chargeA * potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
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
        EnergyFactor energyFactor =
            potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                           scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dudlambdaCharge += energyFactor.dUdlambda;

        GradientFactor gradientFactor =
            scalingCoulombA * chargeA * potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
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
        EnergyFactor energyFactor =
            potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.dudlambdaVDW -= energyFactor.dUdlambda;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA,
                                                           scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.dudlambdaCharge -= energyFactor.dUdlambda;

        GradientFactor gradientFactor =
            scalingCoulombA * chargeA * potentialCoulombGradient(forceField, groupIdA, groupIdB, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMolecule[indexB] -= 2.0 * gradientFactor.gradientFactor * dr;
      }
    }
  }

  return std::optional{energySum};
}


std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionVDW(const ForceField &forceField, const SimulationBox &simulationBox,
                                                                                   double3 posA, size_t typeA, std::span<const Atom> frameworkAtoms)
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
      HessianFactor v = potentialVDWHessian(forceField, 0, groupIdB, 1.0, scalingB, rr, typeA, typeB);

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


std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionCoulomb(const ForceField &forceField, const SimulationBox &simulationBox,
                                                                                       double3 posA, double chargeA, std::span<const Atom> frameworkAtoms)
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
      HessianFactor v = potentialCoulombHessian(forceField, 0, groupIdB, 1.0, scalingB, rr, r, chargeA, chargeB);

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


std::tuple<double, double3, double3x3, double3x3x3> 
  Interactions::calculateThirdDerivativeAtPositionVDW(const ForceField &forceField, const SimulationBox &simulationBox,
                                                      double3 posA, size_t typeA, std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};
  double3x3x3 third_derivative{};

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
      ThirdDerivativeFactor v = potentialVDWThirdDerivative(forceField, 0, groupIdB, 1.0, scalingB, rr, typeA, typeB);

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

      third_derivative.m111 += v.thirdDerivativeFactor * dr.x * dr.x * dr.x + 3.0 * v.secondDerivativeFactor * dr.x;
      third_derivative.m112 += v.thirdDerivativeFactor * dr.x * dr.x * dr.y + v.secondDerivativeFactor * dr.y;
      third_derivative.m113 += v.thirdDerivativeFactor * dr.x * dr.x * dr.z + v.secondDerivativeFactor * dr.z;

      third_derivative.m121 += v.thirdDerivativeFactor * dr.x * dr.y * dr.x + v.secondDerivativeFactor * dr.y;
      third_derivative.m122 += v.thirdDerivativeFactor * dr.x * dr.y * dr.y + v.secondDerivativeFactor * dr.x;
      third_derivative.m123 += v.thirdDerivativeFactor * dr.x * dr.y * dr.z;

      third_derivative.m131 += v.thirdDerivativeFactor * dr.x * dr.z * dr.x + v.secondDerivativeFactor * dr.z;
      third_derivative.m132 += v.thirdDerivativeFactor * dr.x * dr.z * dr.y;
      third_derivative.m133 += v.thirdDerivativeFactor * dr.x * dr.z * dr.z + v.secondDerivativeFactor * dr.x;

      third_derivative.m211 += v.thirdDerivativeFactor * dr.y * dr.x * dr.x + v.secondDerivativeFactor * dr.y;
      third_derivative.m212 += v.thirdDerivativeFactor * dr.y * dr.x * dr.y + v.secondDerivativeFactor * dr.x;
      third_derivative.m213 += v.thirdDerivativeFactor * dr.y * dr.x * dr.z;

      third_derivative.m221 += v.thirdDerivativeFactor * dr.y * dr.y * dr.x + v.secondDerivativeFactor * dr.x;
      third_derivative.m222 += v.thirdDerivativeFactor * dr.y * dr.y * dr.y + 3.0 * v.secondDerivativeFactor * dr.y;
      third_derivative.m223 += v.thirdDerivativeFactor * dr.y * dr.y * dr.z + v.secondDerivativeFactor * dr.z;

      third_derivative.m231 += v.thirdDerivativeFactor * dr.y * dr.z * dr.x;
      third_derivative.m232 += v.thirdDerivativeFactor * dr.y * dr.z * dr.y + v.secondDerivativeFactor * dr.z;
      third_derivative.m233 += v.thirdDerivativeFactor * dr.y * dr.z * dr.z + v.secondDerivativeFactor * dr.y;

      third_derivative.m311 += v.thirdDerivativeFactor * dr.z * dr.x * dr.x + v.secondDerivativeFactor * dr.z;
      third_derivative.m312 += v.thirdDerivativeFactor * dr.z * dr.x * dr.y;
      third_derivative.m313 += v.thirdDerivativeFactor * dr.z * dr.x * dr.z + v.secondDerivativeFactor * dr.x;

      third_derivative.m321 += v.thirdDerivativeFactor * dr.z * dr.y * dr.x;
      third_derivative.m322 += v.thirdDerivativeFactor * dr.z * dr.y * dr.y + v.secondDerivativeFactor * dr.z;
      third_derivative.m323 += v.thirdDerivativeFactor * dr.z * dr.y * dr.z + v.secondDerivativeFactor * dr.y;

      third_derivative.m331 += v.thirdDerivativeFactor * dr.z * dr.z * dr.x + v.secondDerivativeFactor * dr.x;
      third_derivative.m332 += v.thirdDerivativeFactor * dr.z * dr.z * dr.y + v.secondDerivativeFactor * dr.y;
      third_derivative.m333 += v.thirdDerivativeFactor * dr.z * dr.z * dr.z + 3.0 * v.secondDerivativeFactor * dr.z;
    }
  }

  return {energy, first_derivative, second_derivative, third_derivative};
}


std::tuple<double, double3, double3x3, double3x3x3> 
  Interactions::calculateThirdDerivativeAtPositionCoulomb(const ForceField &forceField, const SimulationBox &simulationBox,
                                                          double3 posA, double chargeA, std::span<const Atom> frameworkAtoms)
{
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};
  double3x3x3 third_derivative{};

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

      ThirdDerivativeFactor v = potentialCoulombThirdDerivative(forceField, 0, groupIdB, 1.0, scalingB, rr, r, chargeA, chargeB);

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

      third_derivative.m111 += v.thirdDerivativeFactor * dr.x * dr.x * dr.x + 3.0 * v.secondDerivativeFactor * dr.x;
      third_derivative.m112 += v.thirdDerivativeFactor * dr.x * dr.x * dr.y + v.secondDerivativeFactor * dr.y;
      third_derivative.m113 += v.thirdDerivativeFactor * dr.x * dr.x * dr.z + v.secondDerivativeFactor * dr.z;

      third_derivative.m121 += v.thirdDerivativeFactor * dr.x * dr.y * dr.x + v.secondDerivativeFactor * dr.y;
      third_derivative.m122 += v.thirdDerivativeFactor * dr.x * dr.y * dr.y + v.secondDerivativeFactor * dr.x;
      third_derivative.m123 += v.thirdDerivativeFactor * dr.x * dr.y * dr.z;

      third_derivative.m131 += v.thirdDerivativeFactor * dr.x * dr.z * dr.x + v.secondDerivativeFactor * dr.z;
      third_derivative.m132 += v.thirdDerivativeFactor * dr.x * dr.z * dr.y;
      third_derivative.m133 += v.thirdDerivativeFactor * dr.x * dr.z * dr.z + v.secondDerivativeFactor * dr.x;

      third_derivative.m211 += v.thirdDerivativeFactor * dr.y * dr.x * dr.x + v.secondDerivativeFactor * dr.y;
      third_derivative.m212 += v.thirdDerivativeFactor * dr.y * dr.x * dr.y + v.secondDerivativeFactor * dr.x;
      third_derivative.m213 += v.thirdDerivativeFactor * dr.y * dr.x * dr.z;

      third_derivative.m221 += v.thirdDerivativeFactor * dr.y * dr.y * dr.x + v.secondDerivativeFactor * dr.x;
      third_derivative.m222 += v.thirdDerivativeFactor * dr.y * dr.y * dr.y + 3.0 * v.secondDerivativeFactor * dr.y;
      third_derivative.m223 += v.thirdDerivativeFactor * dr.y * dr.y * dr.z + v.secondDerivativeFactor * dr.z;

      third_derivative.m231 += v.thirdDerivativeFactor * dr.y * dr.z * dr.x;
      third_derivative.m232 += v.thirdDerivativeFactor * dr.y * dr.z * dr.y + v.secondDerivativeFactor * dr.z;
      third_derivative.m233 += v.thirdDerivativeFactor * dr.y * dr.z * dr.z + v.secondDerivativeFactor * dr.y;

      third_derivative.m311 += v.thirdDerivativeFactor * dr.z * dr.x * dr.x + v.secondDerivativeFactor * dr.z;
      third_derivative.m312 += v.thirdDerivativeFactor * dr.z * dr.x * dr.y;
      third_derivative.m313 += v.thirdDerivativeFactor * dr.z * dr.x * dr.z + v.secondDerivativeFactor * dr.x;

      third_derivative.m321 += v.thirdDerivativeFactor * dr.z * dr.y * dr.x;
      third_derivative.m322 += v.thirdDerivativeFactor * dr.z * dr.y * dr.y + v.secondDerivativeFactor * dr.z;
      third_derivative.m323 += v.thirdDerivativeFactor * dr.z * dr.y * dr.z + v.secondDerivativeFactor * dr.y;

      third_derivative.m331 += v.thirdDerivativeFactor * dr.z * dr.z * dr.x + v.secondDerivativeFactor * dr.x;
      third_derivative.m332 += v.thirdDerivativeFactor * dr.z * dr.z * dr.y + v.secondDerivativeFactor * dr.y;
      third_derivative.m333 += v.thirdDerivativeFactor * dr.z * dr.z * dr.z + 3.0 * v.secondDerivativeFactor * dr.z;
    }
  }

  return {energy, first_derivative, second_derivative, third_derivative};
}

