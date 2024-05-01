module;

#ifdef USE_LEGACY_HEADERS
#include <numbers>
#include <optional>
#include <iostream>
#include <algorithm>
#include <vector>
#include <span>
#include <cmath>
#include <thread>
#include <future>
#include <deque>
#include <semaphore>
#include <atomic>
#include <utility>
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
#endif

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
//import threading;
import energy_factor;
import force_factor;
import framework;
import component;


// Used to compute the total framework-molecule energy
//

void Interactions::computeFrameworkMoleculeEnergy(const ForceField &forceField, const SimulationBox &simulationBox, 
                                                  std::span<const Atom> frameworkAtoms, 
                                                  std::span<const Atom> moleculeAtoms, 
                                                  RunningEnergy &energyStatus) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtoms.empty()) return;

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

      if (rr < cutOffVDWSquared)
      {
        EnergyFactor energyFactor = 
          potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energyStatus.frameworkMoleculeVDW += energyFactor.energy;
        energyStatus.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = 
          potentialCoulombEnergy(forceField, groupIdA, groupIdB, scaleCoulombA, scaleCoulombB, r, chargeA, chargeB);

        energyStatus.frameworkMoleculeCharge += energyFactor.energy;
        energyStatus.dudlambdaCharge += energyFactor.dUdlambda;
      }
    }
  }
}

void Interactions::computeFrameworkMoleculeTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                                      std::span<const Atom> frameworkAtoms, 
                                                      std::span<const Atom> moleculeAtoms, 
                                                      RunningEnergy &energyStatus) noexcept
{
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
      energyStatus.tail += scalingVDWA * scalingVDWB * temp;
      energyStatus.dudlambdaVDW += (groupIdA ? scalingVDWB * temp : 0.0)
                                 + (groupIdB ? scalingVDWA * temp : 0.0);
    }
  }
}



// Used in Translation and Rotation
//

[[nodiscard]] std::optional<RunningEnergy> 
Interactions::computeFrameworkMoleculeEnergyDifference(const ForceField &forceField, 
                                                       const SimulationBox &simulationBox, 
                                                       std::span<const Atom> frameworkAtoms, 
                                                       std::span<const Atom> newatoms, 
                                                       std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool noCharges = forceField.noCharges;
  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

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
        EnergyFactor energyFactor = 
          potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW +=  energyFactor.energy;
        energySum.dudlambdaVDW += energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = 
          potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

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
        EnergyFactor energyFactor = 
          potentialVDWEnergy(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.dudlambdaVDW -= energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        EnergyFactor energyFactor = 
          potentialCoulombEnergy(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.dudlambdaCharge -= energyFactor.dUdlambda;
      }
    }
  }

  return std::optional{energySum};
}

[[nodiscard]] RunningEnergy 
Interactions::computeFrameworkMoleculeTailEnergyDifference(const ForceField &forceField,                                    
                                                           const SimulationBox &simulationBox,
                                                           std::span<const Atom> frameworkAtoms,
                                                           std::span<const Atom> newatoms, 
                                                           std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

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

std::pair<ForceFactor, ForceFactor> Interactions::computeFrameworkMoleculeGradient(const ForceField &forceField,                                    
                                                           const SimulationBox &simulationBox,                              
                                                           std::span<Atom> frameworkAtoms,
                                                           std::span<Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB;
  double rr;

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  ForceFactor energyVDW{ 0.0, 0.0, 0.0 };
  ForceFactor energyCoulomb{ 0.0, 0.0, 0.0 };

  if (moleculeAtoms.empty()) return {energyVDW, energyCoulomb};

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

      if (rr < cutOffVDWSquared)
      {
        ForceFactor forceFactor =
          potentialVDWGradient(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energyVDW += forceFactor;

        const double3 f = forceFactor.forceFactor * dr;

        it2->gradient -= f;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        ForceFactor forceFactor =
          potentialCoulombGradient(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energyCoulomb += forceFactor;

        const double3 f = forceFactor.forceFactor * dr;

        it2->gradient -= f;
      }
    }
  }
  return {energyVDW, energyCoulomb};
}


[[nodiscard]] std::pair<EnergyStatus, double3x3> 
Interactions::computeFrameworkMoleculeEnergyStrainDerivative(const ForceField &forceField,
                                                             const std::vector<Framework> &frameworkComponents,
                                                             const std::vector<Component> &components,
                                                             const SimulationBox &simulationBox,
                                                             std::span<Atom> frameworkAtoms,
                                                             std::span<Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  double3x3 strainDerivativeTensor;
  EnergyStatus energy(1, frameworkComponents.size(), components.size());

  bool noCharges = forceField.noCharges;
  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  if (moleculeAtoms.empty()) return std::make_pair(energy, strainDerivativeTensor);

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
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

      if (rr < cutOffVDWSquared)
      {
        ForceFactor forceFactor =
            potentialVDWGradient(forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energy.frameworkComponentEnergy(compA, compB).VanDerWaals += EnergyFactor(forceFactor.energy, 0.0);

        const double3 g = forceFactor.forceFactor * dr;

        it2->gradient -= g;

        strainDerivativeTensor.ax += g.x*dr.x;
        strainDerivativeTensor.bx += g.y*dr.x;
        strainDerivativeTensor.cx += g.z*dr.x;

        strainDerivativeTensor.ay += g.x*dr.y;
        strainDerivativeTensor.by += g.y*dr.y;
        strainDerivativeTensor.cy += g.z*dr.y;

        strainDerivativeTensor.az += g.x*dr.z;
        strainDerivativeTensor.bz += g.y*dr.z;
        strainDerivativeTensor.cz += g.z*dr.z;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        ForceFactor forceFactor =
            potentialCoulombGradient(forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r,
                                     chargeA, chargeB);

        energy.frameworkComponentEnergy(compA, compB).CoulombicReal += EnergyFactor(forceFactor.energy, 0.0);

        const double3 g = forceFactor.forceFactor * dr;

        it2->gradient -= g;

        strainDerivativeTensor.ax += g.x*dr.x;
        strainDerivativeTensor.bx += g.y*dr.x;
        strainDerivativeTensor.cx += g.z*dr.x;

        strainDerivativeTensor.ay += g.x*dr.y;
        strainDerivativeTensor.by += g.y*dr.y;
        strainDerivativeTensor.cy += g.z*dr.y;

        strainDerivativeTensor.az += g.x*dr.z;
        strainDerivativeTensor.bz += g.y*dr.z;
        strainDerivativeTensor.cz += g.z*dr.z;
      }
    }
  }

  return std::make_pair(energy, strainDerivativeTensor);
}
