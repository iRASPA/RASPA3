module;

module interactions_framework_molecule;

import std;

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
import interactions_pair_kernel;
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
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;

  if (!framework.has_value()) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (const Atom& moleculeAtom : moleculeAtoms)
  {
    std::size_t typeB = static_cast<std::size_t>(moleculeAtom.type);
    bool isFractional = static_cast<bool>(moleculeAtom.isFractional);

    if (interpolationGrids[typeB].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      energySum.frameworkMoleculeVDW += interpolationGrids[typeB]->interpolate(moleculeAtom.position);
      if (useCharge)
      {
        energySum.frameworkMoleculeCharge +=
            moleculeAtom.charge * interpolationGrids.back()->interpolate(moleculeAtom.position);
      }
    }
    else
    {
      forEachFrameworkMoleculePair<0>(
          forceField, simulationBox, moleculeAtom, frameworkAtoms,
          [&](std::size_t, const Atom& frameworkAtom, const Potentials::PairDerivatives<0>& factors, const double3&)
          {
            energySum.frameworkMoleculeVDW += factors.energy;
            energySum.addDudlambdaVDW(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingVDW,
                                      frameworkAtom.scalingVDW, factors.dUdlambda);
          },
          [&](std::size_t, const Atom& frameworkAtom, const Potentials::PairDerivatives<0>& factors, const double3&)
          {
            energySum.frameworkMoleculeCharge += factors.energy;
            energySum.addDudlambdaCharge(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingCoulomb,
                                         frameworkAtom.scalingCoulomb, factors.dUdlambda);
          });
    }
  }
  return energySum;
}

RunningEnergy Interactions::computeFrameworkMoleculeTailEnergy(const ForceField& forceField,
                                                               const SimulationBox& simulationBox,
                                                               std::span<const Atom> frameworkAtoms,
                                                               std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, temp);
    }
  }

  return energySum;
}

// Used in Translation and Rotation
//

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (!framework.has_value()) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (auto& atom : newatoms)
  {
    double3 posB = atom.position;
    std::size_t typeB = static_cast<std::size_t>(atom.type);
    std::uint8_t groupIdB = atom.groupId;
    bool isFractional = static_cast<bool>(atom.isFractional);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
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
        std::size_t typeA = static_cast<std::size_t>(it1->type);
        std::uint8_t groupIdA = it1->groupId;
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor =
              Potentials::potentialVDWEnergy(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor =
              Potentials::potentialCoulombEnergy(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);
        }
      }
    }
  }

  for (auto& atom : oldatoms)
  {
    double3 posB = atom.position;
    std::size_t typeB = static_cast<std::size_t>(atom.type);
    std::uint8_t groupIdB = atom.groupId;
    bool isFractional = static_cast<bool>(atom.isFractional);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
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
        std::size_t typeA = static_cast<std::size_t>(it1->type);
        std::uint8_t groupIdA = it1->groupId;
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::EnergyFactor energyFactor =
              Potentials::potentialVDWEnergy(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW -= energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::EnergyFactor energyFactor =
              Potentials::potentialCoulombEnergy(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge -= energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);
        }
      }
    }
  }

  return std::optional{energySum};
}

std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox,
    [[maybe_unused]] const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    [[maybe_unused]] const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(newatoms.begin(), it2));
      double3 posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialVDWEnergy(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialCoulombEnergy(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA * Potentials::potentialCoulombGradient(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeNew[indexB] += gradientFactor.gradientFactor * dr;
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(oldatoms.begin(), it2));
      double3 posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialVDWEnergy(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialCoulombEnergy(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA * Potentials::potentialCoulombGradient(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeOld[indexB] -= gradientFactor.gradientFactor * dr;
      }
    }
  }

  return std::optional{energySum};
}

void Interactions::computeFrameworkMoleculeElectricFieldDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(newatoms.begin(), it2));
      double3 posB = it2->position;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA * Potentials::potentialCoulombGradient(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeNew[indexB] += gradientFactor.gradientFactor * dr;
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(oldatoms.begin(), it2));
      double3 posB = it2->position;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA * Potentials::potentialCoulombGradient(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeOld[indexB] -= gradientFactor.gradientFactor * dr;
      }
    }
  }
}

[[nodiscard]] RunningEnergy Interactions::computeFrameworkMoleculeTailEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;

    for (const Atom& atom : newatoms)
    {
      std::size_t typeB = static_cast<std::size_t>(atom.type);
      std::uint8_t groupIdB = atom.groupId;
      double scalingVDWB = atom.scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, temp);
    }

    for (const Atom& atom : oldatoms)
    {
      std::size_t typeB = static_cast<std::size_t>(atom.type);
      std::uint8_t groupIdB = atom.groupId;
      double scalingVDWB = atom.scalingVDW;

      double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail -= scalingVDWA * scalingVDWB * temp;
      energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -temp);
    }
  }

  return energySum;
}

RunningEnergy Interactions::computeFrameworkMoleculeGradient(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<AtomDynamics> frameworkDynamics) noexcept
{
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const bool flexibleFramework = framework && !framework->rigid && frameworkDynamics.size() == frameworkAtoms.size();

  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(it1 - moleculeAtoms.begin());
    const Atom& moleculeAtom = *it1;
    std::size_t typeA = static_cast<std::size_t>(moleculeAtom.type);
    bool isFractional = static_cast<bool>(moleculeAtom.isFractional);

    if (!flexibleFramework && interpolationGrids[typeA].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      auto [energy_vdw, gradient_vdw] = interpolationGrids[typeA]->interpolateGradient(moleculeAtom.position);
      energySum.frameworkMoleculeVDW += energy_vdw;
      moleculeDynamics[indexA].gradient += gradient_vdw;
      if (useCharge)
      {
        auto [energy_real_ewald, gradient_real_ewald] =
            interpolationGrids.back()->interpolateGradient(moleculeAtom.position);
        energySum.frameworkMoleculeCharge += moleculeAtom.charge * energy_real_ewald;
        moleculeDynamics[indexA].gradient += moleculeAtom.charge * gradient_real_ewald;
      }
    }
    else
    {
      forEachFrameworkMoleculePair<1>(
          forceField, simulationBox, moleculeAtom, frameworkAtoms,
          [&](std::size_t indexB, const Atom& frameworkAtom, const Potentials::PairDerivatives<1>& factors,
              const double3& dr)
          {
            energySum.frameworkMoleculeVDW += factors.energy;
            energySum.addDudlambdaVDW(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingVDW,
                                      frameworkAtom.scalingVDW, factors.dUdlambda);

            const double3 f = factors.firstDerivativeFactor * dr;

            moleculeDynamics[indexA].gradient += f;
            if (flexibleFramework) frameworkDynamics[indexB].gradient -= f;
          },
          [&](std::size_t indexB, const Atom& frameworkAtom, const Potentials::PairDerivatives<1>& factors,
              const double3& dr)
          {
            energySum.frameworkMoleculeCharge += factors.energy;
            energySum.addDudlambdaCharge(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingCoulomb,
                                         frameworkAtom.scalingCoulomb, factors.dUdlambda);

            const double3 g = factors.firstDerivativeFactor * dr;

            moleculeDynamics[indexA].gradient += g;
            if (flexibleFramework) frameworkDynamics[indexB].gradient -= g;
          });
    }
  }
  return energySum;
}

[[nodiscard]] std::pair<EnergyStatus, double3x3> Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
    const ForceField& forceField, const std::optional<Framework>& framework,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<Component>& components, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics) noexcept
{
  double3x3 strainDerivativeTensor;
  EnergyStatus energy(1, framework.has_value() ? 1 : 0, components.size());

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  if (moleculeAtoms.empty()) return std::make_pair(energy, strainDerivativeTensor);

  // Framework atoms are fixed, so their gradient is never consumed and is not accumulated.
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(it1 - moleculeAtoms.begin());
    const double3 posA = it1->position;
    std::size_t compA = static_cast<std::size_t>(it1->componentId);
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    bool isFractional = static_cast<bool>(it1->isFractional);
    double scalingVDWA = it1->scalingVDW;
    double chargeA = it1->charge;

    if (interpolationGrids[typeA].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      auto [energy_vdw, gradient_vdw] = interpolationGrids[typeA]->interpolateGradient(posA);
      energy.frameworkComponentEnergy(0, compA).VanDerWaals += Potentials::EnergyFactor(energy_vdw, 0.0);
      const double3 f = gradient_vdw;

      moleculeDynamics[indexA].gradient += f;
      accumulateStrainDerivative(strainDerivativeTensor, f, posA);

      if (useCharge)
      {
        auto [energy_real_ewald, gradient_real_ewald] = interpolationGrids.back()->interpolateGradient(posA);
        energy.frameworkComponentEnergy(0, compA).CoulombicReal +=
            Potentials::EnergyFactor(chargeA * energy_real_ewald, 0.0);
        const double3 g = chargeA * gradient_real_ewald;

        moleculeDynamics[indexA].gradient += g;
        accumulateStrainDerivative(strainDerivativeTensor, g, posA);
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it2 = frameworkAtoms.begin(); it2 != frameworkAtoms.end(); ++it2)
      {
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        double scalingVDWB = it2->scalingVDW;

        Potentials::EnergyFactor temp(
            preFactor * scalingVDWA * scalingVDWB * forceField(typeB, typeA).tailCorrectionEnergy, 0.0);
        energy.frameworkComponentEnergy(0, compA).VanDerWaalsTailCorrection += 2.0 * temp;

        const auto accumulateGradientAndStrain = [&](const double3& g, const double3& dr)
        {
          moleculeDynamics[indexA].gradient += g;
          accumulateStrainDerivative(strainDerivativeTensor, g, dr);
        };

        evaluatePair<1>(
            forceField, simulationBox, *it1, *it2, cutOffFrameworkVDWSquared, cutOffChargeSquared, useCharge,
            [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
            {
              energy.frameworkComponentEnergy(0, compA).VanDerWaals += Potentials::EnergyFactor(factors.energy, 0.0);
              accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
            },
            [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
            {
              energy.frameworkComponentEnergy(0, compA).CoulombicReal += Potentials::EnergyFactor(factors.energy, 0.0);
              accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
            });
      }
    }
  }

  return std::make_pair(energy, strainDerivativeTensor);
}

void Interactions::computeFrameworkMoleculeElectrostaticPotential(const ForceField& forceField,
                                                                  const SimulationBox& simulationBox,
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

        double potential = Potentials::potentialElectrostatics(forceField, 1.0, r, 1.0);

        std::size_t indexB = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it2));
        electricPotentialMolecules[indexB] += scalingCoulombA * chargeA * potential;
      }
    }
  }
}

RunningEnergy Interactions::computeFrameworkMoleculeElectricField(const ForceField& forceField,
                                                                  const SimulationBox& simulationBox,
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
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialVDWEnergy(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::EnergyFactor energyFactor =
            Potentials::potentialCoulombEnergy(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

        Potentials::GradientFactor gradientFactor =
            scalingCoulombA * chargeA * Potentials::potentialCoulombGradient(forceField, 1.0, 1.0, r, 1.0, 1.0);
        std::size_t index = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it2));
        electricFieldMolecules[index] += gradientFactor.gradientFactor * dr;
      }
    }
  }

  return energySum;
}

std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionVDW(const ForceField& forceField,
                                                                                   const SimulationBox& simulationBox,
                                                                                   double3 posA, std::size_t typeA,
                                                                                   std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    std::size_t typeB = static_cast<std::size_t>(it1->type);
    double scalingB = it1->scalingVDW;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < cutOffFrameworkVDWSquared)
    {
      Potentials::HessianFactor v = Potentials::potentialVDWHessian(forceField, 1.0, scalingB, rr, typeA, typeB);

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
    const ForceField& forceField, const SimulationBox& simulationBox, double3 posA, double chargeA,
    std::span<const Atom> frameworkAtoms)
{
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    double scalingB = it1->scalingVDW;
    double chargeB = it1->charge;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < cutOffChargeSquared)
    {
      double r = std::sqrt(rr);
      Potentials::HessianFactor v =
          Potentials::potentialCoulombHessian(forceField, 1.0, scalingB, rr, r, chargeA, chargeB);

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
