module;

module interactions_intermolecular;

import std;

import energy_status;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;
import potential_correction_vdw;
import potential_electrostatics;
import interactions_pair_kernel;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_dudlambda;
import energy_status_inter;
import running_energy;
import component;
import units;
import threadpool;

// used in volume moves for computing the state at a new box and new, scaled atom positions
RunningEnergy Interactions::computeInterMolecularEnergy(const ForceField& forceField, const SimulationBox& box,
                                                        std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  forEachMoleculeMoleculePair<0>(
      forceField, box, moleculeAtoms,
      [&energySum](std::size_t, std::size_t, const Atom& atomA, const Atom& atomB,
                   const Potentials::PairDerivatives<0>& factors, const double3&)
      {
        energySum.moleculeMoleculeVDW += factors.energy;
        energySum.addDudlambdaVDW(atomA.groupId, atomB.groupId, atomA.scalingVDW, atomB.scalingVDW, factors.dUdlambda);
      },
      [&energySum](std::size_t, std::size_t, const Atom& atomA, const Atom& atomB,
                   const Potentials::PairDerivatives<0>& factors, const double3&)
      {
        energySum.moleculeMoleculeCharge += factors.energy;
        energySum.addDudlambdaCharge(atomA.groupId, atomB.groupId, atomA.scalingCoulomb, atomB.scalingCoulomb,
                                     factors.dUdlambda);
      });

  return energySum;
}

RunningEnergy Interactions::computeInterMolecularTailEnergy(const ForceField& forceField,
                                                            const SimulationBox& simulationBox,
                                                            std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;

    double temp_self = preFactor * forceField(typeA, typeA).tailCorrectionEnergy;
    energySum.tail += scalingVDWA * scalingVDWA * temp_self;
    energySum.addDudlambdaVDW(groupIdA, groupIdA, scalingVDWA, scalingVDWA, temp_self);

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
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

// used in mc_moves_translation.cpp, mc_moves_rotation.cpp,
//         mc_moves_random_translation.cpp, mc_moves_random_rotation.cpp
//         mc_moves_swap_cfcmc.cpp, mc_moves_swap_cfcmc_cbmc.cpp, mc_moves_gibbs_swap_cfcmc.cpp
[[nodiscard]] std::optional<RunningEnergy> Interactions::computeInterMolecularEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    double3 posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (const Atom& atom : newatoms)
    {
      std::size_t molB = static_cast<std::size_t>(atom.moleculeId);

      if (molA != molB)
      {
        double3 posB = atom.position;
        std::size_t typeB = static_cast<std::size_t>(atom.type);
        std::uint8_t groupIdB = atom.groupId;
        double scalingVDWB = atom.scalingVDW;
        double scalingCoulombB = atom.scalingCoulomb;
        double chargeB = atom.charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);
        }
      }
    }

    for (const Atom& atom : oldatoms)
    {
      std::size_t molB = static_cast<std::size_t>(atom.moleculeId);

      if (molA != molB)
      {
        double3 posB = atom.position;
        std::size_t typeB = static_cast<std::size_t>(atom.type);
        std::uint8_t groupIdB = atom.groupId;
        double scalingVDWB = atom.scalingVDW;
        double scalingCoulombB = atom.scalingCoulomb;
        double chargeB = atom.charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW -= energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge -= energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);
        }
      }
    }
  }

  return std::optional{energySum};
}

[[nodiscard]] RunningEnergy Interactions::computeInterMolecularTailEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;

    for (const Atom& atom : newatoms)
    {
      std::size_t molB = static_cast<std::size_t>(atom.moleculeId);

      if (molA != molB)
      {
        std::size_t typeB = static_cast<std::size_t>(atom.type);
        std::uint8_t groupIdB = atom.groupId;
        double scalingVDWB = atom.scalingVDW;

        double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
        energySum.tail += scalingVDWA * scalingVDWB * temp;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, temp);
      }
    }

    for (const Atom& atom : oldatoms)
    {
      std::size_t molB = static_cast<std::size_t>(atom.moleculeId);

      if (molA != molB)
      {
        std::size_t typeB = static_cast<std::size_t>(atom.type);
        std::uint8_t groupIdB = atom.groupId;
        double scalingVDWB = atom.scalingVDW;

        double temp = 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
        energySum.tail -= scalingVDWA * scalingVDWB * temp;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -temp);
      }
    }
  }

  for (const Atom& atomA : newatoms)
  {
    std::size_t typeA = static_cast<std::size_t>(atomA.type);
    std::uint8_t groupIdA = atomA.groupId;
    double scalingVDWA = atomA.scalingVDW;

    for (const Atom& atomB : newatoms)
    {
      std::size_t typeB = static_cast<std::size_t>(atomB.type);
      std::uint8_t groupIdB = atomB.groupId;
      double scalingVDWB = atomB.scalingVDW;

      double temp = preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail += scalingVDWA * scalingVDWB * temp;
      energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, temp);
    }
  }

  for (const Atom& atomA : oldatoms)
  {
    std::size_t typeA = static_cast<std::size_t>(atomA.type);
    std::uint8_t groupIdA = atomA.groupId;
    double scalingVDWA = atomA.scalingVDW;

    for (const Atom& atomB : oldatoms)
    {
      std::size_t typeB = static_cast<std::size_t>(atomB.type);
      std::uint8_t groupIdB = atomB.groupId;
      double scalingVDWB = atomB.scalingVDW;

      double temp = preFactor * forceField(typeA, typeB).tailCorrectionEnergy;
      energySum.tail -= scalingVDWA * scalingVDWB * temp;
      energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -temp);
    }
  }

  return energySum;
}

[[nodiscard]] RunningEnergy Interactions::computeInterMolecularTailEnergyFromTypeCounts(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const std::size_t> typeCounts) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  const std::size_t numberOfPseudoAtoms = forceField.numberOfPseudoAtoms;

  for (std::size_t typeA = 0; typeA < numberOfPseudoAtoms; ++typeA)
  {
    const double countA = static_cast<double>(typeCounts[typeA]);
    for (std::size_t typeB = 0; typeB < numberOfPseudoAtoms; ++typeB)
    {
      if (!forceField.tailCorrections[typeA * numberOfPseudoAtoms + typeB]) continue;

      const double countB = static_cast<double>(typeCounts[typeB]);
      energySum.tail += preFactor * countA * countB * forceField(typeA, typeB).tailCorrectionEnergy;
    }
  }

  return energySum;
}

[[nodiscard]] RunningEnergy Interactions::computeInterMolecularTailEnergyDifferenceAddRemove(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const std::size_t> currentTypeCounts,
    const Component& componentToAdd, const Component& componentToRemove) noexcept
{
  if (forceField.omitInterInteractions) return RunningEnergy{};

  std::vector<std::size_t> newTypeCounts(currentTypeCounts.begin(), currentTypeCounts.end());
  for (const Atom& atom : componentToAdd.atoms)
  {
    newTypeCounts[static_cast<std::size_t>(atom.type)] += 1;
  }
  for (const Atom& atom : componentToRemove.atoms)
  {
    newTypeCounts[static_cast<std::size_t>(atom.type)] -= 1;
  }

  return computeInterMolecularTailEnergyFromTypeCounts(forceField, simulationBox, newTypeCounts) -
         computeInterMolecularTailEnergyFromTypeCounts(forceField, simulationBox, currentTypeCounts);
}

[[nodiscard]] RunningEnergy Interactions::computeInterMolecularTailEnergyDifferenceReaction(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const std::size_t> currentTypeCounts,
    const std::vector<std::size_t>& reactantStoichiometry, const std::vector<std::size_t>& productStoichiometry,
    const std::vector<Component>& components, bool forward) noexcept
{
  if (forceField.omitInterInteractions) return RunningEnergy{};

  std::vector<std::size_t> newTypeCounts(currentTypeCounts.begin(), currentTypeCounts.end());

  if (forward)
  {
    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      for (std::size_t k = 0; k < reactantStoichiometry[componentId]; ++k)
      {
        for (const Atom& atom : components[componentId].atoms)
        {
          newTypeCounts[static_cast<std::size_t>(atom.type)] -= 1;
        }
      }
      for (std::size_t k = 0; k < productStoichiometry[componentId]; ++k)
      {
        for (const Atom& atom : components[componentId].atoms)
        {
          newTypeCounts[static_cast<std::size_t>(atom.type)] += 1;
        }
      }
    }
  }
  else
  {
    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      for (std::size_t k = 0; k < reactantStoichiometry[componentId]; ++k)
      {
        for (const Atom& atom : components[componentId].atoms)
        {
          newTypeCounts[static_cast<std::size_t>(atom.type)] += 1;
        }
      }
      for (std::size_t k = 0; k < productStoichiometry[componentId]; ++k)
      {
        for (const Atom& atom : components[componentId].atoms)
        {
          newTypeCounts[static_cast<std::size_t>(atom.type)] -= 1;
        }
      }
    }
  }

  return computeInterMolecularTailEnergyFromTypeCounts(forceField, simulationBox, newTypeCounts) -
         computeInterMolecularTailEnergyFromTypeCounts(forceField, simulationBox, currentTypeCounts);
}

RunningEnergy Interactions::computeInterMolecularGradient(const ForceField& forceField,
                                                          const SimulationBox& simulationBox,
                                                          std::span<const Atom> moleculeAtoms,
                                                          std::span<AtomDynamics> moleculeDynamics) noexcept
{
  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;

  forEachMoleculeMoleculePair<1>(
      forceField, simulationBox, moleculeAtoms,
      [&](std::size_t indexA, std::size_t indexB, const Atom& atomA, const Atom& atomB,
          const Potentials::PairDerivatives<1>& factors, const double3& dr)
      {
        energySum.moleculeMoleculeVDW += factors.energy;
        energySum.addDudlambdaVDW(atomA.groupId, atomB.groupId, atomA.scalingVDW, atomB.scalingVDW, factors.dUdlambda);

        const double3 f = factors.firstDerivativeFactor * dr;

        moleculeDynamics[indexA].gradient += f;
        moleculeDynamics[indexB].gradient -= f;
      },
      [&](std::size_t indexA, std::size_t indexB, const Atom& atomA, const Atom& atomB,
          const Potentials::PairDerivatives<1>& factors, const double3& dr)
      {
        energySum.moleculeMoleculeCharge += factors.energy;
        energySum.addDudlambdaCharge(atomA.groupId, atomB.groupId, atomA.scalingCoulomb, atomB.scalingCoulomb,
                                     factors.dUdlambda);

        const double3 f = factors.firstDerivativeFactor * dr;

        moleculeDynamics[indexA].gradient += f;
        moleculeDynamics[indexB].gradient -= f;
      });

  return energySum;
}

// Used in smart-MC
void Interactions::computeInterMolecularGradientMolecule(const ForceField& forceField,
                                                         const SimulationBox& simulationBox,
                                                         std::span<const Atom> moleculeAtoms,
                                                         std::span<const Atom> selectedAtoms,
                                                         std::span<AtomDynamics> selectedDynamics) noexcept
{
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  bool useCharge = forceField.useCharge;

  if (forceField.omitInterInteractions) return;
  if (selectedAtoms.empty() || moleculeAtoms.empty()) return;

  for (std::size_t indexA = 0; indexA != selectedAtoms.size(); ++indexA)
  {
    const Atom& atomA = selectedAtoms[indexA];

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      // Skip interactions with atoms of the selected molecule itself. Consistent with the energy-difference
      // kernel, molecules are distinguished by their (globally unique) molecule id; a trial configuration keeps
      // the same molecule id as the stored copy present in 'moleculeAtoms', so it does not interact with itself.
      if (atomA.moleculeId == it2->moleculeId) continue;

      evaluatePair<1>(
          forceField, simulationBox, atomA, *it2, cutOffMoleculeVDWSquared, cutOffChargeSquared, useCharge,
          [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
          { selectedDynamics[indexA].gradient += factors.firstDerivativeFactor * dr; },
          [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
          { selectedDynamics[indexA].gradient += factors.firstDerivativeFactor * dr; });
    }
  }
}

std::pair<EnergyStatus, double3x3> Interactions::computeInterMolecularEnergyStrainDerivative(
    const ForceField& forceField, const std::vector<Component>& components, const SimulationBox& simulationBox,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics) noexcept
{
  bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  EnergyStatus energy(1, 1, components.size());
  double3x3 strainDerivativeTensor{};

  if (forceField.omitInterInteractions) return {energy, strainDerivativeTensor};
  if (moleculeAtoms.empty()) return {energy, strainDerivativeTensor};

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(it1 - moleculeAtoms.begin());
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    std::size_t compA = static_cast<std::size_t>(it1->componentId);
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    double scalingVDWA = it1->scalingVDW;
    energy.componentEnergy(compA, compA).VanDerWaalsTailCorrection += EnergyDuDlambda(
        preFactor * scalingVDWA * scalingVDWA * forceField(typeA, typeA).tailCorrectionEnergy, 0.0);

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(it2 - moleculeAtoms.begin());
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);
      std::size_t compB = static_cast<std::size_t>(it2->componentId);
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      double scalingVDWB = it2->scalingVDW;

      EnergyDuDlambda temp(
          preFactor * scalingVDWA * scalingVDWB * forceField(typeA, typeB).tailCorrectionEnergy, 0.0);
      energy.componentEnergy(compA, compB).VanDerWaalsTailCorrection += 2.0 * temp;

      // skip interactions within the same molecule
      if (molA != molB)
      {
        const auto accumulateGradientAndStrain = [&](const double3& g, const double3& dr)
        {
          moleculeDynamics[indexA].gradient += g;
          moleculeDynamics[indexB].gradient -= g;
          accumulateStrainDerivative(strainDerivativeTensor, g, dr);
        };

        evaluatePair<1>(
            forceField, simulationBox, *it1, *it2, cutOffMoleculeVDWSquared, cutOffChargeSquared, useCharge,
            [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
            {
              energy.componentEnergy(compA, compB).VanDerWaals += 0.5 * EnergyDuDlambda(factors.energy, 0.0);
              energy.componentEnergy(compB, compA).VanDerWaals += 0.5 * EnergyDuDlambda(factors.energy, 0.0);
              accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
            },
            [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
            {
              energy.componentEnergy(compA, compB).CoulombicReal += 0.5 * EnergyDuDlambda(factors.energy, 0.0);
              energy.componentEnergy(compB, compA).CoulombicReal += 0.5 * EnergyDuDlambda(factors.energy, 0.0);
              accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
            });
      }
    }
  }

  return {energy, strainDerivativeTensor};
}

void Interactions::computeInterMolecularElectrostaticPotential(const ForceField& forceField, const SimulationBox& box,
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
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);

      // skip interactions within the same molecule
      if (molA != molB)
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

          std::size_t indexA = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it1));
          electricPotentialMolecules[indexA] += scalingCoulombB * chargeB * potential;

          std::size_t indexB = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it2));
          electricPotentialMolecules[indexB] += scalingCoulombA * chargeA * potential;
        }
      }
    }
  }
}

RunningEnergy Interactions::computeInterMolecularElectricField(const ForceField& forceField, const SimulationBox& box,
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
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);

      // skip interactions within the same molecule
      if (molA != molB)
      {
        posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        std::uint8_t groupIdB = it2->groupId;
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = box.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);
          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

          Potentials::PairDerivatives<1> gradient = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);

          double gradientFactorA = scalingCoulombB * chargeB * gradient.firstDerivativeFactor;
          std::size_t indexA = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it1));
          electricFieldMolecules[indexA] -= gradientFactorA * dr;

          double gradientFactorB = scalingCoulombA * chargeA * gradient.firstDerivativeFactor;
          std::size_t indexB = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it2));
          electricFieldMolecules[indexB] += gradientFactorB * dr;
        }
      }
    }
  }

  return energySum;
}

std::optional<RunningEnergy> Interactions::computeInterMolecularElectricFieldDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<double3> electricFieldMolecules,
    std::span<double3> electricFieldMolecule, std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;
  if (forceField.omitInterPolarization) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it1));
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    double3 posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(newatoms.begin(), it2));
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);

      if (molA != molB)
      {
        double3 posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        std::uint8_t groupIdB = it2->groupId;
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

          Potentials::PairDerivatives<1> gradient = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);

          double gradientFactorA = scalingCoulombB * chargeB * gradient.firstDerivativeFactor;
          electricFieldMolecules[indexA] -= gradientFactorA * dr;

          double gradientFactorB = scalingCoulombA * chargeA * gradient.firstDerivativeFactor;
          electricFieldMolecule[indexB] += gradientFactorB * dr;
        }
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(oldatoms.begin(), it2));
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);

      if (molA != molB)
      {
        double3 posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        std::uint8_t groupIdB = it2->groupId;
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW -= energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge -= energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);

          Potentials::PairDerivatives<1> gradient = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);

          double gradientFactorA = scalingCoulombB * chargeB * gradient.firstDerivativeFactor;
          electricFieldMolecules[indexA] += gradientFactorA * dr;

          double gradientFactorB = scalingCoulombA * chargeA * gradient.firstDerivativeFactor;
          electricFieldMolecule[indexB] -= gradientFactorB * dr;
        }
      }
    }
  }

  return std::optional{energySum};
}

std::optional<RunningEnergy> Interactions::computeInterMolecularPolarizationElectricFieldDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<double3> electricFieldNeighborDelta,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr;
  double rr;

  RunningEnergy energySum{};

  if (forceField.omitInterInteractions) return energySum;
  if (forceField.omitInterPolarization) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  // it1 runs over every atom currently in the system (the "neighbors"); atoms belonging to the moved molecule are
  // skipped through the moleculeId comparison so their neighbor-delta entries remain zero.
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it1));
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    double3 posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(newatoms.begin(), it2));
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);

      if (molA != molB)
      {
        double3 posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        std::uint8_t groupIdB = it2->groupId;
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.moleculeMoleculeVDW += energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge += energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

          Potentials::PairDerivatives<1> gradient = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);

          // field on the neighbor atom due to the moved atom at its new position (added contribution)
          double gradientFactorA = scalingCoulombB * chargeB * gradient.firstDerivativeFactor;
          electricFieldNeighborDelta[indexA] -= gradientFactorA * dr;

          // field on the moved atom (new position) due to the neighbor
          double gradientFactorB = scalingCoulombA * chargeA * gradient.firstDerivativeFactor;
          electricFieldMoleculeNew[indexB] += gradientFactorB * dr;
        }
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(oldatoms.begin(), it2));
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);

      if (molA != molB)
      {
        double3 posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        std::uint8_t groupIdB = it2->groupId;
        double scalingVDWB = it2->scalingVDW;
        double scalingCoulombB = it2->scalingCoulomb;
        double chargeB = it2->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffMoleculeVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.moleculeMoleculeVDW -= energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.moleculeMoleculeCharge -= energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);

          Potentials::PairDerivatives<1> gradient = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);

          // remove the contribution of the moved atom at its old position from the neighbor field
          double gradientFactorA = scalingCoulombB * chargeB * gradient.firstDerivativeFactor;
          electricFieldNeighborDelta[indexA] += gradientFactorA * dr;

          // field on the moved atom (old position) due to the neighbor. NOTE: the "Old" buffer follows the same
          // convention as computeFrameworkMoleculeEnergyDifference and energyDifferenceEwaldFourier, which store
          // the negated old field (it is only ever used squared in computePolarizationEnergyDifference); use -=
          // here as well so the framework, reciprocal and inter-molecular contributions stay sign-consistent.
          double gradientFactorB = scalingCoulombA * chargeA * gradient.firstDerivativeFactor;
          electricFieldMoleculeOld[indexB] -= gradientFactorB * dr;
        }
      }
    }
  }

  return std::optional{energySum};
}
