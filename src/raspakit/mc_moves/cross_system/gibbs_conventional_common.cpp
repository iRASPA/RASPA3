module;

module mc_moves_gibbs_conventional_common;

import std;

import randomnumbers;
import running_energy;
import system;
import molecule;
import atom;
import component;
import property_lambda_probability_histogram;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;
import cbmc;
import cbmc_chain_data;

namespace
{

enum class GibbsMoveKind : std::uint8_t
{
  LambdaChange = 0,
  Insert = 1,
  Delete = 2
};

std::size_t lambdaBinFromValue(const PropertyLambdaProbabilityHistogram& lambda, double lambdaValue)
{
  if (lambda.numberOfSamplePoints == 0)
  {
    return 0;
  }
  std::size_t bin = static_cast<std::size_t>(static_cast<double>(lambda.numberOfSamplePoints) * lambdaValue);
  if (bin == lambda.numberOfSamplePoints)
  {
    --bin;
  }
  return bin;
}

void applyTargetScaling(std::span<Atom> atoms, GibbsMoveKind moveKind, double lambda)
{
  for (Atom& atom : atoms)
  {
    switch (moveKind)
    {
      case GibbsMoveKind::Insert:
        atom.setScalingToInteger();
        break;
      case GibbsMoveKind::Delete:
        atom.setScalingToInactiveFractional();
        break;
      case GibbsMoveKind::LambdaChange:
        atom.setScaling(lambda);
        break;
    }
  }
}

std::optional<RunningEnergy> computeMoleculeEnergyDifference(System& system, std::span<const Atom> trialAtoms,
                                                             std::span<const Atom> oldAtoms)
{
  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), trialAtoms, oldAtoms);
  if (!frameworkDifference.has_value())
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialAtoms, oldAtoms);
  if (!moleculeDifference.has_value())
  {
    return std::nullopt;
  }

  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, trialAtoms, oldAtoms, system.netCharge);

  RunningEnergy tailDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                            system.spanOfMoleculeAtoms(), trialAtoms, oldAtoms) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), trialAtoms, oldAtoms);

  if (system.hasExternalField)
  {
    std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        trialAtoms, oldAtoms);
    if (!externalFieldDifference.has_value())
    {
      return std::nullopt;
    }
    return frameworkDifference.value() + moleculeDifference.value() + ewaldDifference + tailDifference +
           externalFieldDifference.value();
  }

  return frameworkDifference.value() + moleculeDifference.value() + ewaldDifference + tailDifference;
}

void acceptInsertStep(System& system, std::size_t selectedComponent, std::size_t fractionalMoleculeIndex,
                      const Molecule& molecule, std::vector<Atom> atoms, double remainderLambda)
{
  Component& component = system.components[selectedComponent];
  std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  for (Atom& atom : fractionalMolecule)
  {
    atom.setScalingToInteger();
  }

  for (Atom& atom : atoms)
  {
    atom.componentId = static_cast<std::uint8_t>(selectedComponent);
    atom.setScalingToFractional(remainderLambda, component.lambdaGibbs.dUdlambdaGroupId);
  }

  system.insertMolecule(selectedComponent, molecule, atoms);

  // Re-fetch spans: insertMolecule may reallocate atom storage.
  const std::size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
  std::span<Atom> fractionalMoleculeAfterInsert =
      system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
  std::swap_ranges(fractionalMoleculeAfterInsert.begin(), fractionalMoleculeAfterInsert.end(), lastMolecule.begin());
  std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, fractionalMoleculeIndex)],
            system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeId)]);

  component.lambdaGibbs.setCurrentBin(lambdaBinFromValue(component.lambdaGibbs, remainderLambda));
  system.updateMoleculeAtomInformation();
}

void acceptDeleteStep(System& system, std::size_t selectedComponent, std::size_t fractionalMoleculeIndex,
                      std::size_t selectedIntegerMoleculeIndex, double remainderLambda)
{
  Component& component = system.components[selectedComponent];
  std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  std::span<Atom> selectedMolecule = system.spanOfMolecule(selectedComponent, selectedIntegerMoleculeIndex);

  for (Atom& atom : fractionalMolecule)
  {
    atom.setScalingToInactiveFractional();
  }

  std::vector<Atom> savedFractional(fractionalMolecule.begin(), fractionalMolecule.end());
  std::swap_ranges(selectedMolecule.begin(), selectedMolecule.end(), fractionalMolecule.begin());
  std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedIntegerMoleculeIndex)],
            system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, fractionalMoleculeIndex)]);

  std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  const std::uint8_t groupId = component.lambdaGibbs.dUdlambdaGroupId;
  for (Atom& atom : newFractionalMolecule)
  {
    atom.setScalingToFractional(remainderLambda, groupId);
  }

  system.deleteMolecule(selectedComponent, selectedIntegerMoleculeIndex, savedFractional);

  component.lambdaGibbs.setCurrentBin(lambdaBinFromValue(component.lambdaGibbs, remainderLambda));
  system.updateMoleculeAtomInformation();
}

RunningEnergy growEwaldTailDifference(System& system, std::span<const Atom> growAtoms)
{
  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, growAtoms, {}, system.netCharge);
  RunningEnergy tailDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), growAtoms, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), growAtoms, {});
  return ewaldDifference + tailDifference;
}

RunningEnergy retraceEwaldTailDifference(System& system, std::span<const Atom> retraceAtoms)
{
  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, {}, retraceAtoms, system.netCharge);
  RunningEnergy tailDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                            system.spanOfMoleculeAtoms(), {}, retraceAtoms) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), {}, retraceAtoms);
  return ewaldDifference + tailDifference;
}

}  // namespace

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsConventionalCommon::gibbsConventionalMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent, bool useCBMC)
{
  const Move::Types move =
      useCBMC ? Move::Types::GibbsConventionalCBCFCMC : Move::Types::GibbsConventionalCFCMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];

  componentA.mc_moves_statistics.addTrial(move);

  PropertyLambdaProbabilityHistogram& lambdaA = componentA.lambdaGibbs;
  PropertyLambdaProbabilityHistogram& lambdaB = componentB.lambdaGibbs;

  const double lambdaOldA = lambdaA.lambdaValue();
  const double lambdaOldB = lambdaB.lambdaValue();
  const std::size_t oldBinA = lambdaBinFromValue(lambdaA, lambdaOldA);
  const std::size_t oldBinB = lambdaBinFromValue(lambdaB, lambdaOldB);
  const double biasOldA = lambdaA.biasFactor[oldBinA];
  const double biasOldB = lambdaB.biasFactor[oldBinB];

  const double maxChange = componentA.mc_moves_statistics.getMaxChange(move, 2);
  const double vNew = (2.0 * random.uniform() - 1.0) * maxChange;

  double lambdaNewA = lambdaOldA + vNew;
  double lambdaNewB = lambdaOldB - vNew;

  GibbsMoveKind moveKindA = GibbsMoveKind::LambdaChange;
  if (lambdaNewA > 1.0)
  {
    moveKindA = GibbsMoveKind::Insert;
    lambdaNewA -= 1.0;
  }
  else if (lambdaNewA < 0.0)
  {
    moveKindA = GibbsMoveKind::Delete;
    lambdaNewA += 1.0;
  }

  GibbsMoveKind moveKindB = GibbsMoveKind::LambdaChange;
  if (lambdaNewB > 1.0)
  {
    moveKindB = GibbsMoveKind::Insert;
    lambdaNewB -= 1.0;
  }
  else if (lambdaNewB < 0.0)
  {
    moveKindB = GibbsMoveKind::Delete;
    lambdaNewB += 1.0;
  }

  if (moveKindA == GibbsMoveKind::Delete &&
      systemA.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
  {
    return std::nullopt;
  }
  if (moveKindB == GibbsMoveKind::Delete &&
      systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
  {
    return std::nullopt;
  }

  const std::size_t newBinA = lambdaBinFromValue(lambdaA, lambdaNewA);
  const std::size_t newBinB = lambdaBinFromValue(lambdaB, lambdaNewB);
  const double biasNewA = lambdaA.biasFactor[newBinA];
  const double biasNewB = lambdaB.biasFactor[newBinB];

  const std::size_t indexFractionalA =
      systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, selectedComponent);
  const std::size_t indexFractionalB =
      systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, selectedComponent);
  std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalA);
  std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalB);

  const auto originalStoredEikA = systemA.storedEik;
  const auto originalTotalEikA = systemA.totalEik;
  const auto originalStoredEikB = systemB.storedEik;
  const auto originalTotalEikB = systemB.totalEik;
  auto restoreEwaldState = [&]()
  {
    systemA.storedEik = originalStoredEikA;
    systemA.totalEik = originalTotalEikA;
    systemB.storedEik = originalStoredEikB;
    systemB.totalEik = originalTotalEikB;
  };

  // First step: the fractional molecule of each box takes its target scaling
  std::vector<Atom> trialFractionalA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
  applyTargetScaling(trialFractionalA, moveKindA, lambdaNewA);
  std::optional<RunningEnergy> energyFirstStepA =
      computeMoleculeEnergyDifference(systemA, trialFractionalA, fractionalMoleculeA);
  if (!energyFirstStepA.has_value())
  {
    restoreEwaldState();
    return std::nullopt;
  }

  std::vector<Atom> trialFractionalB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
  applyTargetScaling(trialFractionalB, moveKindB, lambdaNewB);
  std::optional<RunningEnergy> energyFirstStepB =
      computeMoleculeEnergyDifference(systemB, trialFractionalB, fractionalMoleculeB);
  if (!energyFirstStepB.has_value())
  {
    restoreEwaldState();
    return std::nullopt;
  }

  std::optional<std::vector<Atom>> savedFractionalA{};
  std::optional<std::vector<Atom>> savedFractionalB{};
  if (moveKindA != GibbsMoveKind::LambdaChange)
  {
    savedFractionalA = std::vector<Atom>(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    applyTargetScaling(fractionalMoleculeA, moveKindA, lambdaNewA);
    // The second boundary operation starts from the structure factor produced by
    // the first fractional-scaling operation, so reciprocal-space cross terms are retained.
    systemA.storedEik = systemA.totalEik;
  }
  if (moveKindB != GibbsMoveKind::LambdaChange)
  {
    savedFractionalB = std::vector<Atom>(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    applyTargetScaling(fractionalMoleculeB, moveKindB, lambdaNewB);
    systemB.storedEik = systemB.totalEik;
  }

  auto restoreFractionals = [&]()
  {
    if (savedFractionalA.has_value())
    {
      std::copy(savedFractionalA->begin(), savedFractionalA->end(), fractionalMoleculeA.begin());
    }
    if (savedFractionalB.has_value())
    {
      std::copy(savedFractionalB->begin(), savedFractionalB->end(), fractionalMoleculeB.begin());
    }
    restoreEwaldState();
  };

  struct BoundaryTrial
  {
    RunningEnergy energy{};
    double acceptanceFactor{1.0};
    std::pair<Molecule, std::vector<Atom>> conventionalInsert{};
    std::optional<ChainGrowData> cbmcInsert{};
    std::size_t selectedInteger{0};
  };

  auto constructBoundaryTrial = [&](System& system, Component& component, GibbsMoveKind moveKind, double lambdaNew,
                                    BoundaryTrial& trial) -> bool
  {
    const double integerCount =
        static_cast<double>(system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
    const double volume = system.simulationBox.volume;
    const CBMC::GrowContext context{
        system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
        system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
        system.spanOfMoleculeAtoms(), system.beta, system.forceField.cutOffFrameworkVDW,
        system.forceField.cutOffMoleculeVDW, system.forceField.cutOffCoulomb};

    if (moveKind == GibbsMoveKind::LambdaChange)
    {
      return true;
    }

    if (moveKind == GibbsMoveKind::Insert)
    {
      if (useCBMC)
      {
        trial.cbmcInsert = CBMC::growMoleculeSwapInsertion(
            random, context, component, selectedComponent, component.growType, system.numberOfMolecules(), lambdaNew,
            component.lambdaGibbs.dUdlambdaGroupId, true);
        if (!trial.cbmcInsert.has_value() || system.insideBlockedPockets(component, trial.cbmcInsert->atoms))
        {
          return false;
        }

        RunningEnergy ewaldTail = growEwaldTailDifference(system, trial.cbmcInsert->atoms);
        trial.energy = trial.cbmcInsert->energies + ewaldTail;
        const double idealGas = component.idealGasRosenbluthWeight.value_or(1.0);
        trial.acceptanceFactor = (trial.cbmcInsert->RosenbluthWeight / idealGas) *
                                 std::exp(-system.beta * ewaldTail.potentialEnergy()) *
                                 volume / (integerCount + 1.0);
        return true;
      }

      trial.conventionalInsert =
          component.equilibratedMoleculeRandomInBox(random, selectedComponent, system.simulationBox);
      std::vector<Atom> trialAtoms = trial.conventionalInsert.second;
      const std::size_t upcomingMoleculeId = system.numberOfMolecules();
      for (Atom& atom : trialAtoms)
      {
        atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeId);
        atom.componentId = static_cast<std::uint8_t>(selectedComponent);
        atom.setScalingToFractional(lambdaNew, component.lambdaGibbs.dUdlambdaGroupId);
      }
      std::optional<RunningEnergy> insertionEnergy =
          computeMoleculeEnergyDifference(system, trialAtoms, std::span<const Atom>{});
      if (!insertionEnergy.has_value())
      {
        return false;
      }
      trial.energy = insertionEnergy.value();
      trial.acceptanceFactor =
          std::exp(-system.beta * trial.energy.potentialEnergy()) * volume / (integerCount + 1.0);
      return true;
    }

    trial.selectedInteger = system.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedMolecule = system.spanOfMolecule(selectedComponent, trial.selectedInteger);
    if (useCBMC)
    {
      std::vector<Atom> oldSelectedMolecule(selectedMolecule.begin(), selectedMolecule.end());
      ChainRetraceData retraceData =
          CBMC::retraceMoleculeSwapDeletion(random, context, component, component.growType, selectedMolecule);
      RunningEnergy ewaldTail = retraceEwaldTailDifference(system, selectedMolecule);
      trial.acceptanceFactor = (integerCount / volume) /
                               (retraceData.RosenbluthWeight *
                                std::exp(system.beta * ewaldTail.potentialEnergy()));
      std::copy(oldSelectedMolecule.begin(), oldSelectedMolecule.end(), selectedMolecule.begin());
    }

    std::vector<Atom> trialSelected(selectedMolecule.begin(), selectedMolecule.end());
    for (Atom& atom : trialSelected)
    {
      atom.setScalingToFractional(lambdaNew, component.lambdaGibbs.dUdlambdaGroupId);
    }
    std::optional<RunningEnergy> deletionEnergy =
        computeMoleculeEnergyDifference(system, trialSelected, selectedMolecule);
    if (!deletionEnergy.has_value())
    {
      return false;
    }
    trial.energy = deletionEnergy.value();
    if (!useCBMC)
    {
      trial.acceptanceFactor =
          std::exp(-system.beta * trial.energy.potentialEnergy()) * integerCount / volume;
    }
    return true;
  };

  BoundaryTrial boundaryA{};
  BoundaryTrial boundaryB{};
  if (!constructBoundaryTrial(systemA, componentA, moveKindA, lambdaNewA, boundaryA) ||
      !constructBoundaryTrial(systemB, componentB, moveKindB, lambdaNewB, boundaryB))
  {
    restoreFractionals();
    return std::nullopt;
  }

  componentA.mc_moves_statistics.addConstructed(move);

  const double acceptanceProbability =
      std::exp(-systemA.beta * energyFirstStepA->potentialEnergy() -
               systemB.beta * energyFirstStepB->potentialEnergy()) *
      boundaryA.acceptanceFactor * boundaryB.acceptanceFactor *
      std::exp((biasNewA - biasOldA) + (biasNewB - biasOldB));

  if (random.uniform() >= acceptanceProbability)
  {
    restoreFractionals();
    return std::nullopt;
  }

  componentA.mc_moves_statistics.addAccepted(move);
  Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
  Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

  auto acceptBox = [&](System& system, Component& component, GibbsMoveKind moveKind, double lambdaNew,
                       std::size_t indexFractional, BoundaryTrial& boundary)
  {
    if (moveKind == GibbsMoveKind::LambdaChange)
    {
      for (Atom& atom : system.spanOfMolecule(selectedComponent, indexFractional))
      {
        atom.setScalingToFractional(lambdaNew, component.lambdaGibbs.dUdlambdaGroupId);
      }
      component.lambdaGibbs.setCurrentBin(lambdaBinFromValue(component.lambdaGibbs, lambdaNew));
      return;
    }

    if (moveKind == GibbsMoveKind::Insert)
    {
      if (useCBMC)
      {
        acceptInsertStep(system, selectedComponent, indexFractional, boundary.cbmcInsert->molecule,
                         boundary.cbmcInsert->atoms, lambdaNew);
      }
      else
      {
        acceptInsertStep(system, selectedComponent, indexFractional, boundary.conventionalInsert.first,
                         boundary.conventionalInsert.second, lambdaNew);
      }
      return;
    }
    acceptDeleteStep(system, selectedComponent, indexFractional, boundary.selectedInteger, lambdaNew);
  };

  acceptBox(systemA, componentA, moveKindA, lambdaNewA, indexFractionalA, boundaryA);
  acceptBox(systemB, componentB, moveKindB, lambdaNewB, indexFractionalB, boundaryB);

  RunningEnergy energyDifferenceA = energyFirstStepA.value() + boundaryA.energy;
  RunningEnergy energyDifferenceB = energyFirstStepB.value() + boundaryB.energy;
  return std::make_pair(energyDifferenceA, energyDifferenceB);
}
