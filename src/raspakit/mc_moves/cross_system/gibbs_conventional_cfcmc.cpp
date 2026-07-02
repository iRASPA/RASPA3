module;

module mc_moves_gibbs_conventional_cfcmc;

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
  std::size_t bin = static_cast<std::size_t>(lambda.numberOfSamplePoints * lambdaValue);
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
        atom.setScalingFullyOff();
        atom.groupId = false;
        atom.isFractional = true;
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

std::optional<RunningEnergy> computeFirstStepEnergy(System& system, std::size_t selectedComponent,
                                                      std::size_t fractionalMoleculeIndex, GibbsMoveKind moveKind,
                                                      double lambdaNew)
{
  std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  std::vector<Atom> trialAtoms(fractionalMolecule.begin(), fractionalMolecule.end());
  applyTargetScaling(trialAtoms, moveKind, lambdaNew);
  return computeMoleculeEnergyDifference(system, trialAtoms, fractionalMolecule);
}

std::vector<Atom> applyFirstStepToFractional(System& system, std::size_t selectedComponent,
                                             std::size_t fractionalMoleculeIndex, GibbsMoveKind moveKind,
                                             double lambdaNew)
{
  std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  std::vector<Atom> savedFractional(fractionalMolecule.begin(), fractionalMolecule.end());
  applyTargetScaling(fractionalMolecule, moveKind, lambdaNew);
  return savedFractional;
}

void restoreFractional(System& system, std::size_t selectedComponent, std::size_t fractionalMoleculeIndex,
                       const std::vector<Atom>& savedFractional)
{
  std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  std::copy(savedFractional.begin(), savedFractional.end(), fractionalMolecule.begin());
}

void acceptLambdaScaling(System& system, std::size_t selectedComponent, std::size_t fractionalMoleculeIndex,
                         double lambdaNew)
{
  std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  const bool groupId = system.components[selectedComponent].lambdaGibbs.computeDUdlambda;
  for (Atom& atom : fractionalMolecule)
  {
    atom.setScaling(lambdaNew);
    atom.isFractional = true;
    atom.groupId = groupId;
  }
  Component& component = system.components[selectedComponent];
  component.lambdaGibbs.setCurrentBin(lambdaBinFromValue(component.lambdaGibbs, lambdaNew));
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
    atom.groupId = component.lambdaGibbs.computeDUdlambda;
    atom.setScaling(remainderLambda);
    atom.isFractional = true;
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
    atom.setScalingFullyOff();
    atom.groupId = false;
    atom.isFractional = true;
  }

  std::vector<Atom> savedFractional(fractionalMolecule.begin(), fractionalMolecule.end());
  std::swap_ranges(selectedMolecule.begin(), selectedMolecule.end(), fractionalMolecule.begin());
  std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedIntegerMoleculeIndex)],
            system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, fractionalMoleculeIndex)]);

  std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, fractionalMoleculeIndex);
  const bool groupId = component.lambdaGibbs.computeDUdlambda;
  for (Atom& atom : newFractionalMolecule)
  {
    atom.setScaling(remainderLambda);
    atom.isFractional = true;
    atom.groupId = groupId;
  }

  system.deleteMolecule(selectedComponent, selectedIntegerMoleculeIndex, savedFractional);

  component.lambdaGibbs.setCurrentBin(lambdaBinFromValue(component.lambdaGibbs, remainderLambda));
  system.updateMoleculeAtomInformation();
}

void acceptInsertStep(System& system, std::size_t selectedComponent, std::size_t fractionalMoleculeIndex,
                      const ChainGrowData& growData, double remainderLambda)
{
  acceptInsertStep(system, selectedComponent, fractionalMoleculeIndex, growData.molecule, growData.atom,
                   remainderLambda);
}

RunningEnergy growEwaldTailDifference(System& system, std::span<const Atom> growAtoms)
{
  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
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

std::optional<std::pair<RunningEnergy, RunningEnergy>> gibbsConventionalCFCMCMoveImpl(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent, bool useCBMC)
{
  const Move::Types move =
      useCBMC ? Move::Types::GibbsConventionalCFCMCCBMC : Move::Types::GibbsConventionalCFCMC;
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

  std::optional<RunningEnergy> energyFirstStepA =
      computeFirstStepEnergy(systemA, selectedComponent, indexFractionalA, moveKindA,
                             moveKindA == GibbsMoveKind::LambdaChange ? lambdaNewA : 0.0);
  if (!energyFirstStepA.has_value())
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> energyFirstStepB =
      computeFirstStepEnergy(systemB, selectedComponent, indexFractionalB, moveKindB,
                             moveKindB == GibbsMoveKind::LambdaChange ? lambdaNewB : 0.0);
  if (!energyFirstStepB.has_value())
  {
    return std::nullopt;
  }

  std::optional<std::vector<Atom>> savedFractionalA{};
  std::optional<std::vector<Atom>> savedFractionalB{};
  if (moveKindA != GibbsMoveKind::LambdaChange)
  {
    savedFractionalA =
        applyFirstStepToFractional(systemA, selectedComponent, indexFractionalA, moveKindA, lambdaNewA);
  }
  if (moveKindB != GibbsMoveKind::LambdaChange)
  {
    savedFractionalB =
        applyFirstStepToFractional(systemB, selectedComponent, indexFractionalB, moveKindB, lambdaNewB);
  }

  auto restoreFractionals = [&]()
  {
    if (savedFractionalA.has_value())
    {
      restoreFractional(systemA, selectedComponent, indexFractionalA, savedFractionalA.value());
    }
    if (savedFractionalB.has_value())
    {
      restoreFractional(systemB, selectedComponent, indexFractionalB, savedFractionalB.value());
    }
  };

  RunningEnergy energySecondStepA{};
  RunningEnergy energySecondStepB{};
  double rosenbluthNew = 1.0;
  double rosenbluthOld = 1.0;
  double rosenbluthNewA = 1.0;
  double rosenbluthOldA = 1.0;
  double rosenbluthNewB = 1.0;
  double rosenbluthOldB = 1.0;

  std::pair<Molecule, std::vector<Atom>> trialInsertA{};
  std::pair<Molecule, std::vector<Atom>> trialInsertB{};
  std::optional<ChainGrowData> growDataA{};
  std::optional<ChainGrowData> growDataB{};
  std::size_t selectedIntegerA = 0;
  std::size_t selectedIntegerB = 0;

  const double cutOffFrameworkVDWA = systemA.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDWA = systemA.forceField.cutOffMoleculeVDW;
  const double cutOffCoulombA = systemA.forceField.cutOffCoulomb;
  const double cutOffFrameworkVDWB = systemB.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDWB = systemB.forceField.cutOffMoleculeVDW;
  const double cutOffCoulombB = systemB.forceField.cutOffCoulomb;
  const Component::GrowType growTypeA = componentA.growType;
  const Component::GrowType growTypeB = componentB.growType;

  if (moveKindA == GibbsMoveKind::Insert)
  {
    if (useCBMC)
    {
      growDataA = CBMC::growMoleculeSwapInsertion(
          random, componentA, selectedComponent, systemA.hasExternalField, systemA.forceField, systemA.simulationBox,
          systemA.interpolationGrids, systemA.externalFieldInterpolationGrid, systemA.framework,
          systemA.spanOfFrameworkAtoms(), systemA.spanOfMoleculeAtoms(), systemA.beta, growTypeA, cutOffFrameworkVDWA,
          cutOffMoleculeVDWA, cutOffCoulombA, systemA.numberOfMolecules(), lambdaNewA,
          componentA.lambdaGibbs.computeDUdlambda, true);
      if (!growDataA.has_value() || systemA.insideBlockedPockets(componentA, growDataA->atom))
      {
        restoreFractionals();
        return std::nullopt;
      }

      RunningEnergy growEwaldTailA = growEwaldTailDifference(systemA, growDataA->atom);
      energySecondStepA += growDataA->energies + growEwaldTailA;

      const double idealGasA = componentA.idealGasRosenbluthWeight.value_or(1.0);
      const double correctionGrowA = std::exp(-systemA.beta * growEwaldTailA.potentialEnergy());
      rosenbluthNewA =
          (growDataA->RosenbluthWeight / idealGasA) * correctionGrowA *
          (static_cast<double>(systemB.numberOfMoleculesPerComponent[selectedComponent] - 1) *
           systemA.simulationBox.volume) /
          (static_cast<double>(systemA.numberOfMoleculesPerComponent[selectedComponent]) * systemB.simulationBox.volume);
    }
    else
    {
      trialInsertA = componentA.equilibratedMoleculeRandomInBox(random, selectedComponent, systemA.simulationBox);

      std::vector<Atom> trialAtomsA = trialInsertA.second;
      const bool groupIdA = componentA.lambdaGibbs.computeDUdlambda;
      const std::size_t upcomingMoleculeIdA = systemA.numberOfMolecules();
      for (Atom& atom : trialAtomsA)
      {
        atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeIdA);
        atom.componentId = static_cast<std::uint8_t>(selectedComponent);
        atom.groupId = groupIdA;
        atom.isFractional = true;
        atom.setScaling(lambdaNewA);
      }

      std::optional<RunningEnergy> insertEnergyA =
          computeMoleculeEnergyDifference(systemA, trialAtomsA, std::span<const Atom>{});
      if (!insertEnergyA.has_value())
      {
        restoreFractionals();
        return std::nullopt;
      }
      energySecondStepA += insertEnergyA.value();
    }

    selectedIntegerB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedMoleculeB = systemB.spanOfMolecule(selectedComponent, selectedIntegerB);
    if (useCBMC)
    {
      std::vector<Atom> oldSelectedMoleculeB(selectedMoleculeB.begin(), selectedMoleculeB.end());
      ChainRetraceData retraceDataB = CBMC::retraceMoleculeSwapDeletion(
          random, componentB, systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
          systemB.interpolationGrids, systemB.externalFieldInterpolationGrid, systemB.framework,
          systemB.spanOfFrameworkAtoms(), systemB.spanOfMoleculeAtoms(), systemB.beta, growTypeB, cutOffFrameworkVDWB,
          cutOffMoleculeVDWB, cutOffCoulombB, selectedMoleculeB);
      RunningEnergy retraceEwaldTailB = retraceEwaldTailDifference(systemB, selectedMoleculeB);
      const double correctionRetraceB = std::exp(systemB.beta * retraceEwaldTailB.potentialEnergy());
      rosenbluthOldB = retraceDataB.RosenbluthWeight * correctionRetraceB;
      std::copy(oldSelectedMoleculeB.begin(), oldSelectedMoleculeB.end(), selectedMoleculeB.begin());

      std::vector<Atom> trialSelectedB(selectedMoleculeB.begin(), selectedMoleculeB.end());
      applyTargetScaling(trialSelectedB, GibbsMoveKind::LambdaChange, lambdaNewB);
      for (Atom& atom : trialSelectedB)
      {
        atom.groupId = componentB.lambdaGibbs.computeDUdlambda;
        atom.isFractional = true;
      }
      std::optional<RunningEnergy> promoteEnergyB =
          computeMoleculeEnergyDifference(systemB, trialSelectedB, selectedMoleculeB);
      if (!promoteEnergyB.has_value())
      {
        restoreFractionals();
        return std::nullopt;
      }
      energySecondStepB += promoteEnergyB.value();
    }
    else
    {
      std::vector<Atom> trialSelectedB(selectedMoleculeB.begin(), selectedMoleculeB.end());
      applyTargetScaling(trialSelectedB, GibbsMoveKind::LambdaChange, lambdaNewB);
      for (Atom& atom : trialSelectedB)
      {
        atom.groupId = componentB.lambdaGibbs.computeDUdlambda;
        atom.isFractional = true;
      }

      std::optional<RunningEnergy> retraceEnergyB =
          computeMoleculeEnergyDifference(systemB, trialSelectedB, selectedMoleculeB);
      if (!retraceEnergyB.has_value())
      {
        restoreFractionals();
        return std::nullopt;
      }
      energySecondStepB += retraceEnergyB.value();

      const double totalSecondStep = (energySecondStepA + energySecondStepB).potentialEnergy();
      rosenbluthNew =
          std::exp(-systemA.beta * totalSecondStep) *
          (static_cast<double>(systemB.numberOfMoleculesPerComponent[selectedComponent] - 1) *
           systemA.simulationBox.volume) /
          (static_cast<double>(systemA.numberOfMoleculesPerComponent[selectedComponent]) *
           systemB.simulationBox.volume);
    }
  }
  else if (moveKindA == GibbsMoveKind::Delete)
  {
    selectedIntegerA = systemA.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedMoleculeA = systemA.spanOfMolecule(selectedComponent, selectedIntegerA);
    if (useCBMC)
    {
      std::vector<Atom> oldSelectedMoleculeA(selectedMoleculeA.begin(), selectedMoleculeA.end());
      ChainRetraceData retraceDataA = CBMC::retraceMoleculeSwapDeletion(
          random, componentA, systemA.hasExternalField, systemA.forceField, systemA.simulationBox,
          systemA.interpolationGrids, systemA.externalFieldInterpolationGrid, systemA.framework,
          systemA.spanOfFrameworkAtoms(), systemA.spanOfMoleculeAtoms(), systemA.beta, growTypeA, cutOffFrameworkVDWA,
          cutOffMoleculeVDWA, cutOffCoulombA, selectedMoleculeA);
      RunningEnergy retraceEwaldTailA = retraceEwaldTailDifference(systemA, selectedMoleculeA);
      const double correctionRetraceA = std::exp(systemA.beta * retraceEwaldTailA.potentialEnergy());
      rosenbluthOldA = retraceDataA.RosenbluthWeight * correctionRetraceA;
      std::copy(oldSelectedMoleculeA.begin(), oldSelectedMoleculeA.end(), selectedMoleculeA.begin());

      std::vector<Atom> trialSelectedA(selectedMoleculeA.begin(), selectedMoleculeA.end());
      applyTargetScaling(trialSelectedA, GibbsMoveKind::LambdaChange, lambdaNewA);
      for (Atom& atom : trialSelectedA)
      {
        atom.groupId = componentA.lambdaGibbs.computeDUdlambda;
        atom.isFractional = true;
      }
      std::optional<RunningEnergy> promoteEnergyA =
          computeMoleculeEnergyDifference(systemA, trialSelectedA, selectedMoleculeA);
      if (!promoteEnergyA.has_value())
      {
        restoreFractionals();
        return std::nullopt;
      }
      energySecondStepA += promoteEnergyA.value();
    }
    else
    {
      std::vector<Atom> trialSelectedA(selectedMoleculeA.begin(), selectedMoleculeA.end());
      applyTargetScaling(trialSelectedA, GibbsMoveKind::LambdaChange, lambdaNewA);
      for (Atom& atom : trialSelectedA)
      {
        atom.groupId = componentA.lambdaGibbs.computeDUdlambda;
        atom.isFractional = true;
      }

      std::optional<RunningEnergy> retraceEnergyA =
          computeMoleculeEnergyDifference(systemA, trialSelectedA, selectedMoleculeA);
      if (!retraceEnergyA.has_value())
      {
        restoreFractionals();
        return std::nullopt;
      }
      energySecondStepA += retraceEnergyA.value();
    }

    if (useCBMC)
    {
      growDataB = CBMC::growMoleculeSwapInsertion(
          random, componentB, selectedComponent, systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
          systemB.interpolationGrids, systemB.externalFieldInterpolationGrid, systemB.framework,
          systemB.spanOfFrameworkAtoms(), systemB.spanOfMoleculeAtoms(), systemB.beta, growTypeB, cutOffFrameworkVDWB,
          cutOffMoleculeVDWB, cutOffCoulombB, systemB.numberOfMolecules(), lambdaNewB,
          componentB.lambdaGibbs.computeDUdlambda, true);
      if (!growDataB.has_value() || systemB.insideBlockedPockets(componentB, growDataB->atom))
      {
        restoreFractionals();
        return std::nullopt;
      }

      RunningEnergy growEwaldTailB = growEwaldTailDifference(systemB, growDataB->atom);
      energySecondStepB += growDataB->energies + growEwaldTailB;

      const double idealGasB = componentB.idealGasRosenbluthWeight.value_or(1.0);
      const double correctionGrowB = std::exp(-systemB.beta * growEwaldTailB.potentialEnergy());
      rosenbluthNewB =
          (growDataB->RosenbluthWeight / idealGasB) * correctionGrowB *
          (static_cast<double>(systemA.numberOfMoleculesPerComponent[selectedComponent] - 1) *
           systemB.simulationBox.volume) /
          (static_cast<double>(systemB.numberOfMoleculesPerComponent[selectedComponent]) *
           systemA.simulationBox.volume);
    }
    else
    {
      trialInsertB = componentB.equilibratedMoleculeRandomInBox(random, selectedComponent, systemB.simulationBox);

      std::vector<Atom> trialAtomsB = trialInsertB.second;
      const bool groupIdB = componentB.lambdaGibbs.computeDUdlambda;
      const std::size_t upcomingMoleculeIdB = systemB.numberOfMolecules();
      for (Atom& atom : trialAtomsB)
      {
        atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeIdB);
        atom.componentId = static_cast<std::uint8_t>(selectedComponent);
        atom.groupId = groupIdB;
        atom.isFractional = true;
        atom.setScaling(lambdaNewB);
      }

      std::optional<RunningEnergy> insertEnergyB =
          computeMoleculeEnergyDifference(systemB, trialAtomsB, std::span<const Atom>{});
      if (!insertEnergyB.has_value())
      {
        restoreFractionals();
        return std::nullopt;
      }
      energySecondStepB += insertEnergyB.value();

      const double totalSecondStep = (energySecondStepA + energySecondStepB).potentialEnergy();
      rosenbluthOld =
          std::exp(-systemA.beta * totalSecondStep) *
          (static_cast<double>(systemA.numberOfMoleculesPerComponent[selectedComponent] - 1) *
           systemB.simulationBox.volume) /
          (static_cast<double>(systemB.numberOfMoleculesPerComponent[selectedComponent]) *
           systemA.simulationBox.volume);
    }
  }

  componentA.mc_moves_statistics.addConstructed(move);

  double acceptanceProbability = 0.0;
  if (useCBMC)
  {
    acceptanceProbability =
        std::exp(-systemA.beta * (energyFirstStepA->potentialEnergy() + energyFirstStepB->potentialEnergy())) *
        (rosenbluthNewA / rosenbluthOldB) * (rosenbluthNewB / rosenbluthOldA) *
        std::exp((biasNewA - biasOldA) + (biasNewB - biasOldB));
  }
  else
  {
    acceptanceProbability =
        std::exp(-systemA.beta * (energyFirstStepA->potentialEnergy() + energyFirstStepB->potentialEnergy())) *
        rosenbluthNew * rosenbluthOld * std::exp((biasNewA - biasOldA) + (biasNewB - biasOldB));
  }

  if (random.uniform() >= acceptanceProbability)
  {
    restoreFractionals();
    return std::nullopt;
  }

  componentA.mc_moves_statistics.addAccepted(move);
  Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
  Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

  if (moveKindA == GibbsMoveKind::LambdaChange)
  {
    acceptLambdaScaling(systemA, selectedComponent, indexFractionalA, lambdaNewA);
    acceptLambdaScaling(systemB, selectedComponent, indexFractionalB, lambdaNewB);
  }
  else if (moveKindA == GibbsMoveKind::Insert)
  {
    if (useCBMC)
    {
      acceptInsertStep(systemA, selectedComponent, indexFractionalA, growDataA.value(), lambdaNewA);
    }
    else
    {
      acceptInsertStep(systemA, selectedComponent, indexFractionalA, trialInsertA.first, trialInsertA.second,
                       lambdaNewA);
    }
    acceptDeleteStep(systemB, selectedComponent, indexFractionalB, selectedIntegerB, lambdaNewB);
  }
  else
  {
    acceptDeleteStep(systemA, selectedComponent, indexFractionalA, selectedIntegerA, lambdaNewA);
    if (useCBMC)
    {
      acceptInsertStep(systemB, selectedComponent, indexFractionalB, growDataB.value(), lambdaNewB);
    }
    else
    {
      acceptInsertStep(systemB, selectedComponent, indexFractionalB, trialInsertB.first, trialInsertB.second,
                       lambdaNewB);
    }
  }

  RunningEnergy energyDifferenceA = energyFirstStepA.value() + energySecondStepA;
  RunningEnergy energyDifferenceB = energyFirstStepB.value() + energySecondStepB;
  return std::make_pair(energyDifferenceA, energyDifferenceB);
}

}  // namespace

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsConventionalCFCMCMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent)
{
  return gibbsConventionalCFCMCMoveImpl(random, systemA, systemB, selectedComponent, false);
}

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsConventionalCFCMCCBMCMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent)
{
  return gibbsConventionalCFCMCMoveImpl(random, systemA, systemB, selectedComponent, true);
}
