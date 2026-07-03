module;

module mc_moves_pair_deletion_cbmc;

import std;

import double3;
import component;
import atom;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::pairDeletionMoveCBMC(RandomNumber& random, System& system,
                                                                                 std::size_t selectedComponent,
                                                                                 std::size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Component& componentA = system.components[selectedComponent];

  if (!componentA.pairComponentId.has_value() || !componentA.maximumPairDistance.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::size_t componentB = componentA.pairComponentId.value();
  if (selectedComponent > componentB)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  Component& componentBRef = system.components[componentB];
  const double R_max = componentA.maximumPairDistance.value();
  if (R_max <= 0.0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  componentA.mc_moves_statistics.addTrial(Move::Types::PairSwapCBMC, 1);

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  if (system.numberOfIntegerMoleculesPerComponent[componentB] == 0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  std::span<Atom> moleculeAAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  const double3 positionA = moleculeAAtoms[system.components[selectedComponent].startingBead].position;

  // integer molecules are stored after the fractional molecules of a component
  const std::size_t firstIntegerMoleculeB = system.numberOfFractionalMoleculesPerComponent[componentB];

  // the partner is the closest integer molecule of the paired component within R_max
  std::optional<std::size_t> selectedPartner;
  double closestDistance = R_max;
  for (std::size_t moleculeB = firstIntegerMoleculeB;
       moleculeB < system.numberOfMoleculesPerComponent[componentB]; ++moleculeB)
  {
    std::span<Atom> moleculeBAtoms = system.spanOfMolecule(componentB, moleculeB);
    const double distance =
        (positionA - moleculeBAtoms[system.components[componentB].startingBead].position).length();
    if (distance <= closestDistance)
    {
      closestDistance = distance;
      selectedPartner = moleculeB;
    }
  }

  if (!selectedPartner.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::size_t selectedMoleculeB = selectedPartner.value();
  std::span<Atom> moleculeA = system.spanOfMolecule(selectedComponent, selectedMolecule);
  std::span<Atom> moleculeB = system.spanOfMolecule(componentB, selectedMoleculeB);

  const double r = (moleculeA[componentA.startingBead].position - moleculeB[componentBRef.startingBead].position).length();
  if (r <= 0.0 || r > R_max)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;

  std::vector<Atom> backgroundWithoutPair;
  backgroundWithoutPair.reserve(system.spanOfMoleculeAtoms().size());

  for (std::size_t component = 0; component < system.components.size(); ++component)
  {
    // include all molecules (fractional ones as well) except the pair that is being deleted
    for (std::size_t molecule = 0; molecule < system.numberOfMoleculesPerComponent[component]; ++molecule)
    {
      if ((component == selectedComponent && molecule == selectedMolecule) ||
          (component == componentB && molecule == selectedMoleculeB))
      {
        continue;
      }
      std::span<Atom> atoms = system.spanOfMolecule(component, molecule);
      backgroundWithoutPair.insert(backgroundWithoutPair.end(), atoms.begin(), atoms.end());
    }
  }

  std::vector<Atom> backgroundForB = backgroundWithoutPair;
  backgroundForB.insert(backgroundForB.end(), moleculeA.begin(), moleculeA.end());

  time_begin = std::chrono::system_clock::now();
  ChainRetraceData retraceDataA = CBMC::retraceMoleculeSwapDeletion(
      random, componentA, system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
      system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(), backgroundWithoutPair,
      system.beta, componentA.growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, moleculeA);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC]["NonEwald"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC]["NonEwald"] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  ChainRetraceData retraceDataB = CBMC::retraceMoleculePairSecondSwapDeletion(
      componentBRef, system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
      system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(), backgroundForB,
      system.beta, componentBRef.growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, moleculeB);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC]["NonEwald"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC]["NonEwald"] += (time_end - time_begin);

  const std::span<const Atom> oldMoleculeA = std::span<const Atom>(moleculeA.data(), moleculeA.size());
  const std::span<const Atom> oldMoleculeB = std::span<const Atom>(moleculeB.data(), moleculeB.size());

  const auto savedStoredEik = system.storedEik;
  const auto savedTotalEik = system.totalEik;

  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, {}, oldMoleculeB, system.netCharge);
  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, {}, oldMoleculeA,
      system.netCharge - system.components[componentB].netCharge);
  system.storedEik = savedStoredEik;
  system.totalEik = savedTotalEik;
  RunningEnergy energyFourierDifference = energyFourierDifferenceA + energyFourierDifferenceB;
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC]["Ewald"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC]["Ewald"] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), {}, oldMoleculeA) +
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), {}, oldMoleculeB) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), {}, oldMoleculeA) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), {}, oldMoleculeB);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC]["Tail"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC]["Tail"] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    std::vector<double3> oldElectricFieldA = std::vector<double3>(moleculeA.size());
    std::vector<double3> oldElectricFieldB = std::vector<double3>(moleculeB.size());

    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), {}, oldElectricFieldA,
                                                                  {}, oldMoleculeA);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, {}, oldElectricFieldA, {}, oldMoleculeA);

    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), {}, oldElectricFieldB,
                                                                  {}, oldMoleculeB);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, {}, oldElectricFieldB, {}, oldMoleculeB);

    polarizationDifference =
        Interactions::computePolarizationEnergyDifference(system.forceField, {}, oldElectricFieldA, {}, moleculeA) +
        Interactions::computePolarizationEnergyDifference(system.forceField, {}, oldElectricFieldB, {}, moleculeB);
  }

  const double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  const double distanceBias = (R_max * R_max) / (3.0 * r * r);

  const double fugacityA = componentA.molFraction * componentA.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double fugacityB = componentBRef.molFraction * componentBRef.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasA = componentA.idealGasRosenbluthWeight.value_or(1.0);
  const double idealGasB = componentBRef.idealGasRosenbluthWeight.value_or(1.0);

  const double N_A = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
  const double N_B = double(system.numberOfIntegerMoleculesPerComponent[componentB]);
  const double rosenbluthWeight = retraceDataA.RosenbluthWeight * retraceDataB.RosenbluthWeight;

  const double preFactor = correctionFactorEwald * (N_A * N_B) /
                           (system.beta * fugacityA * fugacityB * system.simulationBox.volume * distanceBias);

  const double Pacc = preFactor * (idealGasA * idealGasB) / rosenbluthWeight;

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

  if (system.tmmc.doTMMC)
  {
    const std::size_t newN = oldN - 1;
    if (newN < system.tmmc.minMacrostate)
    {
      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
    }
  }

  componentA.mc_moves_statistics.addConstructed(Move::Types::PairSwapCBMC, 1);

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    componentA.mc_moves_statistics.addAccepted(Move::Types::PairSwapCBMC, 1);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, {}, oldMoleculeB);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.deleteMolecule(componentB, selectedMoleculeB, moleculeB);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, {}, oldMoleculeA);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.deleteMolecule(selectedComponent, selectedMolecule, moleculeA);

    return {retraceDataA.energies + retraceDataB.energies - energyFourierDifference - tailEnergyDifference -
                polarizationDifference,
            double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::pairDeletionMove(RandomNumber& random, System& system,
                                                                             std::size_t selectedComponent,
                                                                             std::size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Component& componentA = system.components[selectedComponent];

  if (!componentA.pairComponentId.has_value() || !componentA.maximumPairDistance.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::size_t componentB = componentA.pairComponentId.value();
  if (selectedComponent > componentB)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  Component& componentBRef = system.components[componentB];
  const double R_max = componentA.maximumPairDistance.value();
  if (R_max <= 0.0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  componentA.mc_moves_statistics.addTrial(Move::Types::PairSwap, 1);

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  if (system.numberOfIntegerMoleculesPerComponent[componentB] == 0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  std::span<Atom> moleculeAAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  const double3 positionA = moleculeAAtoms[system.components[selectedComponent].startingBead].position;

  // integer molecules are stored after the fractional molecules of a component
  const std::size_t firstIntegerMoleculeB = system.numberOfFractionalMoleculesPerComponent[componentB];

  // the partner is the closest integer molecule of the paired component within R_max
  std::optional<std::size_t> selectedPartner;
  double closestDistance = R_max;
  for (std::size_t moleculeB = firstIntegerMoleculeB;
       moleculeB < system.numberOfMoleculesPerComponent[componentB]; ++moleculeB)
  {
    std::span<Atom> moleculeBAtoms = system.spanOfMolecule(componentB, moleculeB);
    const double distance =
        (positionA - moleculeBAtoms[system.components[componentB].startingBead].position).length();
    if (distance <= closestDistance)
    {
      closestDistance = distance;
      selectedPartner = moleculeB;
    }
  }

  if (!selectedPartner.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::size_t selectedMoleculeB = selectedPartner.value();
  std::span<Atom> moleculeA = system.spanOfMolecule(selectedComponent, selectedMolecule);
  std::span<Atom> moleculeB = system.spanOfMolecule(componentB, selectedMoleculeB);

  const double r = (moleculeA[componentA.startingBead].position - moleculeB[componentBRef.startingBead].position).length();
  if (r <= 0.0 || r > R_max)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;

  std::vector<Atom> backgroundWithoutPair;
  backgroundWithoutPair.reserve(system.spanOfMoleculeAtoms().size());

  for (std::size_t component = 0; component < system.components.size(); ++component)
  {
    // include all molecules (fractional ones as well) except the pair that is being deleted
    for (std::size_t molecule = 0; molecule < system.numberOfMoleculesPerComponent[component]; ++molecule)
    {
      if ((component == selectedComponent && molecule == selectedMolecule) ||
          (component == componentB && molecule == selectedMoleculeB))
      {
        continue;
      }
      std::span<Atom> atoms = system.spanOfMolecule(component, molecule);
      backgroundWithoutPair.insert(backgroundWithoutPair.end(), atoms.begin(), atoms.end());
    }
  }

  std::vector<Atom> backgroundForB = backgroundWithoutPair;
  backgroundForB.insert(backgroundForB.end(), moleculeA.begin(), moleculeA.end());

  time_begin = std::chrono::system_clock::now();
  ChainRetraceData retraceDataA = CBMC::retraceMoleculeSwapDeletion(
      random, componentA, system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
      system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(), backgroundWithoutPair,
      system.beta, componentA.growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, moleculeA);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap]["NonEwald"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap]["NonEwald"] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  ChainRetraceData retraceDataB = CBMC::retraceMoleculePairSecondSwapDeletion(
      componentBRef, system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
      system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(), backgroundForB,
      system.beta, componentBRef.growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, moleculeB);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap]["NonEwald"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap]["NonEwald"] += (time_end - time_begin);

  const std::span<const Atom> oldMoleculeA = std::span<const Atom>(moleculeA.data(), moleculeA.size());
  const std::span<const Atom> oldMoleculeB = std::span<const Atom>(moleculeB.data(), moleculeB.size());

  const auto savedStoredEik = system.storedEik;
  const auto savedTotalEik = system.totalEik;

  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, {}, oldMoleculeB, system.netCharge);
  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, {}, oldMoleculeA,
      system.netCharge - system.components[componentB].netCharge);
  system.storedEik = savedStoredEik;
  system.totalEik = savedTotalEik;
  RunningEnergy energyFourierDifference = energyFourierDifferenceA + energyFourierDifferenceB;
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap]["Ewald"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap]["Ewald"] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), {}, oldMoleculeA) +
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), {}, oldMoleculeB) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), {}, oldMoleculeA) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), {}, oldMoleculeB);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap]["Tail"] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap]["Tail"] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    std::vector<double3> oldElectricFieldA = std::vector<double3>(moleculeA.size());
    std::vector<double3> oldElectricFieldB = std::vector<double3>(moleculeB.size());

    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), {}, oldElectricFieldA,
                                                                  {}, oldMoleculeA);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, {}, oldElectricFieldA, {}, oldMoleculeA);

    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), {}, oldElectricFieldB,
                                                                  {}, oldMoleculeB);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, {}, oldElectricFieldB, {}, oldMoleculeB);

    polarizationDifference =
        Interactions::computePolarizationEnergyDifference(system.forceField, {}, oldElectricFieldA, {}, moleculeA) +
        Interactions::computePolarizationEnergyDifference(system.forceField, {}, oldElectricFieldB, {}, moleculeB);
  }

  const double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  const double distanceBias = 1.0;

  const double fugacityA = componentA.molFraction * componentA.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double fugacityB = componentBRef.molFraction * componentBRef.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasA = componentA.idealGasRosenbluthWeight.value_or(1.0);
  const double idealGasB = componentBRef.idealGasRosenbluthWeight.value_or(1.0);

  const double N_A = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
  const double N_B = double(system.numberOfIntegerMoleculesPerComponent[componentB]);
  const double rosenbluthWeight = retraceDataA.RosenbluthWeight * retraceDataB.RosenbluthWeight;

  const double preFactor = correctionFactorEwald * (N_A * N_B) /
                           (system.beta * fugacityA * fugacityB * system.simulationBox.volume * distanceBias);

  const double Pacc = preFactor * (idealGasA * idealGasB) / rosenbluthWeight;

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

  if (system.tmmc.doTMMC)
  {
    const std::size_t newN = oldN - 1;
    if (newN < system.tmmc.minMacrostate)
    {
      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
    }
  }

  componentA.mc_moves_statistics.addConstructed(Move::Types::PairSwap, 1);

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    componentA.mc_moves_statistics.addAccepted(Move::Types::PairSwap, 1);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, {}, oldMoleculeB);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.deleteMolecule(componentB, selectedMoleculeB, moleculeB);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, {}, oldMoleculeA);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.deleteMolecule(selectedComponent, selectedMolecule, moleculeA);

    return {retraceDataA.energies + retraceDataB.energies - energyFourierDifference - tailEnergyDifference -
                polarizationDifference,
            double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
}
