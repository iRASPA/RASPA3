module;

module mc_moves_pair_insertion_cbmc;

import std;

import double3;
import component;
import molecule;
import atom;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import energy_status;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::pairInsertionMoveCBMC(RandomNumber& random, System& system,
                                                                                  std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Component& componentA = system.components[selectedComponent];

  if (!componentA.pairComponentId.has_value() || !componentA.maximumPairDistance.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::size_t componentB = componentA.pairComponentId.value();
  if (componentB >= system.components.size())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Only the lower-index component performs pair moves to avoid double counting.
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

  componentA.mc_moves_statistics.addTrial(Move::Types::PairSwapCBMC, 0);

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;

  const std::size_t selectedMoleculeA = system.numberOfMolecules();
  const std::size_t selectedMoleculeB = system.numberOfMolecules() + 1;

  time_begin = std::chrono::system_clock::now();
  std::optional<ChainGrowData> growDataA = CBMC::growMoleculeSwapInsertion(
      random, componentA, selectedComponent, system.hasExternalField, system.forceField, system.simulationBox,
      system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), system.beta, componentA.growType, cutOffFrameworkVDW, cutOffMoleculeVDW,
      cutOffCoulomb, selectedMoleculeA, 1.0, false, false);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::NonEwald] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growDataA) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  if (system.insideBlockedPockets(componentA, std::span<const Atom>(growDataA->atom)))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const double r = R_max * random.uniform();
  const double3 direction = random.UnitSphere();
  const double3 firstBeadPositionA = growDataA->atom[componentA.startingBead].position;
  const double3 fixedFirstBeadPositionB = firstBeadPositionA + r * direction;

  std::vector<Atom> moleculeAtomDataWithTrialA(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());
  moleculeAtomDataWithTrialA.insert(moleculeAtomDataWithTrialA.end(), growDataA->atom.begin(), growDataA->atom.end());

  time_begin = std::chrono::system_clock::now();
  std::optional<ChainGrowData> growDataB = CBMC::growMoleculePairSecondSwapInsertion(
      random, componentBRef, componentB, system.hasExternalField, system.forceField, system.simulationBox,
      system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
      moleculeAtomDataWithTrialA, system.beta, componentBRef.growType, cutOffFrameworkVDW, cutOffMoleculeVDW,
      cutOffCoulomb, selectedMoleculeB, fixedFirstBeadPositionB, 1.0, false, false);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::NonEwald] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growDataB) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  if (system.insideBlockedPockets(componentBRef, std::span<const Atom>(growDataB->atom)))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  componentA.mc_moves_statistics.addConstructed(Move::Types::PairSwapCBMC, 0);

  const std::span<const Atom> newMoleculeA = std::span(growDataA->atom.begin(), growDataA->atom.end());
  const std::span<const Atom> newMoleculeB = std::span(growDataB->atom.begin(), growDataB->atom.end());
  std::vector<double3> newElectricFieldA = std::vector<double3>(newMoleculeA.size());
  std::vector<double3> newElectricFieldB = std::vector<double3>(newMoleculeB.size());

  const auto savedStoredEik = system.storedEik;
  const auto savedTotalEik = system.totalEik;

  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMoleculeA, {}, system.netCharge);
  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMoleculeB, {},
      system.netCharge + system.components[selectedComponent].netCharge);
  system.storedEik = savedStoredEik;
  system.totalEik = savedTotalEik;
  RunningEnergy energyFourierDifference = energyFourierDifferenceA + energyFourierDifferenceB;
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::Ewald] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::Ewald] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMoleculeA, {}) +
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMoleculeB, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMoleculeA, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMoleculeB, {});
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::Tail] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwapCBMC][Move::Timing::Tail] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), newElectricFieldA, {},
                                                                  growDataA->atom, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, newElectricFieldA, {}, growDataA->atom, {});

    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), newElectricFieldB, {},
                                                                  growDataB->atom, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, newElectricFieldB, {}, growDataB->atom, {});

    polarizationDifference =
        Interactions::computePolarizationEnergyDifference(system.forceField, newElectricFieldA, {}, growDataA->atom) +
        Interactions::computePolarizationEnergyDifference(system.forceField, newElectricFieldB, {}, growDataB->atom);
  }

  const double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  const double distanceBias = 3.0 * r * r / (R_max * R_max);

  const double fugacityA = componentA.molFraction * componentA.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double fugacityB = componentBRef.molFraction * componentBRef.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasA = componentA.idealGasRosenbluthWeight.value_or(1.0);
  const double idealGasB = componentBRef.idealGasRosenbluthWeight.value_or(1.0);

  const double N_A = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
  const double N_B = double(system.numberOfIntegerMoleculesPerComponent[componentB]);

  const double preFactor = correctionFactorEwald * system.beta * fugacityA * fugacityB * system.simulationBox.volume /
                           ((N_A + 1.0) * (N_B + 1.0)) * distanceBias;

  const double Pacc = preFactor * (growDataA->RosenbluthWeight / idealGasA) * (growDataB->RosenbluthWeight / idealGasB);

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  if (system.tmmc.doTMMC)
  {
    const std::size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    componentA.mc_moves_statistics.addAccepted(Move::Types::PairSwapCBMC, 0);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, newMoleculeA,
                                               {});
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMoleculePolarization(selectedComponent, growDataA->molecule, growDataA->atom, newElectricFieldA);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, newMoleculeB,
                                               {});
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMoleculePolarization(componentB, growDataB->molecule, growDataB->atom, newElectricFieldB);

    return {growDataA->energies + growDataB->energies + energyFourierDifference + tailEnergyDifference +
                polarizationDifference,
            double3(0.0, 1.0 - Pacc, Pacc)};
  }

  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::pairInsertionMove(RandomNumber& random, System& system,
                                                                              std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Component& componentA = system.components[selectedComponent];

  if (!componentA.pairComponentId.has_value() || !componentA.maximumPairDistance.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::size_t componentB = componentA.pairComponentId.value();
  if (componentB >= system.components.size())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Only the lower-index component performs pair moves to avoid double counting.
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

  componentA.mc_moves_statistics.addTrial(Move::Types::PairSwap, 0);

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;

  const std::size_t selectedMoleculeA = system.numberOfMolecules();
  const std::size_t selectedMoleculeB = system.numberOfMolecules() + 1;

  time_begin = std::chrono::system_clock::now();
  std::optional<ChainGrowData> growDataA = CBMC::growMoleculeSwapInsertion(
      random, componentA, selectedComponent, system.hasExternalField, system.forceField, system.simulationBox,
      system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), system.beta, componentA.growType, cutOffFrameworkVDW, cutOffMoleculeVDW,
      cutOffCoulomb, selectedMoleculeA, 1.0, false, false);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::NonEwald] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growDataA) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  if (system.insideBlockedPockets(componentA, std::span<const Atom>(growDataA->atom)))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const double r = R_max * std::cbrt(random.uniform());
  const double3 direction = random.UnitSphere();
  const double3 firstBeadPositionA = growDataA->atom[componentA.startingBead].position;
  const double3 fixedFirstBeadPositionB = firstBeadPositionA + r * direction;

  std::vector<Atom> moleculeAtomDataWithTrialA(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());
  moleculeAtomDataWithTrialA.insert(moleculeAtomDataWithTrialA.end(), growDataA->atom.begin(), growDataA->atom.end());

  time_begin = std::chrono::system_clock::now();
  std::optional<ChainGrowData> growDataB = CBMC::growMoleculePairSecondSwapInsertion(
      random, componentBRef, componentB, system.hasExternalField, system.forceField, system.simulationBox,
      system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
      moleculeAtomDataWithTrialA, system.beta, componentBRef.growType, cutOffFrameworkVDW, cutOffMoleculeVDW,
      cutOffCoulomb, selectedMoleculeB, fixedFirstBeadPositionB, 1.0, false, false);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::NonEwald] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growDataB) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  if (system.insideBlockedPockets(componentBRef, std::span<const Atom>(growDataB->atom)))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  componentA.mc_moves_statistics.addConstructed(Move::Types::PairSwap, 0);

  const std::span<const Atom> newMoleculeA = std::span(growDataA->atom.begin(), growDataA->atom.end());
  const std::span<const Atom> newMoleculeB = std::span(growDataB->atom.begin(), growDataB->atom.end());
  std::vector<double3> newElectricFieldA = std::vector<double3>(newMoleculeA.size());
  std::vector<double3> newElectricFieldB = std::vector<double3>(newMoleculeB.size());

  const auto savedStoredEik = system.storedEik;
  const auto savedTotalEik = system.totalEik;

  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMoleculeA, {}, system.netCharge);
  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMoleculeB, {},
      system.netCharge + system.components[selectedComponent].netCharge);
  system.storedEik = savedStoredEik;
  system.totalEik = savedTotalEik;
  RunningEnergy energyFourierDifference = energyFourierDifferenceA + energyFourierDifferenceB;
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::Ewald] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::Ewald] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMoleculeA, {}) +
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMoleculeB, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMoleculeA, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMoleculeB, {});
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::Tail] += (time_end - time_begin);
  componentA.mc_moves_cputime[Move::Types::PairSwap][Move::Timing::Tail] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), newElectricFieldA, {},
                                                                  growDataA->atom, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, newElectricFieldA, {}, growDataA->atom, {});

    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), newElectricFieldB, {},
                                                                  growDataB->atom, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, newElectricFieldB, {}, growDataB->atom, {});

    polarizationDifference =
        Interactions::computePolarizationEnergyDifference(system.forceField, newElectricFieldA, {}, growDataA->atom) +
        Interactions::computePolarizationEnergyDifference(system.forceField, newElectricFieldB, {}, growDataB->atom);
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

  const double preFactor = correctionFactorEwald * system.beta * fugacityA * fugacityB * system.simulationBox.volume /
                           ((N_A + 1.0) * (N_B + 1.0)) * distanceBias;

  const double Pacc = preFactor * (growDataA->RosenbluthWeight / idealGasA) * (growDataB->RosenbluthWeight / idealGasB);

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  if (system.tmmc.doTMMC)
  {
    const std::size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    componentA.mc_moves_statistics.addAccepted(Move::Types::PairSwap, 0);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, newMoleculeA,
                                               {});
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMoleculePolarization(selectedComponent, growDataA->molecule, growDataA->atom, newElectricFieldA);

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik,
                                               system.totalEik, system.forceField, system.simulationBox, newMoleculeB,
                                               {});
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMoleculePolarization(componentB, growDataB->molecule, growDataB->atom, newElectricFieldB);

    return {growDataA->energies + growDataB->energies + energyFourierDifference + tailEnergyDifference +
                polarizationDifference,
            double3(0.0, 1.0 - Pacc, Pacc)};
  }

  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
