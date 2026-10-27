module;

module mc_moves_swap_multi_trial_cbmc;

import std;

import component;
import molecule;
import atom;
import double3;
import simulationbox;
import cbmc;
import cbmc_multi_trial;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_ewald;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_multi_trial_common;

// Multiple-try CBMC swap. Insertion:
//
//   P_acc = min(1, beta f V / (N+1) x W-bar_new / W^IG),   W-bar_new = (1/M) sum_i W_i e^{-beta dU_lr,i}
//
// with the trial selected with probability W~_i / sum_j W~_j. Deletion:
//
//   P_acc = min(1, N / (beta f V) x W^IG / W-bar_old),     W-bar_old = (1/M) ( W~_old + sum_{j<M} W~'_j )
//
// where the M-1 reference molecules are grown into the system WITHOUT the molecule being deleted, and every
// long-range term is measured relative to that state. The long-range factors are inside the mean weights:
// no separate Ewald correction factor multiplies the acceptance (with M = 1 the two are the same thing).

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMoveMultiTrialCBMC(RandomNumber& random,
                                                                                       System& system,
                                                                                       std::size_t selectedComponent)
{
  const Move::Types move = Move::Types::SwapMultiTrialCBMC;
  Component& component = system.components[selectedComponent];
  checkMultiTrialSupported(system, "Swap multi-trial (CBMC)");

  // The trial molecule gets an id beyond every existing molecule; 'insertMolecule' places it after the last
  // molecule of the component.
  const std::size_t trialMoleculeId = system.numberOfMolecules();
  component.mc_moves_statistics.addTrial(move, 0);

  const CBMC::GrowContext growContext = system.makeGrowContext();
  const CBMC::MultiTrialGrower& grower = multiTrialGrower(system);
  const std::size_t numberOfTrials = std::max<std::size_t>(1, component.numberOfTrialMolecules);

  // Insertion into the unchanged system: the long-range energy of each trial is dU(trial, {}).
  const CBMC::LongRangeCorrection longRange = longRangeCorrection(system, {}, true, 0.0);

  std::optional<CBMC::MultiTrialGrowResult> growData =
      timed(system, component, move, Move::Timing::NonEwald,
            [&]
            {
              return grower.growNew(random, growContext, component,
                                    {.componentId = selectedComponent, .moleculeId = trialMoleculeId}, numberOfTrials,
                                    longRange);
            });
  if (!growData) return {std::nullopt, double3(0.0, 1.0, 0.0)};
  component.mc_moves_statistics.addConstructed(move, 0);

  // Re-evaluate the selected trial so that 'trialEik' holds its structure factor (the callback left the last
  // trial's) and the running energies get the separate Fourier and tail terms.
  std::span<const Atom> newMolecule = std::span(growData->selected.atoms.begin(), growData->selected.atoms.end());
  const LongRangeDifference selectedLongRange =
      timed(system, component, move, Move::Timing::Ewald,
            [&] { return longRangeDifference(system, newMolecule, {}, true); });

  const double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const double logPreFactor =
      std::log(system.beta * fugacity * system.simulationBox.volume / static_cast<double>(oldN + 1));
  // 'logProposalCorrection' is zero for the CPU grower; a device grower that proposes with its own (single
  // precision) potential adds the ratio of the target and proposal Boltzmann factors of the selected trial here.
  const double Pacc = std::exp(logPreFactor + growData->logMeanWeight - std::log(idealGasRosenbluthWeight) +
                               growData->logProposalCorrection);

  if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN >= system.tmmc.maxMacrostate)
  {
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }
  const std::size_t newN = oldN == std::numeric_limits<std::size_t>::max() ? oldN : oldN + 1;
  const double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    component.mc_moves_statistics.addAccepted(move, 0);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    system.insertMolecule(selectedComponent, growData->selected.molecule, growData->selected.atoms);
    return {growData->selected.energies + selectedLongRange.total(), double3(0.0, 1.0 - Pacc, Pacc)};
  }
  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMoveMultiTrialCBMC(RandomNumber& random,
                                                                                      System& system,
                                                                                      std::size_t selectedComponent,
                                                                                      std::size_t selectedMolecule)
{
  const Move::Types move = Move::Types::SwapMultiTrialCBMC;
  Component& component = system.components[selectedComponent];
  checkMultiTrialSupported(system, "Swap multi-trial (CBMC)");

  component.mc_moves_statistics.addTrial(move, 1);
  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // A cross-linked molecule can not leave the system while it is bonded to another molecule.
  if (system.crossLinks.moleculeIsLinked(selectedComponent, selectedMolecule))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);
  const std::span<const Atom> oldMolecule = molecule;
  const std::size_t oldGlobalMoleculeId = system.moleculeIndexOfComponent(selectedComponent, selectedMolecule);

  const CBMC::GrowContext context = system.makeGrowContext();
  const CBMC::MultiTrialGrower& grower = multiTrialGrower(system);
  const std::size_t numberOfTrials = std::max<std::size_t>(1, component.numberOfTrialMolecules);

  // The reference state is the system without the molecule: the long-range energy of a reference trial is
  // dU(trial, old) - dU({}, old), and that of the old molecule itself -dU({}, old) (its own contribution).
  const double baseline =
      timed(system, component, move, Move::Timing::Ewald,
            [&] { return longRangeDifference(system, {}, oldMolecule, true).potentialEnergy(); });
  const CBMC::LongRangeCorrection longRange = longRangeCorrection(system, oldMolecule, true, baseline);

  // The reference molecules carry the id of the molecule being deleted, so the same-id rule of the grow
  // context excludes it from their background.
  CBMC::MultiTrialRetraceResult retraceData =
      timed(system, component, move, Move::Timing::NonEwald,
            [&]
            {
              return grower.retraceNew(random, context, component, oldMolecule,
                                       {.componentId = selectedComponent, .moleculeId = oldGlobalMoleculeId},
                                       numberOfTrials, longRange);
            });
  component.mc_moves_statistics.addConstructed(move, 1);

  // Restore 'trialEik' to the structure factor without the molecule, and get the separate terms for the
  // running energies.
  const LongRangeDifference oldLongRange =
      timed(system, component, move, Move::Timing::Ewald, [&] { return longRangeDifference(system, {}, oldMolecule, true); });

  const double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const double logPreFactor =
      std::log(static_cast<double>(oldN) / (system.beta * fugacity * system.simulationBox.volume));
  const double Pacc = std::exp(logPreFactor + std::log(idealGasRosenbluthWeight) - retraceData.logMeanWeight +
                               retraceData.logProposalCorrection);

  if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
  {
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }
  const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
  const double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    component.mc_moves_statistics.addAccepted(move, 1);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    system.deleteMolecule(selectedComponent, selectedMolecule, molecule);
    return {retraceData.old.energies - oldLongRange.total(), double3(Pacc, 1.0 - Pacc, 0.0)};
  }
  return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
}
