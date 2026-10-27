module;

module mc_moves_reinsertion_multi_trial_cbmc;

import std;

import component;
import atom;
import molecule;
import cbmc;
import cbmc_multi_trial;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import interactions_ewald;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_multi_trial_common;

// Multiple-try reinsertion. Every weight is W~ = W e^{-beta U_F(state)} with U_F the Fourier energy of the
// whole state the trial would produce; the common factor e^{-beta U_F(current)} cancels between the two
// sets, so the energies are measured relative to the CURRENT state throughout. Forward set: M regrowths k
// with W~_k = W_k e^{-beta dU_F(k, o)}; the selected trial is i. Reverse set: the old configuration o itself
// (dU_F = 0) and M-1 regrowths j with W'_j e^{-beta dU_F(j, o)}. The acceptance is the ratio of the two mean
// weights; with M = 1 it reduces to e^{-beta dU_F(i, o)} W_i / W_o, the ordinary reinsertion rule. The tail
// corrections do not change under a reinsertion.

std::optional<RunningEnergy> MC_Moves::reinsertionMoveMultiTrialCBMC(RandomNumber& random, System& system,
                                                                     std::size_t selectedComponent,
                                                                     std::size_t selectedMolecule)
{
  const Move::Types move = Move::Types::ReinsertionMultiTrialCBMC;
  Component& component = system.components[selectedComponent];
  checkMultiTrialSupported(system, "Reinsertion multi-trial (CBMC)");

  component.mc_moves_statistics.addTrial(move);
  if (system.numberOfMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

  // A cross-linked molecule is held by its links; the plain reinsertion handles that case.
  if (system.crossLinks.moleculeIsLinked(selectedComponent, selectedMolecule)) return std::nullopt;

  std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  const std::span<const Atom> oldMolecule = moleculeAtoms;
  Molecule& molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  const CBMC::GrowContext context = system.makeGrowContext();
  const CBMC::MultiTrialGrower& grower = multiTrialGrower(system);
  const std::size_t numberOfTrials = std::max<std::size_t>(1, component.numberOfTrialMolecules);

  // Both sets: the long-range energy of a regrowth k is dU_F(k, o); that of the old configuration is zero.
  const CBMC::LongRangeCorrection longRange = longRangeCorrection(system, oldMolecule, false, 0.0);

  std::optional<CBMC::MultiTrialGrowResult> growData =
      timed(system, component, move, Move::Timing::NonEwald,
            [&] { return grower.regrow(random, context, component, molecule, oldMolecule, numberOfTrials, longRange); });
  if (!growData) return std::nullopt;
  component.mc_moves_statistics.addConstructed(move);

  CBMC::MultiTrialRetraceResult retraceData =
      timed(system, component, move, Move::Timing::NonEwald,
            [&]
            {
              return grower.retraceRegrow(random, context, component, molecule, oldMolecule,
                                          growData->selected.firstBeadStoredR, numberOfTrials, longRange);
            });

  // Re-evaluate the selected trial so that 'trialEik' holds its structure factor.
  std::span<const Atom> newMolecule = std::span(growData->selected.atoms.begin(), growData->selected.atoms.end());
  const LongRangeDifference selectedLongRange =
      timed(system, component, move, Move::Timing::Ewald,
            [&] { return longRangeDifference(system, newMolecule, oldMolecule, false); });

  const double logAcceptance = growData->logMeanWeight - retraceData.logMeanWeight +
                               growData->logProposalCorrection + retraceData.logProposalCorrection;
  if (random.uniform() < std::exp(logAcceptance))
  {
    component.mc_moves_statistics.addAccepted(move);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    std::copy(newMolecule.begin(), newMolecule.end(), moleculeAtoms.begin());
    system.cellListAtomsMoved(moleculeAtoms);
    // the grown record carries the new center of mass and orientation; the position of the molecule in the atom
    // array is unchanged and is not part of a grow result
    const std::size_t atomIndex = molecule.atomIndex;
    molecule = growData->selected.molecule;
    molecule.atomIndex = atomIndex;
    return (growData->selected.energies - retraceData.old.energies) + selectedLongRange.total();
  }
  return std::nullopt;
}
