module;

module mc_moves_gibbs_swap_multi_trial_cbmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;
import atom;
import cbmc;
import cbmc_multi_trial;
import forcefield;
import component;
import simulationbox;
import interactions_ewald;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_multi_trial_common;

namespace
{
// The transition-matrix bookkeeping of a two-box transfer (mirrors the one of the plain Gibbs swap).
class DualTMMCTrial
{
 public:
  DualTMMCTrial(System& systemA, System& systemB, std::size_t selectedComponent)
      : systemA_(systemA),
        systemB_(systemB),
        oldNA_(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]),
        oldNB_(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent])
  {
  }

  ~DualTMMCTrial()
  {
    if (!recorded_)
    {
      systemA_.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldNA_);
      systemB_.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldNB_);
    }
  }

  bool transferIsInBounds() const
  {
    return (!systemA_.tmmc.doTMMC || !systemA_.tmmc.rejectOutOfBound || oldNA_ < systemA_.tmmc.maxMacrostate) &&
           (!systemB_.tmmc.doTMMC || !systemB_.tmmc.rejectOutOfBound ||
            (oldNB_ > 0 && oldNB_ - 1 >= systemB_.tmmc.minMacrostate));
  }

  double biasFactor() const
  {
    return systemA_.tmmc.biasFactor(oldNA_ + 1, oldNA_) * systemB_.tmmc.biasFactor(oldNB_ - 1, oldNB_);
  }

  void recordTransfer(double physicalAcceptance)
  {
    systemA_.tmmc.updateMatrix(double3(0.0, 1.0 - physicalAcceptance, physicalAcceptance), oldNA_);
    systemB_.tmmc.updateMatrix(double3(physicalAcceptance, 1.0 - physicalAcceptance, 0.0), oldNB_);
    recorded_ = true;
  }

 private:
  System& systemA_;
  System& systemB_;
  std::size_t oldNA_;
  std::size_t oldNB_;
  bool recorded_{false};
};
}  // namespace

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMoveMultiTrialCBMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent)
{
  const Move::Types move = Move::Types::GibbsSwapMultiTrialCBMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];
  checkMultiTrialSupported(systemA, "Gibbs swap multi-trial (CBMC)");
  DualTMMCTrial tmmcTrial(systemA, systemB, selectedComponent);

  if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

  componentA.mc_moves_statistics.addTrial(move);
  componentB.mc_moves_statistics.addTrial(move);

  const std::size_t numberOfTrials = std::max<std::size_t>(1, componentA.numberOfTrialMolecules);

  // Insertion side (system A): M trials grown into the unchanged system A.
  const std::size_t newMoleculeIndex = systemA.numberOfMolecules();
  const CBMC::GrowContext growContext = systemA.makeGrowContext();
  const CBMC::LongRangeCorrection longRangeA = longRangeCorrection(systemA, {}, true, 0.0);

  std::optional<CBMC::MultiTrialGrowResult> growData =
      timed(systemA, componentA, move, Move::Timing::NonEwald,
            [&]
            {
              return multiTrialGrower(systemA).growNew(random, growContext, componentA,
                                                       {.componentId = selectedComponent, .moleculeId = newMoleculeIndex},
                                                       numberOfTrials, longRangeA);
            });
  if (!growData) return std::nullopt;
  componentA.mc_moves_statistics.addConstructed(move);

  std::span<const Atom> newMolecule = std::span(growData->selected.atoms.begin(), growData->selected.atoms.end());
  const LongRangeDifference longRangeDifferenceA =
      timed(systemA, componentA, move, Move::Timing::Ewald,
            [&] { return longRangeDifference(systemA, newMolecule, {}, true); });

  // Deletion side (system B): the molecule retraced and M-1 references grown into B without it.
  const std::size_t selectedMolecule = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
  std::span<Atom> molecule = systemB.spanOfMolecule(selectedComponent, selectedMolecule);
  const std::span<const Atom> oldMolecule = molecule;
  const std::size_t oldGlobalMoleculeId = systemB.moleculeIndexOfComponent(selectedComponent, selectedMolecule);

  const CBMC::GrowContext retraceContext = systemB.makeGrowContext();
  const double baselineB =
      timed(systemA, componentA, move, Move::Timing::Ewald,
            [&] { return longRangeDifference(systemB, {}, oldMolecule, true).potentialEnergy(); });
  const CBMC::LongRangeCorrection longRangeB = longRangeCorrection(systemB, oldMolecule, true, baselineB);

  CBMC::MultiTrialRetraceResult retraceData =
      timed(systemA, componentA, move, Move::Timing::NonEwald,
            [&]
            {
              return multiTrialGrower(systemB).retraceNew(
                  random, retraceContext, componentB, oldMolecule,
                  {.componentId = selectedComponent, .moleculeId = oldGlobalMoleculeId}, numberOfTrials, longRangeB);
            });
  componentB.mc_moves_statistics.addConstructed(move);

  const LongRangeDifference longRangeDifferenceB =
      timed(systemA, componentA, move, Move::Timing::Ewald,
            [&] { return longRangeDifference(systemB, {}, oldMolecule, true); });

  const double physicalAcceptance = std::exp(
      growData->logMeanWeight - retraceData.logMeanWeight + growData->logProposalCorrection +
      retraceData.logProposalCorrection +
      std::log((static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) *
                systemA.simulationBox.volume) /
               ((1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent])) *
                systemB.simulationBox.volume)));

  if (!tmmcTrial.transferIsInBounds())
  {
    tmmcTrial.recordTransfer(physicalAcceptance);
    return std::nullopt;
  }
  const double tmmcBias = tmmcTrial.biasFactor();
  tmmcTrial.recordTransfer(physicalAcceptance);

  if (random.uniform() < tmmcBias * physicalAcceptance)
  {
    componentA.mc_moves_statistics.addAccepted(move);
    Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
    systemA.insertMolecule(selectedComponent, growData->selected.molecule, growData->selected.atoms);

    componentB.mc_moves_statistics.addAccepted(move);
    Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);
    systemB.deleteMolecule(selectedComponent, selectedMolecule, molecule);

    return std::make_pair(growData->selected.energies + longRangeDifferenceA.total(),
                          -(retraceData.old.energies - longRangeDifferenceB.total()));
  }
  return std::nullopt;
}
