module;

module mc_moves_gibbs_swap_cbcfcmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;
import molecule;
import atom;
import cbmc;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import component;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_gibbs_swap_serial_common;

// All systems have a fractional molecule, only one of these is 'active', the others are switched off with 'lambda=0'.
// Serial CFCMC: at most one simulation box has containsTheFractionalMolecule==true; systemA must be that box.
// The interchange retraces the fractional molecule and grows a new integer molecule in the active box, and retraces
// an integer molecule and grows a new fractional molecule in the other box; the shuffle retraces the fractional
// molecule in the active box and grows it in the other box. Every CBMC operation has its mirror image in the reverse
// move, so the acceptance rules are ratios of Rosenbluth weights.

using namespace MC_Moves::GibbsSwapSerialCommon;

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

  // A transfer attempt that exits before its acceptance probability is known (overlap, empty box) is a
  // rejected transfer for both boxes. The lambda-change and reinsertion sub-moves leave N unchanged in both
  // boxes and record nothing (see the collection-matrix policy in mc_moves.cpp).
  ~DualTMMCTrial()
  {
    if (transferAttempted_ && !recorded_)
    {
      systemA_.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldNA_);
      systemB_.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldNB_);
    }
  }

  void beginTransfer() { transferAttempted_ = true; }

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
  bool transferAttempted_{false};
  bool recorded_{false};
};

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CBCFCMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent,
    [[maybe_unused]] std::size_t& fractionalMoleculeSystem)
{
  Move::Types move = Move::Types::GibbsSwapCBCFCMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];
  DualTMMCTrial tmmcTrial(systemA, systemB, selectedComponent);

  PropertyLambdaProbabilityHistogram& lambdaA = componentA.lambdaGC;
  PropertyLambdaProbabilityHistogram& lambdaB = componentB.lambdaGC;
  std::size_t oldBin = lambdaA.currentBin;
  double deltaLambda = lambdaA.delta;
  double oldLambda = componentA.lambdaGC.lambdaValue();

  double maxChange = componentA.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambdaA.selectNewBin(random, maxChange);

  double switchValue = random.uniform();

  if (!systemA.containsTheFractionalMolecule || systemB.containsTheFractionalMolecule)
  {
    return std::nullopt;
  }

  if (switchValue < 0.25)
  {
    // Lambda interchange: the fractional molecule of A becomes an integer molecule and an integer molecule of B
    // becomes the fractional molecule, both through CBMC. In A the fractional molecule is retraced at lambda and a
    // new integer molecule is grown; in B a randomly selected integer molecule is retraced and a new fractional
    // molecule is grown at the same lambda. The reverse move (with the roles of the boxes exchanged) regenerates
    // the original configuration, so the acceptance rule is the ratio of the four Rosenbluth weights times the
    // Boltzmann factor of the terms CBMC does not see (Ewald, tail, other-component flag swaps). Every grow and
    // retrace is evaluated against a background of integer molecules only (the fractional molecules of the other
    // components are off in an inactive box), as the matching operation of the reverse move is.

    componentA.mc_moves_statistics.addTrial(move, 0);
    tmmcTrial.beginTransfer();

    if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::size_t indexFractionalMoleculeB =
        systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    const std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    const std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    const FractionalSnapshot snapshotA = saveFractionalMolecules(systemA);
    const FractionalSnapshot snapshotB = saveFractionalMolecules(systemB);

    // Snapshot effective tail-correction counts for both systems; threaded through the sequential sub-steps.
    std::vector<double> tailEffA = systemA.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupA =
        systemA.fractionalPseudoAtomCountsPerGroup;
    std::vector<double> tailEffB = systemB.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupB =
        systemB.fractionalPseudoAtomCountsPerGroup;
    // Coulomb-scaled net charge of each box, threaded across the chained Ewald sub-steps like the tail counts.
    double netChargeA = systemA.netCharge;
    double netChargeB = systemB.netCharge;

    // System A, step 1: Ewald and tail contributions of removing the fractional molecule at lambda (the first Ewald
    // sub-step of the move in A: storedEik -> trialEik). Its short-range and intra-molecular energies come from the
    // CBMC retrace below.
    RunningEnergy EwaldFourierDifferenceRetraceA = timed(
        systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
        [&]
        {
          return Interactions::energyDifferenceEwaldFourier(
              systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.trialEik,
              systemA.forceField, systemA.simulationBox, systemA.components, {}, oldFractionalMoleculeA, netChargeA);
        });
    netChargeA += Interactions::scaledNetChargeDifference({}, oldFractionalMoleculeA);

    RunningEnergy tailEnergyDifferenceRetraceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&] { return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, {}, oldFractionalMoleculeA); });

    // System A, step 2: the fractional molecules of the other components switch off (unbiased). Evaluated while the
    // selected fractional molecule is still in memory at lambda, so that its cross terms with them are counted here;
    // the retrace and grow below then see a background of integer molecules only.
    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return computeFlagSwapOtherComponentsEnergyDifference(systemA, tailEffA, tailGroupA, netChargeA,
                                                                      snapshotA, selectedComponent, false, systemA);
              });
    if (!flagSwapOtherComponentsDifferenceA.has_value())
    {
      return std::nullopt;
    }

    // System A, step 3: the fractional slot takes the inactive placeholder of B, all fractional molecules of A take
    // the inactive scaling in memory, and the old fractional molecule is retraced at lambda against that background
    // (its own slot is excluded by molecule id).
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.isFractional = true;
      atom.moleculeId =
          static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }
    applyFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation of B enters A; the old
    // fractional conformation leaves A through the retrace energies. Not part of the acceptance rule.
    RunningEnergy intraEnergyPlaceholderA = componentA.intraMolecularPotentials.computeInternalEnergies(
        systemA.forceField, systemA.simulationBox, fractionalMoleculeA);

    const CBMC::GrowContext retraceContextA = systemA.makeGrowContext();
    CBMC::RetraceResult retraceDataA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&] { return CBMC::retraceMolecule(random, retraceContextA, componentA, oldFractionalMoleculeA); });

    // System A, step 4: CBMC grow a new integer molecule against the same background.
    // Trial global molecule id (unique, not yet in the system); insertMolecule assigns the final id.
    std::size_t newMoleculeIndex = systemA.numberOfMolecules();
    const CBMC::GrowContext growContextA = systemA.makeGrowContext();
    std::optional<CBMC::GrowResult> growDataA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return CBMC::growNewMolecule(random, growContextA, componentA,
                                             {.componentId = selectedComponent, .moleculeId = newMoleculeIndex});
              });

    if (!growDataA)
    {
      restoreFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceGrowA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.trialEik, systemA.trialEik,
                    systemA.forceField, systemA.simulationBox, systemA.components,
                    std::span(growDataA->atoms.begin(), growDataA->atoms.end()), {}, netChargeA);
              });
    netChargeA += Interactions::scaledNetChargeDifference(growDataA->atoms, {});

    RunningEnergy tailEnergyDifferenceGrowA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemA, tailEffA, tailGroupA,
                                                   std::span(growDataA->atoms.begin(), growDataA->atoms.end()), {});
              });

    // Terms of A that enter the acceptance rule explicitly (everything the retrace and the grow do not see).
    RunningEnergy unbiasedDifferenceA = EwaldFourierDifferenceRetraceA + tailEnergyDifferenceRetraceA +
                                        flagSwapOtherComponentsDifferenceA.value() + EwaldFourierDifferenceGrowA +
                                        tailEnergyDifferenceGrowA;

    RunningEnergy energyDifferenceA =
        -(retraceDataA.energies - EwaldFourierDifferenceRetraceA - tailEnergyDifferenceRetraceA) +
        flagSwapOtherComponentsDifferenceA.value() + intraEnergyPlaceholderA + growDataA->energies +
        EwaldFourierDifferenceGrowA + tailEnergyDifferenceGrowA;

    // System B, step 1: CBMC retrace a randomly selected integer molecule (B is inactive: its fractional molecules
    // are off, so the background is the integer molecules only).
    std::size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);

    // A cross-linked molecule can not become the fractional molecule (it leaves the box when lambda reaches zero,
    // and only integer molecules carry links).
    if (systemB.crossLinks.moleculeIsLinked(selectedComponent, indexSelectedIntegerMoleculeB))
    {
      restoreFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);
    const std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(),
                                                        selectedIntegerMoleculeB.end());

    auto restoreB = [&]()
    {
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      systemB.cellListAtomsMoved(fractionalMoleculeB);
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      systemB.cellListAtomsMoved(selectedIntegerMoleculeB);
    };

    const CBMC::GrowContext retraceContextB = systemB.makeGrowContext();
    CBMC::RetraceResult retraceDataB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&] { return CBMC::retraceMolecule(random, retraceContextB, componentB, selectedIntegerMoleculeB); });

    RunningEnergy EwaldFourierDifferenceRetraceB = timed(
        systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
        [&]
        {
          return Interactions::energyDifferenceEwaldFourier(
              systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
              systemB.forceField, systemB.simulationBox, systemB.components, {}, selectedIntegerMoleculeB, netChargeB);
        });
    netChargeB += Interactions::scaledNetChargeDifference({}, selectedIntegerMoleculeB);

    RunningEnergy tailEnergyDifferenceRetraceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&] { return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, {}, selectedIntegerMoleculeB); });

    // System B, step 2: the retraced molecule leaves the background for the grow of the new fractional molecule.
    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingOff();
    }
    systemB.cellListAtomsMoved(selectedIntegerMoleculeB);

    // System B, step 3: CBMC grow the new fractional molecule at lambda into the fractional slot (the slot's own
    // placeholder is excluded by molecule id).
    const std::size_t globalFractionalMoleculeIndexB =
        systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB);
    const CBMC::GrowContext growContextB = systemB.makeGrowContext();
    std::optional<CBMC::GrowResult> growDataB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return CBMC::growNewMolecule(random, growContextB, componentB,
                                             {.componentId = selectedComponent,
                                              .moleculeId = globalFractionalMoleculeIndexB,
                                              .scaling = oldLambda,
                                              .groupId = componentB.lambdaGC.dUdlambdaGroupId,
                                              .isFractional = true});
              });

    if (!growDataB)
    {
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();
      return std::nullopt;
    }

    std::copy(growDataB->atoms.begin(), growDataB->atoms.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<std::uint32_t>(globalFractionalMoleculeIndexB);
    }
    systemB.cellListAtomsMoved(fractionalMoleculeB);

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation leaves B; the grown
    // conformation enters B through the grow energies. Not part of the acceptance rule.
    RunningEnergy intraEnergyPlaceholderB = componentB.intraMolecularPotentials.computeInternalEnergies(
        systemB.forceField, systemB.simulationBox, oldFractionalMoleculeB);

    RunningEnergy EwaldFourierDifferenceGrowB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.trialEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, fractionalMoleculeB,
                    oldFractionalMoleculeB, netChargeB);
              });
    netChargeB += Interactions::scaledNetChargeDifference(fractionalMoleculeB, oldFractionalMoleculeB);

    RunningEnergy tailEnergyDifferenceGrowB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB,
                                                   oldFractionalMoleculeB);
              });

    // System B, step 4: the fractional molecules of the other components switch on (unbiased), with the new
    // fractional molecule present so that their cross terms with it are counted here.
    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return computeFlagSwapOtherComponentsEnergyDifference(systemB, tailEffB, tailGroupB, netChargeB,
                                                                      snapshotB, selectedComponent, true, systemA);
              });
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();
      return std::nullopt;
    }

    // Terms of B that enter the acceptance rule explicitly (everything the retrace and the grow do not see).
    RunningEnergy unbiasedDifferenceB = EwaldFourierDifferenceRetraceB + tailEnergyDifferenceRetraceB +
                                        EwaldFourierDifferenceGrowB + tailEnergyDifferenceGrowB +
                                        flagSwapOtherComponentsDifferenceB.value();

    RunningEnergy energyDifferenceB =
        -(retraceDataB.energies - EwaldFourierDifferenceRetraceB - tailEnergyDifferenceRetraceB) + growDataB->energies +
        EwaldFourierDifferenceGrowB + tailEnergyDifferenceGrowB + flagSwapOtherComponentsDifferenceB.value() -
        intraEnergyPlaceholderB;

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];
    // One integer molecule moves from B to A and the fractional molecule from A to B: the volume factors of the
    // two transfers cancel, leaving the molecule-count factor of the integer transfer.
    double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));

    componentA.mc_moves_statistics.addConstructed(move, 0);

    // Rosenbluth weights through their exact logarithms (raw weights of long chains underflow to zero). The grow and
    // retrace of an integer molecule (A and B) and of the fractional molecule at the same lambda (B and A) pair up,
    // so the ideal-gas reference weights cancel; their energies are in the weights and must not be counted again.
    const double physicalAcceptance =
        std::exp(std::log(preFactor) + growDataA->logRosenbluthWeight - retraceDataA.logRosenbluthWeight +
                 growDataB->logRosenbluthWeight - retraceDataB.logRosenbluthWeight -
                 systemA.beta * (unbiasedDifferenceA.potentialEnergy() + unbiasedDifferenceB.potentialEnergy()));
    if (!tmmcTrial.transferIsInBounds())
    {
      tmmcTrial.recordTransfer(physicalAcceptance);
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();
      return std::nullopt;
    }
    const double tmmcBias = tmmcTrial.biasFactor();
    tmmcTrial.recordTransfer(physicalAcceptance);

    if (random.uniform() < physicalAcceptance * std::exp(biasTerm) * tmmcBias)
    {
      componentA.mc_moves_statistics.addAccepted(move, 0);

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      // System A: the placeholder record of B moves into the fractional slot (its atoms are already in place), the
      // grown integer molecule is appended.
      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
      systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)] =
          systemB.moleculeData[globalFractionalMoleculeIndexB];
      systemA.insertMolecule(selectedComponent, growDataA->molecule, growDataA->atoms);

      // System B: the grown fractional molecule takes the slot (its atoms are already in place), the retraced
      // integer molecule is deleted.
      systemB.moleculeData[globalFractionalMoleculeIndexB] = growDataB->molecule;
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);

      // all fractional molecules take the scaling of their box's new activity state
      const double activeLambda = componentB.lambdaGC.lambdaValue();
      applyFlagSwapScaling(systemA, selectedComponent, false, activeLambda, systemA);
      applyFlagSwapScaling(systemB, selectedComponent, true, activeLambda, systemA);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      systemA.computeTailCorrectionCounts();
      systemB.computeTailCorrectionCounts();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    restoreFractionalMolecules(systemA, snapshotA);
    restoreB();

    return std::nullopt;
  }
  else if (switchValue < 0.5)
  {
    // Shuffle: move the active fractional molecule from A to B. The fractional molecule is CBMC-retraced in A and
    // CBMC-grown in B, both at the current lambda, so the acceptance rule is the ratio of the two Rosenbluth weights
    // times the Boltzmann factor of the terms CBMC does not see (Ewald, tail, other-component flag swaps). The
    // retrace and the grow are evaluated against the same kind of background (integer molecules only; the fractional
    // molecules of the other components are off in an inactive box), which keeps the move symmetric with its reverse.

    componentA.mc_moves_statistics.addTrial(move, 1);

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::size_t indexFractionalMoleculeB =
        systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    const FractionalSnapshot snapshotA = saveFractionalMolecules(systemA);
    const FractionalSnapshot snapshotB = saveFractionalMolecules(systemB);

    // Snapshot effective tail-correction counts for both systems; threaded through the sequential sub-steps.
    std::vector<double> tailEffA = systemA.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupA =
        systemA.fractionalPseudoAtomCountsPerGroup;
    std::vector<double> tailEffB = systemB.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupB =
        systemB.fractionalPseudoAtomCountsPerGroup;
    // Coulomb-scaled net charge of each box, threaded across the chained Ewald sub-steps like the tail counts.
    double netChargeA = systemA.netCharge;
    double netChargeB = systemB.netCharge;

    // System A, step 1: Ewald and tail contributions of removing the fractional molecule at lambda (the first Ewald
    // sub-step of the move in A: storedEik -> trialEik). Its short-range and intra-molecular energies come from the
    // CBMC retrace below.
    RunningEnergy EwaldFourierDifferenceRetraceA = timed(
        systemA, componentA, move, Move::Timing::LambdaShuffleEwald,
        [&]
        {
          return Interactions::energyDifferenceEwaldFourier(
              systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.trialEik,
              systemA.forceField, systemA.simulationBox, systemA.components, {}, oldFractionalMoleculeA, netChargeA);
        });
    netChargeA += Interactions::scaledNetChargeDifference({}, oldFractionalMoleculeA);

    RunningEnergy tailEnergyDifferenceRetraceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleTail,
              [&] { return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, {}, oldFractionalMoleculeA); });

    // System A, step 2: the fractional molecules of the other components switch off (unbiased). Evaluated while the
    // selected fractional molecule is still in memory at lambda, so that its cross terms with them are counted here;
    // the retrace below then sees a background of integer molecules only (as the grow of the reverse move does).
    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return computeFlagSwapOtherComponentsEnergyDifference(systemA, tailEffA, tailGroupA, netChargeA,
                                                                      snapshotA, selectedComponent, false, systemA);
              });
    if (!flagSwapOtherComponentsDifferenceA.has_value())
    {
      return std::nullopt;
    }

    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());

    // swap_ranges copies moleculeId across systems; restore the ids valid in each system so that
    // self-exclusion in the span-based energy differences and the CBMC grow works correctly.
    // System A becomes inactive: its fractional must not contribute dUdlambda (groupId 0).
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.isFractional = true;
      atom.groupId = std::uint8_t{0};
      atom.moleculeId =
          static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.isFractional = true;
      atom.groupId = componentB.lambdaGC.dUdlambdaGroupId;
      atom.moleculeId =
          static_cast<std::uint32_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
    }
    // keep the cell lists (used by the CBMC retrace and grow below) in step with the edited slots
    systemA.cellListAtomsMoved(fractionalMoleculeA);
    systemB.cellListAtomsMoved(fractionalMoleculeB);

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation of B enters A; the old
    // fractional conformation leaves A through the retrace energies. Not part of the acceptance rule.
    RunningEnergy intraEnergyPlaceholderA = componentA.intraMolecularPotentials.computeInternalEnergies(
        systemA.forceField, systemA.simulationBox, fractionalMoleculeA);

    // System A, step 3: all fractional molecules of A take the inactive scaling in memory, and the old fractional
    // molecule is retraced at lambda against that background (its own slot is excluded by molecule id).
    applyFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    const CBMC::GrowContext retraceContextA = systemA.makeGrowContext();
    CBMC::RetraceResult retraceDataA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&] { return CBMC::retraceMolecule(random, retraceContextA, componentA, oldFractionalMoleculeA); });

    // Terms of A that enter the acceptance rule explicitly (everything the retrace does not see).
    RunningEnergy unbiasedDifferenceA =
        EwaldFourierDifferenceRetraceA + tailEnergyDifferenceRetraceA + flagSwapOtherComponentsDifferenceA.value();

    RunningEnergy energyDifferenceA =
        -(retraceDataA.energies - EwaldFourierDifferenceRetraceA - tailEnergyDifferenceRetraceA) +
        flagSwapOtherComponentsDifferenceA.value() + intraEnergyPlaceholderA;

    const std::size_t globalFractionalMoleculeIndexB =
        systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB);
    const CBMC::GrowContext growContextB = systemB.makeGrowContext();
    std::optional<CBMC::GrowResult> growData =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return CBMC::growNewMolecule(random, growContextB, componentB,
                                             {.componentId = selectedComponent,
                                              .moleculeId = globalFractionalMoleculeIndexB,
                                              .scaling = oldLambda,
                                              .groupId = componentB.lambdaGC.dUdlambdaGroupId,
                                              .isFractional = true});
              });

    if (!growData)
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      systemB.cellListAtomsMoved(fractionalMoleculeB);
      return std::nullopt;
    }

    std::copy(growData->atoms.begin(), growData->atoms.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<std::uint32_t>(globalFractionalMoleculeIndexB);
    }
    systemB.cellListAtomsMoved(fractionalMoleculeB);

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation leaves B; the grown
    // conformation enters B through the grow energies. Not part of the acceptance rule.
    RunningEnergy intraEnergyPlaceholderB = componentB.intraMolecularPotentials.computeInternalEnergies(
        systemB.forceField, systemB.simulationBox, oldFractionalMoleculeB);

    // The short-range (framework, inter-molecular) and intra-molecular energies of the grown molecule are in
    // 'growData' (energies and Rosenbluth weight); only the Ewald, tail, and other-component terms remain.
    RunningEnergy EwaldEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, fractionalMoleculeB,
                    oldFractionalMoleculeB, netChargeB);
              });
    netChargeB += Interactions::scaledNetChargeDifference(fractionalMoleculeB, oldFractionalMoleculeB);

    RunningEnergy tailEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB,
                                                   oldFractionalMoleculeB);
              });

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB = computeFlagSwapOtherComponentsEnergyDifference(
        systemB, tailEffB, tailGroupB, netChargeB, snapshotB, selectedComponent, true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      systemB.cellListAtomsMoved(fractionalMoleculeB);
      return std::nullopt;
    }

    // Terms of B that enter the acceptance rule explicitly (everything the grow does not see).
    RunningEnergy unbiasedDifferenceB =
        EwaldEnergyDifferenceB + tailEnergyDifferenceB + flagSwapOtherComponentsDifferenceB.value();

    RunningEnergy energyDifferenceB = growData->energies + EwaldEnergyDifferenceB + tailEnergyDifferenceB +
                                      flagSwapOtherComponentsDifferenceB.value() - intraEnergyPlaceholderB;

    componentA.mc_moves_statistics.addConstructed(move, 1);

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];
    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;

    // Rosenbluth weights through their exact logarithms (raw weights of long chains underflow to zero). The grow in B
    // and the retrace in A are both of the fractional molecule at the same lambda, so the ideal-gas reference weight
    // cancels; their energies are in the weights and must not be counted again here.
    if (random.uniform() <
        std::exp(std::log(preFactor) + growData->logRosenbluthWeight - retraceDataA.logRosenbluthWeight -
                 systemA.beta * (unbiasedDifferenceA.potentialEnergy() + unbiasedDifferenceB.potentialEnergy()) +
                 biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);

      std::swap(systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
                systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)]);
      systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)] =
          growData->molecule;
      std::copy(growData->atoms.begin(), growData->atoms.end(), fractionalMoleculeB.begin());

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      const double activeLambda = componentB.lambdaGC.lambdaValue();
      applyFlagSwapScaling(systemA, selectedComponent, false, activeLambda, systemA);
      applyFlagSwapScaling(systemB, selectedComponent, true, activeLambda, systemA);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      systemA.computeTailCorrectionCounts();
      systemB.computeTailCorrectionCounts();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    restoreFractionalMolecules(systemA, snapshotA);
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
    systemB.cellListAtomsMoved(fractionalMoleculeB);

    return std::nullopt;
  }
  else  // lambda move
  {
    componentA.mc_moves_statistics.addTrial(move, 2);

    if (selectedNewBin < 0) return std::nullopt;
    if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfSamplePoints)) return std::nullopt;

    std::size_t newBin = static_cast<std::size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);

    std::vector<Atom> trialPositions(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::transform(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), trialPositions.begin(),
                   [&](Atom a)
                   {
                     a.setScaling(newLambda);
                     return a;
                   });

    // Snapshot effective tail-correction counts; threaded through the sub-step.
    std::vector<double> tailEffA = systemA.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupA =
        systemA.fractionalPseudoAtomCountsPerGroup;

    std::optional<RunningEnergy> frameworkEnergyDifference =
        timed(systemA, componentA, move, Move::Timing::LambdaChangeNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
                    systemA.frameworkCellList(), systemA.spanOfFrameworkAtoms(), trialPositions, fractionalMoleculeA);
              });

    if (!frameworkEnergyDifference.has_value()) return std::nullopt;

    std::optional<RunningEnergy> interEnergyDifference =
        timed(systemA, componentA, move, Move::Timing::LambdaChangeNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                           systemA.spanOfMoleculeAtoms(),
                                                                           trialPositions, fractionalMoleculeA);
              });

    if (!interEnergyDifference.has_value()) return std::nullopt;

    RunningEnergy EwaldFourierDifference =
        timed(systemA, componentA, move, Move::Timing::LambdaChangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.trialEik,
                    systemA.forceField, systemA.simulationBox, systemA.components, trialPositions, fractionalMoleculeA,
                    systemA.netCharge);
              });

    RunningEnergy tailEnergyDifference = timed(
        systemA, componentA, move, Move::Timing::LambdaChangeTail, [&]
        { return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, trialPositions, fractionalMoleculeA); });

    RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() +
                                     EwaldFourierDifference + tailEnergyDifference;

    componentA.mc_moves_statistics.addConstructed(move, 2);

    double biasTerm = lambdaA.biasFactor[newBin] - lambdaA.biasFactor[oldBin];

    if (random.uniform() < std::exp(-systemA.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);

      componentA.mc_moves_statistics.addAccepted(move, 2);

      std::copy(trialPositions.begin(), trialPositions.end(), fractionalMoleculeA.begin());
      systemA.cellListAtomsMoved(fractionalMoleculeA);
      systemA.computeTailCorrectionCounts();

      componentA.lambdaGC.setCurrentBin(newBin);

      return std::make_pair(energyDifference, RunningEnergy());
    };

    return std::nullopt;
  }

  return std::nullopt;
}
