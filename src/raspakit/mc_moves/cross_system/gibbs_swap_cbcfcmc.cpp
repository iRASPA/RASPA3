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
// Integer insertions and deletions use CBMC grow and retrace.

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
    // Lambda interchange: fractional molecule exchange between boxes with CBMC grow (A) and retrace (B).

    componentA.mc_moves_statistics.addTrial(move, 0);

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

    // System A: copy inactive fractional molecule from B into the active fractional slot
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.moleculeId =
          static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }

    RunningEnergy intraFractionalDifferenceA =
        componentA.intraMolecularPotentials.computeInternalEnergies(systemA.forceField, systemA.simulationBox, fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(systemA.forceField, systemA.simulationBox, oldFractionalMoleculeA);

    RunningEnergy energyDifferenceA = intraFractionalDifferenceA;

    // Combined flag-swap energy for all fractional molecules in A (selected: new position + off; others: off),
    // computed and applied to memory before the grow so the grow sees the post-flag-swap background.
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return computeFlagSwapFractionalEnergyDifference(systemA, tailEffA, tailGroupA, snapshotA,
                                                                       selectedComponent, false, oldLambda, systemA);
              });
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }
    energyDifferenceA += flagSwapFractionalDifferenceA.value();

    applyFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    // System A: CBMC grow a new integer molecule
    // Trial global molecule id (unique, not yet in the system); insertMolecule assigns the final id.
    std::size_t newMoleculeIndex = systemA.numberOfMolecules();
    const CBMC::GrowContext growContextA = systemA.makeGrowContext();
    std::optional<CBMC::GrowResult> growData =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return CBMC::growNewMolecule(random, growContextA, componentA,
                                             {.componentId = selectedComponent, .moleculeId = newMoleculeIndex});
              });

    if (!growData)
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
                    std::span(growData->atoms.begin(), growData->atoms.end()), {}, systemA.netCharge);
              });

    RunningEnergy tailEnergyDifferenceGrowA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemA, tailEffA, tailGroupA,
                                                   std::span(growData->atoms.begin(), growData->atoms.end()), {});
              });

    energyDifferenceA += growData->energies + EwaldFourierDifferenceGrowA + tailEnergyDifferenceGrowA;

    double correctionFactorEwaldGrowA = std::exp(
        -systemA.beta * (EwaldFourierDifferenceGrowA.potentialEnergy() + tailEnergyDifferenceGrowA.potentialEnergy()));

    // System B: CBMC retrace a randomly selected integer molecule
    std::size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);
    std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

    const CBMC::GrowContext retraceContextB = systemB.makeGrowContext();
    CBMC::RetraceResult retraceData =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&] { return CBMC::retraceMolecule(random, retraceContextB, componentB, selectedIntegerMoleculeB); });
    RunningEnergy EwaldFourierDifferenceRetraceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, {}, selectedIntegerMoleculeB, systemB.netCharge);
              });

    RunningEnergy tailEnergyDifferenceRetraceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&] { return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, {}, selectedIntegerMoleculeB); });

    double correctionFactorEwaldRetraceB = std::exp(systemB.beta * (EwaldFourierDifferenceRetraceB.potentialEnergy() +
                                                                    tailEnergyDifferenceRetraceB.potentialEnergy()));

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingOff();
      atom.position = systemB.simulationBox.randomPosition(random);
    }

    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId =
          static_cast<std::uint32_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
      atom.setScalingToFractional(oldLambda, componentB.lambdaGC.dUdlambdaGroupId);
    }

    RunningEnergy intraFractionalDifferenceB =
        componentB.intraMolecularPotentials.computeInternalEnergies(systemB.forceField, systemB.simulationBox, fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(systemB.forceField, systemB.simulationBox, oldFractionalMoleculeB);

    std::optional<RunningEnergy> frameworkDifferenceB2 = timed(
        systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
        [&]
        {
          return Interactions::computeFrameworkMoleculeEnergyDifference(
              systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
              systemB.frameworkCellList(), systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
        });
    if (!frameworkDifferenceB2.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    std::optional<RunningEnergy> moleculeDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                           systemB.spanOfMoleculeAtoms(),
                                                                           fractionalMoleculeB, oldFractionalMoleculeB);
              });
    if (!moleculeDifferenceB2.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.trialEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, fractionalMoleculeB, oldFractionalMoleculeB,
                    systemB.netCharge);
              });

    RunningEnergy tailEnergyDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB,
                                                   oldFractionalMoleculeB);
              });

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        computeFlagSwapOtherComponentsEnergyDifference(systemB, tailEffB, tailGroupB, snapshotB, selectedComponent,
                                                             true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceB =
        -(retraceData.energies - EwaldFourierDifferenceRetraceB - tailEnergyDifferenceRetraceB) +
        frameworkDifferenceB2.value() + moleculeDifferenceB2.value() + EwaldFourierDifferenceB2 +
        tailEnergyDifferenceB2 + intraFractionalDifferenceB + flagSwapOtherComponentsDifferenceB.value();

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];
    double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));
    double idealGasRosenbluthWeight = componentA.idealGasRosenbluthWeight.value_or(1.0);

    componentA.mc_moves_statistics.addConstructed(move, 0);

    // Rosenbluth weights through their exact logarithms (raw weights of long chains underflow to zero).
    const double physicalAcceptance =
        std::exp(std::log(preFactor) + growData->logRosenbluthWeight - std::log(idealGasRosenbluthWeight) -
                 retraceData.logRosenbluthWeight + std::log(correctionFactorEwaldGrowA) +
                 std::log(correctionFactorEwaldRetraceB) -
                 systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()));
    if (!tmmcTrial.transferIsInBounds())
    {
      tmmcTrial.recordTransfer(physicalAcceptance);
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }
    const double tmmcBias = tmmcTrial.biasFactor();
    tmmcTrial.recordTransfer(physicalAcceptance);

    if (random.uniform() < physicalAcceptance * std::exp(biasTerm) * tmmcBias)
    {
      componentA.mc_moves_statistics.addAccepted(move, 0);

      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
      for (Atom& atom : fractionalMoleculeA)
      {
        atom.moleculeId =
            static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
        atom.isFractional = true;
      }

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
      systemA.insertMolecule(selectedComponent, growData->molecule, growData->atoms);
      systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)] =
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)];

      std::swap(
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)],
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexSelectedIntegerMoleculeB)]);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);
      systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);

      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
      for (Atom& atom : fractionalMoleculeB)
      {
        atom.moleculeId =
            static_cast<std::uint32_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
        atom.isFractional = true;
      }

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
    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

    return std::nullopt;
  }
  else if (switchValue < 0.5)
  {
    // Shuffle: move active fractional molecule from A to B using CBMC grow in B.

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

    RunningEnergy intraEnergyDifferenceA =
        componentA.intraMolecularPotentials.computeInternalEnergies(systemA.forceField, systemA.simulationBox, fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(systemA.forceField, systemA.simulationBox, oldFractionalMoleculeA);

    // Combined flag-swap energy for all fractional molecules in A (selected: new position + off; others: off),
    // computed in one batch to count the selected-other fractional cross-terms exactly once, and applied to
    // memory so the remainder of the move sees the post-flag-swap background.
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return computeFlagSwapFractionalEnergyDifference(systemA, tailEffA, tailGroupA, snapshotA,
                                                                       selectedComponent, false, oldLambda, systemA);
              });
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceA = intraEnergyDifferenceA + flagSwapFractionalDifferenceA.value();

    applyFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

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
      return std::nullopt;
    }

    std::copy(growData->atoms.begin(), growData->atoms.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<std::uint32_t>(globalFractionalMoleculeIndexB);
    }

    RunningEnergy intraEnergyDifference =
        componentB.intraMolecularPotentials.computeInternalEnergies(systemB.forceField, systemB.simulationBox, fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(systemB.forceField, systemB.simulationBox, oldFractionalMoleculeB);

    std::optional<RunningEnergy> frameworkDifferenceB = timed(
        systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
        [&]
        {
          return Interactions::computeFrameworkMoleculeEnergyDifference(
              systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
              systemB.frameworkCellList(), systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
        });
    if (!frameworkDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    std::optional<RunningEnergy> moleculeDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                           systemB.spanOfMoleculeAtoms(),
                                                                           fractionalMoleculeB, oldFractionalMoleculeB);
              });
    if (!moleculeDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, fractionalMoleculeB, oldFractionalMoleculeB,
                    systemB.netCharge);
              });

    RunningEnergy tailEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB,
                                                   oldFractionalMoleculeB);
              });

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        computeFlagSwapOtherComponentsEnergyDifference(systemB, tailEffB, tailGroupB, snapshotB, selectedComponent,
                                                             true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() +
                                      EwaldEnergyDifferenceB + tailEnergyDifferenceB + intraEnergyDifference +
                                      flagSwapOtherComponentsDifferenceB.value();

    componentA.mc_moves_statistics.addConstructed(move, 1);

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];
    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;
    double idealGasRosenbluthWeight = componentB.idealGasRosenbluthWeight.value_or(1.0);

    // Rosenbluth weight through its exact logarithm (raw weights of long chains underflow to zero).
    if (random.uniform() <
        std::exp(std::log(preFactor) + growData->logRosenbluthWeight - std::log(idealGasRosenbluthWeight) -
                 systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()) +
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
                    systemA.forceField, systemA.simulationBox, systemA.components, trialPositions, fractionalMoleculeA, systemA.netCharge);
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
      systemA.computeTailCorrectionCounts();

      componentA.lambdaGC.setCurrentBin(newBin);

      return std::make_pair(energyDifference, RunningEnergy());
    };

    return std::nullopt;
  }

  return std::nullopt;
}
