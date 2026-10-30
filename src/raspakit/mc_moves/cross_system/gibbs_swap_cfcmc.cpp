module;

module mc_moves_gibbs_swap_cfcmc;

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
import running_energy;
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
import intra_molecular_potentials;

// All systems have a fractional molecule, only one of these is 'active', the others are switched off with 'lambda=0'.
// Serial CFCMC: at most one simulation box has containsTheFractionalMolecule==true; systemA must be that box. The flag
// is a property of the box, so when a move transfers the activity, the fractional molecules of the other swappable
// components follow (see mc_moves_gibbs_swap_serial_common). Integer insertions use random placement.
// Implementation advantage: the number of fractional molecules per system remains constant.

using namespace MC_Moves::GibbsSwapSerialCommon;

namespace
{
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

// systemA contains the fractional molecule
std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CFCMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent,
    [[maybe_unused]] std::size_t& fractionalMoleculeSystem)
{
  Move::Types move = Move::Types::GibbsSwapCFCMC;
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
    // Swap move:
    // Changing the fractional molecule into a whole molecule, keeping its position fixed
    // Changing a randomly selected molecule in the other simulation box into a fractional molecule (at same lambda)

    componentA.mc_moves_statistics.addTrial(move, 0);

    if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, selectedComponent);
    std::size_t indexFractionalMoleculeB =
        systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    // make copy of old fractional molecules for reference and restoring
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

    // System A: the fractional slot takes the inactive placeholder of B; all fractional molecules of A switch off
    //=============================================================================================================

    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.moleculeId =
          static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation enters A, the old
    // fractional conformation stays in A as the new integer molecule.
    RunningEnergy intraFractionalDifferenceA = componentA.intraMolecularPotentials.computeInternalEnergies(
        systemA.forceField, systemA.simulationBox, fractionalMoleculeA);

    // Combined flag-swap energy for all fractional molecules in A (selected: placeholder + off; others: off),
    // computed in one batch to count the cross-terms exactly once, and applied to memory so the remainder of the
    // move sees the post-flag-swap background.
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return computeFlagSwapFractionalEnergyDifference(systemA, tailEffA, tailGroupA, netChargeA, snapshotA,
                                                                 selectedComponent, false, oldLambda, systemA);
              });
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    applyFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    // System A, step 2: the old fractional molecule becomes an integer molecule at its position
    std::vector<Atom> newMolecule(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
    for (Atom& atom : newMolecule)
    {
      atom.setScalingToInteger();
      atom.moleculeId = static_cast<std::uint32_t>(systemA.numberOfMolecules());
    }

    std::optional<RunningEnergy> frameworkDifferenceA2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
                    systemA.frameworkCellList(), systemA.spanOfFrameworkAtoms(), newMolecule, {});
              });

    if (!frameworkDifferenceA2.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    std::optional<RunningEnergy> moleculeDifferenceA2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(
                    systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(), newMolecule, {});
              });

    if (!moleculeDifferenceA2.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceA2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.trialEik, systemA.trialEik,
                    systemA.forceField, systemA.simulationBox, systemA.components, newMolecule, {}, netChargeA);
              });
    netChargeA += Interactions::scaledNetChargeDifference(newMolecule, {});

    RunningEnergy tailEnergyDifferenceA2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&] { return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, newMolecule, {}); });

    RunningEnergy energyDifferenceA = intraFractionalDifferenceA + flagSwapFractionalDifferenceA.value() +
                                      frameworkDifferenceA2.value() + moleculeDifferenceA2.value() +
                                      EwaldFourierDifferenceA2 + tailEnergyDifferenceA2;

    // System B: Changing a randomly selected molecule in the other simulation box into a fractional molecule
    // (at the same lambda); the fractional molecules of the other components switch on
    //================================================================================================================

    std::size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);

    // make copy of selected molecule for reference and restoring
    const std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(),
                                                        selectedIntegerMoleculeB.end());

    auto restoreB = [&]()
    {
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
    };

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingOff();
    }

    std::optional<RunningEnergy> frameworkDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
                    systemB.frameworkCellList(), systemB.spanOfFrameworkAtoms(), selectedIntegerMoleculeB,
                    oldSelectedIntegerMoleculeB);
              });
    if (!frameworkDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();
      return std::nullopt;
    }

    std::optional<RunningEnergy> moleculeDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(
                    systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(), selectedIntegerMoleculeB,
                    oldSelectedIntegerMoleculeB);
              });
    if (!moleculeDifferenceB.has_value())
    {
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, selectedIntegerMoleculeB,
                    oldSelectedIntegerMoleculeB, netChargeB);
              });
    netChargeB += Interactions::scaledNetChargeDifference(selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);

    RunningEnergy tailEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, selectedIntegerMoleculeB,
                                                   oldSelectedIntegerMoleculeB);
              });

    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId =
          static_cast<std::uint32_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
      atom.setScalingToFractional(oldLambda, componentB.lambdaGC.dUdlambdaGroupId);
    }

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingOff();
      atom.position = systemB.simulationBox.randomPosition(random);
    }

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation leaves B, the selected
    // conformation stays in B as the new fractional molecule.
    RunningEnergy intraFractionalDifferenceB = -componentB.intraMolecularPotentials.computeInternalEnergies(
        systemB.forceField, systemB.simulationBox, oldFractionalMoleculeB);

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
      restoreB();
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
      restoreB();
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.trialEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, systemB.components, fractionalMoleculeB,
                    oldFractionalMoleculeB, netChargeB);
              });
    netChargeB += Interactions::scaledNetChargeDifference(fractionalMoleculeB, oldFractionalMoleculeB);

    RunningEnergy tailEnergyDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
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
      restoreB();
      return std::nullopt;
    }

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() +
                                      EwaldFourierDifferenceB + tailEnergyDifferenceB + frameworkDifferenceB2.value() +
                                      moleculeDifferenceB2.value() + EwaldFourierDifferenceB2 + tailEnergyDifferenceB2 +
                                      intraFractionalDifferenceB + flagSwapOtherComponentsDifferenceB.value();

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));

    componentA.mc_moves_statistics.addConstructed(move, 0);

    // The intra-molecular terms are bookkeeping only (conformations are neither created nor destroyed here).
    const double physicalAcceptance =
        preFactor * std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy() -
                                              intraFractionalDifferenceA.potentialEnergy() -
                                              intraFractionalDifferenceB.potentialEnergy()));
    if (!tmmcTrial.transferIsInBounds())
    {
      tmmcTrial.recordTransfer(physicalAcceptance);
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();
      return std::nullopt;
    }
    const double tmmcBias = tmmcTrial.biasFactor();
    tmmcTrial.recordTransfer(physicalAcceptance);

    // apply acceptance/rejection rule
    if (random.uniform() < physicalAcceptance * std::exp(biasTerm) * tmmcBias)
    {
      componentA.mc_moves_statistics.addAccepted(move, 0);

      // restore
      restoreFractionalMolecules(systemA, snapshotA);
      restoreB();

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      // copy the fractional molecule from B to A
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
      for (Atom& atom : fractionalMoleculeA)
      {
        atom.moleculeId =
            static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
        atom.isFractional = true;
      }

      // make old fractional molecule integer
      std::vector<Atom> addedMolecule = std::vector<Atom>(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
      for (Atom& atom : addedMolecule)
      {
        atom.setScalingToInteger();
      }
      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
      systemA.insertMolecule(
          selectedComponent,
          systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
          addedMolecule);
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
    // Move fractional molecule to the other box (from A to a random position in B)

    componentA.mc_moves_statistics.addTrial(move, 1);

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, selectedComponent);
    std::size_t indexFractionalMoleculeB =
        systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    // make copy of old fractional molecule for reference and restoring
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

    // swap the active and the inactive fractional molecule
    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());

    // swap_ranges copies moleculeId across systems; restore the ids valid in each system so that the
    // self-exclusion in the span-based energy differences works correctly (the fractional molecule has a
    // different global index in A and B as soon as the preceding components hold different numbers of
    // molecules). System A becomes inactive: its fractional must not contribute dUdlambda (groupId 0).
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

    // Intra-molecular bookkeeping (zero for rigid molecules): the placeholder conformation of B replaces the old
    // fractional conformation in A.
    RunningEnergy intraEnergyDifferenceA =
        componentA.intraMolecularPotentials.computeInternalEnergies(systemA.forceField, systemA.simulationBox,
                                                                    fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(systemA.forceField, systemA.simulationBox,
                                                                    oldFractionalMoleculeA);

    // Combined flag-swap energy for all fractional molecules in A (selected: placeholder + off; others: off),
    // computed in one batch to count the cross-terms exactly once, and applied to memory so the remainder of the
    // move sees the post-flag-swap background.
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return computeFlagSwapFractionalEnergyDifference(systemA, tailEffA, tailGroupA, netChargeA, snapshotA,
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

    // System B: the fractional molecule arrives at a random position (with an ideal-gas conformation)
    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        systemB.equilibratedIdealGasMoleculeRandomInBox(random, selectedComponent);

    std::transform(fractionalMoleculeB.begin(), fractionalMoleculeB.end(), trialMolecule.second.begin(),
                   fractionalMoleculeB.begin(),
                   [](const Atom& a, const Atom& b)
                   {
                     return Atom(b.position, a.charge, a.scalingVDW, a.scalingCoulomb, a.moleculeId, a.type,
                                 a.componentId, a.groupId, a.isFractional);
                   });

    // Intra-molecular bookkeeping (zero for rigid molecules): the ideal-gas conformation replaces the placeholder
    // conformation in B. It is drawn from the ideal-gas Boltzmann distribution, so it does not enter the acceptance.
    RunningEnergy intraEnergyDifferenceB =
        componentB.intraMolecularPotentials.computeInternalEnergies(systemB.forceField, systemB.simulationBox,
                                                                    fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(systemB.forceField, systemB.simulationBox,
                                                                    oldFractionalMoleculeB);

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
      return std::nullopt;
    }

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() +
                                      EwaldEnergyDifferenceB + tailEnergyDifferenceB + intraEnergyDifferenceB +
                                      flagSwapOtherComponentsDifferenceB.value();

    componentA.mc_moves_statistics.addConstructed(move, 1);

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;

    // apply acceptance/rejection rule (the intra-molecular terms are bookkeeping only, see above)
    if (random.uniform() < preFactor * std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() +
                                                                 energyDifferenceB.potentialEnergy() -
                                                                 intraEnergyDifferenceA.potentialEnergy() -
                                                                 intraEnergyDifferenceB.potentialEnergy()) +
                                                biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);

      std::swap(systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
                systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)]);
      systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)] =
          trialMolecule.first;

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

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

    // reject, set fractional molecules back to old state
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
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, selectedComponent);
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

    // apply acceptance/rejection rule
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
