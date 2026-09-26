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

// All systems have a fractional molecule, only one of these is 'active', the others are switched off with 'lambda=0'.
// Implementation advantage: the number of fractional molecules per system remains constant.

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

// Brick-CFCMC-style aggregated tail-correction difference. The effective type counts are threaded across the
// sequential sub-steps of the move: each call returns the difference for its (newAtoms, oldAtoms) change and then
// advances the running counts so the next sub-step sees the updated background.
static RunningEnergy computeTailEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms)
{
  RunningEnergy result =
      Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
          system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms);
  return result;
}

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

  // if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfBins))
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

    // assert(fractionalMoleculeA.front().groupId == std::uint8_t{ 1 });
    // assert(fractionalMoleculeB.front().groupId == std::uint8_t{ 1 });

    // make copy of old fractional molecule for reference and restoring
    const std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    const std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    std::vector<Atom> oldFractionalMoleculeB2(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

    // Snapshot effective tail-correction counts for both systems; threaded through the sequential sub-steps.
    std::vector<double> tailEffA = systemA.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupA =
        systemA.fractionalPseudoAtomCountsPerGroup;
    std::vector<double> tailEffB = systemB.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupB =
        systemB.fractionalPseudoAtomCountsPerGroup;

    // System A: Changing the fractional molecule into a whole molecule, keeping its position fixed
    //=============================================================================================

    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.moleculeId =
          static_cast<std::uint32_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }

    std::optional<RunningEnergy> frameworkDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
                    systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
              });

    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    std::optional<RunningEnergy> moleculeDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                           systemA.spanOfMoleculeAtoms(),
                                                                           fractionalMoleculeA, oldFractionalMoleculeA);
              });

    if (!moleculeDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.trialEik,
                    systemA.forceField, systemA.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA,
                    systemA.netCharge);
              });

    RunningEnergy tailEnergyDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, fractionalMoleculeA,
                                                   oldFractionalMoleculeA);
              });

    // step 2

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
                    systemA.spanOfFrameworkAtoms(), newMolecule, {});
              });

    if (!frameworkDifferenceA2.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceA2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.trialEik, systemA.trialEik,
                    systemA.forceField, systemA.simulationBox, newMolecule, {}, systemA.netCharge);
              });

    RunningEnergy tailEnergyDifferenceA2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&] { return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, newMolecule, {}); });

    RunningEnergy energyDifferenceA = frameworkDifferenceA.value() + moleculeDifferenceA.value() +
                                      EwaldFourierDifferenceA + tailEnergyDifferenceA + frameworkDifferenceA2.value() +
                                      moleculeDifferenceA2.value() + EwaldFourierDifferenceA2 + tailEnergyDifferenceA2;

    // System B: Changing a randomly selected molecule in the other simulation box into a fractional molecule
    // (at the same lambda)
    //================================================================================================================

    std::size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);

    // make copy of selected molecule for reference and restoring
    std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());
    std::vector<Atom> oldSelectedIntegerMoleculeB2(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

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
                    systemB.spanOfFrameworkAtoms(), selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
              });
    if (!frameworkDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
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
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldFourierDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB,
                    systemB.netCharge);
              });

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
          static_cast<std::uint16_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
      atom.setScalingToFractional(oldLambda, componentB.lambdaGC.dUdlambdaGroupId);
    }

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingOff();
      atom.position = systemA.simulationBox.randomPosition(random);
    }

    std::optional<RunningEnergy> frameworkDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
                    systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
              });
    if (!frameworkDifferenceB2.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
                    systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB,
                    systemB.netCharge);
              });

    RunningEnergy tailEnergyDifferenceB2 =
        timed(systemA, componentA, move, Move::Timing::LambdaInterchangeTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB,
                                                   oldFractionalMoleculeB);
              });

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() +
                                      EwaldFourierDifferenceB + tailEnergyDifferenceB + frameworkDifferenceB2.value() +
                                      moleculeDifferenceB2.value() + EwaldFourierDifferenceB2 + tailEnergyDifferenceB2;

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));

    componentA.mc_moves_statistics.addConstructed(move, 0);

    const double physicalAcceptance =
        preFactor *
        std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()));
    if (!tmmcTrial.transferIsInBounds())
    {
      tmmcTrial.recordTransfer(physicalAcceptance);
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }
    const double tmmcBias = tmmcTrial.biasFactor();
    tmmcTrial.recordTransfer(physicalAcceptance);

    // apply acceptance/rejection rule
    if (random.uniform() < physicalAcceptance * std::exp(biasTerm) * tmmcBias)
    {
      componentA.mc_moves_statistics.addAccepted(move, 0);

      // restore
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());

      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      // copy the fractional molecule from B to A
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
      for (Atom& atom : fractionalMoleculeA)
      {
        atom.moleculeId =
            static_cast<std::uint16_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
      }

      // make old fractional molecule integer
      std::vector<Atom> addedMolecule = std::vector<Atom>(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
      for (Atom& atom : addedMolecule)
      {
        atom.setScalingToInteger();
      }
      systemA.insertMolecule(
          selectedComponent,
          systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
          addedMolecule);
      systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)] =
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)];

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);

      std::swap(
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)],
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexSelectedIntegerMoleculeB)]);
      systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);

      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
      for (Atom& atom : fractionalMoleculeB)
      {
        atom.moleculeId =
            static_cast<std::uint16_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
        atom.setScalingToFractional(oldLambda, componentB.lambdaGC.dUdlambdaGroupId);
      }

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      systemA.computeTailCorrectionCounts();
      systemB.computeTailCorrectionCounts();

      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

    return std::nullopt;
  }
  // else if (selectedNewBin < 0)
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
    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

    // Snapshot effective tail-correction counts for both systems; threaded through the sequential sub-steps.
    std::vector<double> tailEffA = systemA.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupA =
        systemA.fractionalPseudoAtomCountsPerGroup;
    std::vector<double> tailEffB = systemB.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupB =
        systemB.fractionalPseudoAtomCountsPerGroup;

    // swap the active and the inactive fractional molecule
    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());

    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        systemB.equilibratedIdealGasMoleculeRandomInBox(random, selectedComponent);

    // std::copy(trialMolecule.second.begin(), trialMolecule.second.end(), fractionalMoleculeB.begin());
    std::transform(fractionalMoleculeB.begin(), fractionalMoleculeB.end(), trialMolecule.second.begin(),
                   fractionalMoleculeB.begin(),
                   [](const Atom& a, const Atom& b)
                   {
                     return Atom(b.position, a.charge, a.scalingVDW, a.scalingCoulomb, a.moleculeId, a.type,
                                 a.componentId, a.groupId, a.isFractional);
                   });

    std::optional<RunningEnergy> frameworkDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
                    systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
              });

    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    std::optional<RunningEnergy> moleculeDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                           systemA.spanOfMoleculeAtoms(),
                                                                           fractionalMoleculeA, oldFractionalMoleculeA);
              });

    if (!moleculeDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldEnergyDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.trialEik,
                    systemA.forceField, systemA.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA,
                    systemA.netCharge);
              });

    RunningEnergy tailEnergyDifferenceA =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleTail,
              [&]
              {
                return computeTailEnergyDifference(systemA, tailEffA, tailGroupA, fractionalMoleculeA,
                                                   oldFractionalMoleculeA);
              });

    RunningEnergy energyDifferenceA =
        frameworkDifferenceA.value() + moleculeDifferenceA.value() + EwaldEnergyDifferenceA + tailEnergyDifferenceA;

    std::optional<RunningEnergy> frameworkDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleNonEwald,
              [&]
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
                    systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
              });

    if (!frameworkDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy EwaldEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleEwald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
                    systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB,
                    systemB.netCharge);
              });

    RunningEnergy tailEnergyDifferenceB =
        timed(systemA, componentA, move, Move::Timing::LambdaShuffleTail,
              [&]
              {
                return computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB,
                                                   oldFractionalMoleculeB);
              });

    RunningEnergy energyDifferenceB =
        frameworkDifferenceB.value() + moleculeDifferenceB.value() + EwaldEnergyDifferenceB + tailEnergyDifferenceB;

    componentA.mc_moves_statistics.addConstructed(move, 1);

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;

    // apply acceptance/rejection rule
    if (random.uniform() < preFactor * std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() +
                                                                 energyDifferenceB.potentialEnergy()) +
                                                biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);

      std::swap(systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
                systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)]);
      systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)] =
          trialMolecule.first;

      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      systemA.computeTailCorrectionCounts();
      systemB.computeTailCorrectionCounts();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    // reject, set fractional molecule back to old state
    std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
                    systemA.spanOfFrameworkAtoms(), trialPositions, fractionalMoleculeA);
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
                    systemA.forceField, systemA.simulationBox, trialPositions, fractionalMoleculeA, systemA.netCharge);
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
