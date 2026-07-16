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
import cbmc_chain_data;
import cbmc_interactions;
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

// All systems have a fractional molecule, only one of these is 'active', the others are switched off with 'lambda=0'.
// Serial CFCMC: at most one simulation box has containsTheFractionalMolecule==true; systemA must be that box.
// Integer insertions and deletions use CBMC grow and retrace.

using GibbsSwapFractionalSnapshot = std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::vector<Atom>>>;

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

GibbsSwapFractionalSnapshot saveGibbsSwapFractionalMolecules(const System& system)
{
  GibbsSwapFractionalSnapshot snapshot;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }
    const std::size_t index = system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, componentId);
    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    snapshot.emplace_back(std::make_pair(componentId, index),
                          std::vector<Atom>(fractionalMolecule.begin(), fractionalMolecule.end()));
  }
  return snapshot;
}

void restoreGibbsSwapFractionalMolecules(System& system, const GibbsSwapFractionalSnapshot& snapshot)
{
  for (const auto& [indices, moleculeAtoms] : snapshot)
  {
    const auto& [componentId, index] = indices;
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    std::copy(moleculeAtoms.begin(), moleculeAtoms.end(), fractionalMolecule.begin());
  }
}

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

std::optional<RunningEnergy> computeGibbsSwapFractionalMoleculesEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    std::span<const Atom> newAtomsCombined, std::span<const Atom> oldAtomsCombined) noexcept
{
  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtomsCombined, oldAtomsCombined);
  if (!frameworkDifference.has_value())
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newAtomsCombined, oldAtomsCombined);
  if (!moleculeDifference.has_value())
  {
    return std::nullopt;
  }

  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newAtomsCombined, oldAtomsCombined, system.netCharge);

  RunningEnergy tailDifference =
      computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, newAtomsCombined, oldAtomsCombined);

  return frameworkDifference.value() + moleculeDifference.value() + ewaldDifference + tailDifference;
}

std::optional<RunningEnergy> computeSerialFlagSwapFractionalEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    const GibbsSwapFractionalSnapshot& beforeSnapshot, std::size_t selectedComponent, bool systemBecomesActive,
    double selectedActiveLambda, const System& lambdaSourceSystem) noexcept
{
  std::vector<Atom> oldAtomsCombined;
  for (const auto& [indices, moleculeAtoms] : beforeSnapshot)
  {
    (void)indices;
    oldAtomsCombined.insert(oldAtomsCombined.end(), moleculeAtoms.begin(), moleculeAtoms.end());
  }

  // trial state: all fractional molecules take the scaling of the target activity state
  std::vector<Atom> newAtomsCombined;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t indexFractionalMolecule =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, componentId);
    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, indexFractionalMolecule);
    for (Atom atom : fractionalMolecule)
    {
      if (systemBecomesActive)
      {
        const double lambda = (componentId == selectedComponent)
                                  ? selectedActiveLambda
                                  : lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
        atom.setScalingToFractional(lambda, system.components[componentId].lambdaGC.dUdlambdaGroupId);
      }
      else
      {
        // inactive fractional molecules must not contribute dUdlambda
        atom.setScalingToInactiveFractional();
      }
      newAtomsCombined.push_back(atom);
    }
  }

  return computeGibbsSwapFractionalMoleculesEnergyDifference(system, tailEffectiveCounts, tailGroupCounts,
                                                             newAtomsCombined, oldAtomsCombined);
}

std::optional<RunningEnergy> computeSerialFlagSwapOtherComponentsEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    const GibbsSwapFractionalSnapshot& beforeSnapshot, std::size_t selectedComponent, bool systemBecomesActive,
    const System& lambdaSourceSystem) noexcept
{
  std::vector<Atom> oldAtomsCombined;
  std::vector<Atom> newAtomsCombined;

  for (const auto& [indices, moleculeAtoms] : beforeSnapshot)
  {
    const auto& [componentId, index] = indices;
    (void)index;
    if (componentId == selectedComponent)
    {
      continue;
    }
    oldAtomsCombined.insert(oldAtomsCombined.end(), moleculeAtoms.begin(), moleculeAtoms.end());

    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    for (Atom atom : fractionalMolecule)
    {
      if (systemBecomesActive)
      {
        const double lambda = lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
        atom.setScalingToFractional(lambda, system.components[componentId].lambdaGC.dUdlambdaGroupId);
      }
      else
      {
        // inactive fractional molecules must not contribute dUdlambda
        atom.setScalingToInactiveFractional();
      }
      newAtomsCombined.push_back(atom);
    }
  }

  if (oldAtomsCombined.empty())
  {
    return RunningEnergy{};
  }

  return computeGibbsSwapFractionalMoleculesEnergyDifference(system, tailEffectiveCounts, tailGroupCounts,
                                                             newAtomsCombined, oldAtomsCombined);
}

void applySerialGibbsFlagSwapScaling(System& system, std::size_t selectedComponent, bool systemBecomesActive,
                                     double selectedActiveLambda, const System& lambdaSourceSystem)
{
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t indexFractionalMolecule =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, componentId);
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(componentId, indexFractionalMolecule);
    for (Atom& atom : fractionalMolecule)
    {
      if (systemBecomesActive)
      {
        const double lambda = (componentId == selectedComponent)
                                  ? selectedActiveLambda
                                  : lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
        atom.setScalingToFractional(lambda, system.components[componentId].lambdaGC.dUdlambdaGroupId);
      }
      else
      {
        // inactive fractional molecules must not contribute dUdlambda
        atom.setScalingToInactiveFractional();
      }
    }
  }
}

void syncOtherComponentLambdaBins(System& systemBecomingInactive, System& systemBecomingActive,
                                  std::size_t selectedComponent)
{
  for (std::size_t componentId = 0; componentId < systemBecomingInactive.components.size(); ++componentId)
  {
    if (componentId == selectedComponent ||
        systemBecomingInactive.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    systemBecomingActive.components[componentId].lambdaGC.setCurrentBin(
        systemBecomingInactive.components[componentId].lambdaGC.currentBin);
  }
}

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CBCFCMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent,
    [[maybe_unused]] std::size_t& fractionalMoleculeSystem)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
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

  // Determine cutoff distances based on whether dual cutoff is used (only used for the CBMC grow/retrace).
  double cutOffFrameworkVDWA =
      systemA.forceField.useDualCutOff ? systemA.forceField.dualCutOff : systemA.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDWA =
      systemA.forceField.useDualCutOff ? systemA.forceField.dualCutOff : systemA.forceField.cutOffMoleculeVDW;
  double cutOffCoulombA =
      systemA.forceField.useDualCutOff ? systemA.forceField.dualCutOff : systemA.forceField.cutOffCoulomb;
  double cutOffFrameworkVDWB =
      systemB.forceField.useDualCutOff ? systemB.forceField.dualCutOff : systemB.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDWB =
      systemB.forceField.useDualCutOff ? systemB.forceField.dualCutOff : systemB.forceField.cutOffMoleculeVDW;
  double cutOffCoulombB =
      systemB.forceField.useDualCutOff ? systemB.forceField.dualCutOff : systemB.forceField.cutOffCoulomb;
  Component::GrowType growType = componentA.growType;

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
    const GibbsSwapFractionalSnapshot snapshotA = saveGibbsSwapFractionalMolecules(systemA);
    const GibbsSwapFractionalSnapshot snapshotB = saveGibbsSwapFractionalMolecules(systemB);

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
        componentA.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeA);

    RunningEnergy energyDifferenceA = intraFractionalDifferenceA;

    // Combined flag-swap energy for all fractional molecules in A (selected: new position + off; others: off),
    // computed and applied to memory before the grow so the grow sees the post-flag-swap background.
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA = computeSerialFlagSwapFractionalEnergyDifference(
        systemA, tailEffA, tailGroupA, snapshotA, selectedComponent, false, oldLambda, systemA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }
    energyDifferenceA += flagSwapFractionalDifferenceA.value();

    applySerialGibbsFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    // System A: CBMC grow a new integer molecule
    // Trial global molecule id (unique, not yet in the system); insertMolecule assigns the final id.
    std::size_t newMoleculeIndex = systemA.numberOfMolecules();
    const CBMC::GrowContext growContextA{systemA.hasExternalField, systemA.forceField, systemA.simulationBox,
                                         systemA.interpolationGrids, systemA.externalFieldInterpolationGrid,
                                         systemA.framework, systemA.spanOfFrameworkAtoms(),
                                         systemA.spanOfMoleculeAtoms(), systemA.beta, cutOffFrameworkVDWA,
                                         cutOffMoleculeVDWA, cutOffCoulombA};
    time_begin = std::chrono::steady_clock::now();
    std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
        random, growContextA, componentA, selectedComponent, growType, newMoleculeIndex, 1.0, false, false);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);

    if (!growData)
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    if (systemA.insideBlockedPockets(componentA, growData->atoms))
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    if (systemA.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
      // cut-offs, so that Rosenbluth weight and energies behave as if grown at the full cut-offs.
      std::optional<RunningEnergy> correctionNew =
          CBMC::computeDualCutOffCorrection(growContextA, componentA, growData->atoms);
      if (!correctionNew.has_value())
      {
        restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
        return std::nullopt;
      }

      growData->energies += correctionNew.value();
      growData->RosenbluthWeight *= std::exp(-systemA.beta * correctionNew->potentialEnergy());
    }

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldFourierDifferenceGrowA = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.totalEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, std::span(growData->atoms.begin(), growData->atoms.end()), {},
        systemA.netCharge);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceGrowA = computeTailEnergyDifference(
        systemA, tailEffA, tailGroupA, std::span(growData->atoms.begin(), growData->atoms.end()), {});
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeTail] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeTail] += (time_end - time_begin);

    energyDifferenceA += growData->energies + EwaldFourierDifferenceGrowA + tailEnergyDifferenceGrowA;

    double correctionFactorEwaldGrowA = std::exp(
        -systemA.beta * (EwaldFourierDifferenceGrowA.potentialEnergy() + tailEnergyDifferenceGrowA.potentialEnergy()));

    // System B: CBMC retrace a randomly selected integer molecule
    std::size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);
    std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

    const CBMC::GrowContext retraceContextB{systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
                                            systemB.interpolationGrids, systemB.externalFieldInterpolationGrid,
                                            systemB.framework, systemB.spanOfFrameworkAtoms(),
                                            systemB.spanOfMoleculeAtoms(), systemB.beta, cutOffFrameworkVDWB,
                                            cutOffMoleculeVDWB, cutOffCoulombB};
    time_begin = std::chrono::steady_clock::now();
    ChainRetraceData retraceData =
        CBMC::retraceMoleculeSwapDeletion(random, retraceContextB, componentB, growType, selectedIntegerMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);

    if (systemB.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
      // cut-offs, so that Rosenbluth weight and energies behave as if retraced at the full cut-offs.
      std::optional<RunningEnergy> correctionOld =
          CBMC::computeDualCutOffCorrection(retraceContextB, componentB, oldSelectedIntegerMoleculeB);
      if (!correctionOld.has_value())
      {
        restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
        return std::nullopt;
      }

      retraceData.energies += correctionOld.value();
      retraceData.RosenbluthWeight *= std::exp(-systemB.beta * correctionOld->potentialEnergy());
    }

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldFourierDifferenceRetraceB = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, {}, selectedIntegerMoleculeB, systemB.netCharge);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceRetraceB =
        computeTailEnergyDifference(systemB, tailEffB, tailGroupB, {}, selectedIntegerMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeTail] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeTail] += (time_end - time_begin);

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
          static_cast<std::uint16_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
      atom.setScalingToFractional(oldLambda, componentB.lambdaGC.dUdlambdaGroupId);
    }

    RunningEnergy intraFractionalDifferenceB =
        componentB.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeB);

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB2 = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    if (!frameworkDifferenceB2.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB2 = Interactions::computeInterMolecularEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(), fractionalMoleculeB,
        oldFractionalMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeNonEwald] += (time_end - time_begin);
    if (!moleculeDifferenceB2.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldFourierDifferenceB2 = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.totalEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, systemB.netCharge);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceB2 =
        computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeTail] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaInterchangeTail] += (time_end - time_begin);

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        computeSerialFlagSwapOtherComponentsEnergyDifference(systemB, tailEffB, tailGroupB, snapshotB, selectedComponent,
                                                             true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
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

    const double physicalAcceptance =
        preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) / retraceData.RosenbluthWeight *
        correctionFactorEwaldGrowA * correctionFactorEwaldRetraceB *
        std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()));
    if (!tmmcTrial.transferIsInBounds())
    {
      tmmcTrial.recordTransfer(physicalAcceptance);
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
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
            static_cast<std::uint16_t>(systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
        atom.isFractional = true;
      }

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
      systemA.insertMolecule(selectedComponent, growData->molecule, growData->atoms);
      systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)] =
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)];

      std::swap(
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)],
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexSelectedIntegerMoleculeB)]);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);
      systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);

      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
      for (Atom& atom : fractionalMoleculeB)
      {
        atom.moleculeId =
            static_cast<std::uint16_t>(systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
        atom.isFractional = true;
      }

      const double activeLambda = componentB.lambdaGC.lambdaValue();
      applySerialGibbsFlagSwapScaling(systemA, selectedComponent, false, activeLambda, systemA);
      applySerialGibbsFlagSwapScaling(systemB, selectedComponent, true, activeLambda, systemA);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      systemA.computeTailCorrectionCounts();
      systemB.computeTailCorrectionCounts();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
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
    const GibbsSwapFractionalSnapshot snapshotA = saveGibbsSwapFractionalMolecules(systemA);
    const GibbsSwapFractionalSnapshot snapshotB = saveGibbsSwapFractionalMolecules(systemB);

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
        componentA.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeA);

    // Combined flag-swap energy for all fractional molecules in A (selected: new position + off; others: off),
    // computed in one batch to count the selected-other fractional cross-terms exactly once, and applied to
    // memory so the remainder of the move sees the post-flag-swap background.
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA = computeSerialFlagSwapFractionalEnergyDifference(
        systemA, tailEffA, tailGroupA, snapshotA, selectedComponent, false, oldLambda, systemA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceA = intraEnergyDifferenceA + flagSwapFractionalDifferenceA.value();

    applySerialGibbsFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    const std::size_t globalFractionalMoleculeIndexB =
        systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB);
    const CBMC::GrowContext growContextB{systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
                                         systemB.interpolationGrids, systemB.externalFieldInterpolationGrid,
                                         systemB.framework, systemB.spanOfFrameworkAtoms(),
                                         systemB.spanOfMoleculeAtoms(), systemB.beta, cutOffFrameworkVDWB,
                                         cutOffMoleculeVDWB, cutOffCoulombB};
    time_begin = std::chrono::steady_clock::now();
    std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
        random, growContextB, componentB, selectedComponent, growType, globalFractionalMoleculeIndexB, oldLambda,
        componentB.lambdaGC.dUdlambdaGroupId, true);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);

    if (!growData)
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    if (systemB.insideBlockedPockets(componentB, growData->atoms))
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    if (systemB.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
      // cut-offs, so that Rosenbluth weight and energies behave as if grown at the full cut-offs.
      std::optional<RunningEnergy> correctionNew =
          CBMC::computeDualCutOffCorrection(growContextB, componentB, growData->atoms);
      if (!correctionNew.has_value())
      {
        restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
        std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
        return std::nullopt;
      }

      growData->energies += correctionNew.value();
      growData->RosenbluthWeight *= std::exp(-systemB.beta * correctionNew->potentialEnergy());
    }

    std::copy(growData->atoms.begin(), growData->atoms.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<std::uint32_t>(globalFractionalMoleculeIndexB);
    }

    RunningEnergy intraEnergyDifference =
        componentB.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeB);

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    if (!frameworkDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB = Interactions::computeInterMolecularEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(), fractionalMoleculeB,
        oldFractionalMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaShuffleNonEwald] += (time_end - time_begin);
    if (!moleculeDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldEnergyDifferenceB = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, systemB.netCharge);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaShuffleEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaShuffleEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceB =
        computeTailEnergyDifference(systemB, tailEffB, tailGroupB, fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaShuffleTail] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaShuffleTail] += (time_end - time_begin);

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        computeSerialFlagSwapOtherComponentsEnergyDifference(systemB, tailEffB, tailGroupB, snapshotB, selectedComponent,
                                                             true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
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

    if (random.uniform() <
        preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) *
            std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()) +
                     biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

      std::swap(systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
                systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)]);
      systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)] =
          growData->molecule;
      std::copy(growData->atoms.begin(), growData->atoms.end(), fractionalMoleculeB.begin());

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      const double activeLambda = componentB.lambdaGC.lambdaValue();
      applySerialGibbsFlagSwapScaling(systemA, selectedComponent, false, activeLambda, systemA);
      applySerialGibbsFlagSwapScaling(systemB, selectedComponent, true, activeLambda, systemA);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      systemA.computeTailCorrectionCounts();
      systemB.computeTailCorrectionCounts();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
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

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
        systemA.spanOfFrameworkAtoms(), trialPositions, fractionalMoleculeA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaChangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaChangeNonEwald] += (time_end - time_begin);

    if (!frameworkEnergyDifference.has_value()) return std::nullopt;

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = Interactions::computeInterMolecularEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(), trialPositions, fractionalMoleculeA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaChangeNonEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaChangeNonEwald] += (time_end - time_begin);

    if (!interEnergyDifference.has_value()) return std::nullopt;

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldFourierDifference = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, trialPositions, fractionalMoleculeA, systemA.netCharge);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaChangeEwald] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaChangeEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifference =
        computeTailEnergyDifference(systemA, tailEffA, tailGroupA, trialPositions, fractionalMoleculeA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaChangeTail] += (time_end - time_begin);
    systemA.mc_moves_cputime[move][Move::Timing::LambdaChangeTail] += (time_end - time_begin);

    RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() +
                                     EwaldFourierDifference + tailEnergyDifference;

    componentA.mc_moves_statistics.addConstructed(move, 2);

    double biasTerm = lambdaA.biasFactor[newBin] - lambdaA.biasFactor[oldBin];

    if (random.uniform() < std::exp(-systemA.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);

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
