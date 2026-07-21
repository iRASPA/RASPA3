module;

module mc_moves_group_swap_cfcmc;

import std;

import component;
import molecule;
import atom;
import double3;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
import randomnumbers;
import system;
import property_lambda_probability_histogram;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;
import scaling;

// The group CFCMC move generalizes the ion-pair CFCMC move to K = 1 + #satellites molecules that
// share a single lambda coordinate (the lambdaGroupSwap / lambdaGroupSwapCB histogram of the
// central component). Energy differences are computed sequentially, one molecule at a time: each
// molecule is updated in place before the difference of the next molecule is computed, so that
// every intra-group cross-interaction is counted exactly once.

namespace
{

struct GroupMember
{
  std::size_t componentId;
  std::size_t fractionalMoleculeIndex;  // fractional slot index within the component
};

// central molecule first, then the satellites in the order listed in GroupComponents; a component
// occurring multiple times receives consecutive sub-slots
std::vector<GroupMember> buildGroupMembers(const System& system, std::size_t selectedComponent, Move::Types move)
{
  const Component& central = system.components[selectedComponent];
  std::vector<GroupMember> members;
  members.reserve(1 + central.groupComponentIds.size());
  std::vector<std::size_t> occurrence(system.components.size(), 0);
  members.push_back(
      {selectedComponent, system.indexOfFractionalMoleculeForMove(move, selectedComponent, occurrence[selectedComponent])});
  ++occurrence[selectedComponent];
  for (std::size_t satelliteComponentId : central.groupComponentIds)
  {
    members.push_back({satelliteComponentId, system.indexOfFractionalMoleculeForMove(
                                                 move, satelliteComponentId, occurrence[satelliteComponentId])});
    ++occurrence[satelliteComponentId];
  }
  return members;
}

// Accumulate the change of the inter-molecular electric field on every atom when 'oldAtoms' is replaced by
// 'newAtoms' at fixed position (typically a charge-scaling / lambda change).
void accumulateInterMolecularFieldDelta(System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms,
                                        std::span<double3> electricFieldNeighborDelta)
{
  std::vector<double3> fieldNew(newAtoms.size());
  std::vector<double3> fieldOld(oldAtoms.size());
  [[maybe_unused]] std::optional<RunningEnergy> energy =
      Interactions::computeInterMolecularPolarizationElectricFieldDifference(
          system.forceField, system.simulationBox, electricFieldNeighborDelta, fieldNew, fieldOld,
          system.spanOfMoleculeAtoms(), newAtoms, oldAtoms);
}

// difference in scaled net charge when 'oldAtoms' is replaced by 'newAtoms'
double scaledChargeDifference(std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms)
{
  double charge = 0.0;
  for (const Atom& atom : newAtoms) charge += atom.scalingCoulomb * atom.charge;
  for (const Atom& atom : oldAtoms) charge -= atom.scalingCoulomb * atom.charge;
  return charge;
}

// per-group summed (unscaled) charge of the dU/dlambda group-tagged atoms of a set of spans;
// used as the 'external' charge derivative for the net-charge correction when one member of the
// group is updated (integer or switched-off atoms carry no tag and drop out automatically)
std::array<double, maximumNumberOfDUDlambdaGroups> groupChargeSumOfSpans(
    const std::vector<std::span<const Atom>>& spans)
{
  std::array<double, maximumNumberOfDUDlambdaGroups> charge{};
  for (std::span<const Atom> atoms : spans)
  {
    for (const Atom& atom : atoms)
    {
      if (atom.groupId != 0) charge[atom.groupId - 1] += atom.charge;
    }
  }
  return charge;
}

// external-field + framework + intermolecular energy difference for replacing 'oldAtoms' by 'newAtoms'
std::optional<RunningEnergy> computeNonEwaldEnergyDifference(System& system, std::span<const Atom> newAtoms,
                                                             std::span<const Atom> oldAtoms,
                                                             std::span<const Atom> moleculeAtomData)
{
  std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid, newAtoms,
      oldAtoms);
  if (!externalFieldDifference.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  if (!frameworkDifference.has_value()) return std::nullopt;

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, moleculeAtomData, newAtoms, oldAtoms);
  if (!moleculeDifference.has_value()) return std::nullopt;

  return externalFieldDifference.value() + frameworkDifference.value() + moleculeDifference.value();
}

// Brick-CFCMC-style aggregated tail-correction difference; the effective type counts are threaded
// across the sequential sub-steps of the move.
RunningEnergy computeTailEnergyDifference(System& system, std::vector<double>& tailEffectiveCounts,
                                          std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
                                          std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms)
{
  RunningEnergy result =
      Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
          system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms);
  return result;
}

std::pair<std::optional<RunningEnergy>, double3> groupSwapMoveCFCMCImplementation(RandomNumber& random, System& system,
                                                                                  std::size_t selectedComponent,
                                                                                  bool useCBMC)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  const Move::Types move = useCBMC ? Move::Types::GroupSwapCBCFCMC : Move::Types::GroupSwapCFCMC;
  Component& central = system.components[selectedComponent];

  if (central.groupComponentIds.empty())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  for (std::size_t satelliteComponentId : central.groupComponentIds)
  {
    if (satelliteComponentId >= system.components.size())
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
  }

  const std::vector<GroupMember> members = buildGroupMembers(system, selectedComponent, move);
  const std::size_t groupSize = members.size();

  // all group fractional molecules are coupled to the lambda histogram of the central component
  PropertyLambdaProbabilityHistogram& lambda = useCBMC ? central.lambdaGroupSwapCB : central.lambdaGroupSwap;
  const std::uint8_t dUdlambdaGroupId = lambda.dUdlambdaGroupId;
  const std::size_t oldBin = lambda.currentBin;
  const double deltaLambda = lambda.delta;
  const double maxChange = central.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  // current spans of all group fractional molecules
  std::vector<std::span<Atom>> fractionalMolecules(groupSize);
  for (std::size_t i = 0; i < groupSize; ++i)
  {
    fractionalMolecules[i] = system.spanOfMolecule(members[i].componentId, members[i].fractionalMoleculeIndex);
  }
  std::vector<std::vector<Atom>> oldFractionalMolecules(groupSize);
  for (std::size_t i = 0; i < groupSize; ++i)
  {
    oldFractionalMolecules[i].assign(fractionalMolecules[i].begin(), fractionalMolecules[i].end());
  }

  auto restoreFractionalGroup = [&]()
  {
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      std::copy(oldFractionalMolecules[i].begin(), oldFractionalMolecules[i].end(), fractionalMolecules[i].begin());
    }
  };

  // 'external' group-tagged charge of all group members except 'excluded', in their current state
  auto groupChargeSumOfOthers = [&](std::size_t excluded)
  {
    std::vector<std::span<const Atom>> spans;
    spans.reserve(groupSize);
    for (std::size_t j = 0; j < groupSize; ++j)
    {
      if (j != excluded) spans.push_back(fractionalMolecules[j]);
    }
    return groupChargeSumOfSpans(spans);
  };

  // Determine cutoff distances based on whether dual cutoff is used (only used for the CBMC grow/retrace).
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  // Snapshot the committed effective type counts; threaded across the sequential tail sub-steps of the move.
  std::vector<double> tailEffectiveCounts = system.effectiveNumberOfPseudoAtomsVDW;
  std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupCounts =
      system.fractionalPseudoAtomCountsPerGroup;

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints))  // Insertion move
  {
    // Steps for insertion Lambda_new = 1 + epsilon
    // ===================================================================
    // (1) Unbiased: the fractional group with lambda=lambda_old is made integer (lambda=1)
    // (2) A new fractional group is inserted with lambda_new = epsilon: at independent random
    //     positions (conventional) or grown with CBMC (CB variant)

    std::size_t newBin =
        static_cast<std::size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    central.mc_moves_statistics.addTrial(move, 0);

    RunningEnergy energyDifference{};
    double runningNetCharge = system.netCharge;
    bool firstEwaldCall = true;

    // (1) the fractional molecules of the group become integer, sequentially in place
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      for (Atom& atom : fractionalMolecules[i]) atom.setScalingToInteger();

      if (system.insideBlockedPockets(system.components[members[i].componentId], fractionalMolecules[i]))
      {
        restoreFractionalGroup();
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> nonEwaldDifference = computeNonEwaldEnergyDifference(
          system, fractionalMolecules[i], oldFractionalMolecules[i], system.spanOfMoleculeAtoms());
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
      if (!nonEwaldDifference.has_value())
      {
        restoreFractionalGroup();
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }
      energyDifference += nonEwaldDifference.value();

      time_begin = std::chrono::steady_clock::now();
      energyDifference += Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          firstEwaldCall ? system.storedEik : system.trialEik, system.trialEik, system.forceField,
          system.simulationBox, fractionalMolecules[i], oldFractionalMolecules[i], runningNetCharge,
          groupChargeSumOfOthers(i));
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
      firstEwaldCall = false;
      runningNetCharge += scaledChargeDifference(fractionalMolecules[i], oldFractionalMolecules[i]);

      time_begin = std::chrono::steady_clock::now();
      energyDifference += computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts,
                                                      fractionalMolecules[i], oldFractionalMolecules[i]);
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    }

    // (2) create the new fractional group with lambda_new
    std::vector<Molecule> trialMolecules(groupSize);
    std::vector<std::vector<Atom>> trialAtoms(groupSize);
    std::vector<RunningEnergy> trialGrowEnergies(groupSize);
    double rosenbluthRatio = 1.0;  // product of W/W_ideal (CB variant only; 1 for conventional)

    // accumulated background: existing molecules plus already-created trial molecules
    std::vector<Atom> accumulatedBackground(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());

    // trial energy contributions common to both variants (Ewald and tail; conventional adds its
    // non-Ewald part to 'energyDifference' directly)
    RunningEnergy trialExtraDifference{};

    for (std::size_t i = 0; i < groupSize; ++i)
    {
      const std::size_t componentId = members[i].componentId;
      Component& memberComponent = system.components[componentId];
      const std::size_t upcomingMoleculeId = system.numberOfMolecules() + i;

      if (useCBMC)
      {
        const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                            system.interpolationGrids, system.externalFieldInterpolationGrid,
                                            system.framework, system.spanOfFrameworkAtoms(), accumulatedBackground,
                                            system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

        time_begin = std::chrono::steady_clock::now();
        std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
            random, growContext, memberComponent, componentId, upcomingMoleculeId, newLambda, dUdlambdaGroupId, true);
        time_end = std::chrono::steady_clock::now();
        central.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
        system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);

        if (!growData || system.insideBlockedPockets(memberComponent, growData->atoms))
        {
          restoreFractionalGroup();
          return {std::nullopt, double3(0.0, 1.0, 0.0)};
        }

        if (system.forceField.useDualCutOff)
        {
          // correct the grown molecule from the inner cut-off to the full cut-offs, using the same
          // background as the growth
          std::optional<RunningEnergy> correction =
              CBMC::computeDualCutOffCorrection(growContext, memberComponent, growData->atoms);
          if (!correction.has_value())
          {
            restoreFractionalGroup();
            return {std::nullopt, double3(0.0, 1.0, 0.0)};
          }
          growData->energies += correction.value();
          growData->RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
        }

        rosenbluthRatio *= growData->RosenbluthWeight / memberComponent.idealGasRosenbluthWeight.value_or(1.0);
        trialMolecules[i] = growData->molecule;
        trialAtoms[i] = std::move(growData->atoms);
        trialGrowEnergies[i] = growData->energies;
      }
      else
      {
        std::pair<Molecule, std::vector<Atom>> trialMolecule =
            system.equilibratedIdealGasMoleculeRandomInBox(random, componentId);
        for (Atom& atom : trialMolecule.second)
        {
          atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeId);
          atom.componentId = static_cast<std::uint8_t>(componentId);
          atom.setScalingToFractional(newLambda, dUdlambdaGroupId);
        }

        if (system.insideBlockedPockets(memberComponent, trialMolecule.second))
        {
          restoreFractionalGroup();
          return {std::nullopt, double3(0.0, 1.0, 0.0)};
        }

        time_begin = std::chrono::steady_clock::now();
        std::optional<RunningEnergy> nonEwaldDifference =
            computeNonEwaldEnergyDifference(system, trialMolecule.second, {}, accumulatedBackground);
        time_end = std::chrono::steady_clock::now();
        central.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
        system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
        if (!nonEwaldDifference.has_value())
        {
          restoreFractionalGroup();
          return {std::nullopt, double3(0.0, 1.0, 0.0)};
        }
        energyDifference += nonEwaldDifference.value();

        trialMolecules[i] = trialMolecule.first;
        trialAtoms[i] = std::move(trialMolecule.second);
      }

      accumulatedBackground.insert(accumulatedBackground.end(), trialAtoms[i].begin(), trialAtoms[i].end());
    }

    // Ewald and tail-correction contributions of the new fractional group, threaded in growth order.
    // The 'external' group-tagged charge of a trial member consists of the trial members placed before it.
    time_begin = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      std::vector<std::span<const Atom>> earlierTrials;
      earlierTrials.reserve(i);
      for (std::size_t j = 0; j < i; ++j) earlierTrials.push_back(trialAtoms[j]);

      trialExtraDifference += Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.trialEik, system.forceField,
          system.simulationBox, trialAtoms[i], {}, runningNetCharge, groupChargeSumOfSpans(earlierTrials));
      runningNetCharge += scaledChargeDifference(trialAtoms[i], {});
    }
    time_end = std::chrono::steady_clock::now();
    central.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      trialExtraDifference +=
          computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, trialAtoms[i], {});
    }
    time_end = std::chrono::steady_clock::now();
    central.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    central.mc_moves_statistics.addConstructed(move, 0);

    std::vector<double3> electricFieldNeighborDelta;
    std::vector<std::vector<double3>> trialElectricFields(groupSize);
    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      // Step 1: making each fractional molecule integer rescales its own polarization coupling from
      // scalingCoulomb(lambda_old) to 1 (the field it feels is position-dependent only and unchanged).
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        std::span<double3> storedField =
            system.spanElectricFieldOld(members[i].componentId, members[i].fractionalMoleculeIndex);
        polarizationDifference += Interactions::computePolarizationEnergyDifference(
            system.forceField, storedField, storedField, fractionalMolecules[i], oldFractionalMolecules[i]);
      }

      // field felt by every new trial molecule (framework + Ewald; mutual trial-trial terms below)
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        trialElectricFields[i].assign(trialAtoms[i].size(), double3(0.0, 0.0, 0.0));
        Interactions::computeFrameworkMoleculeElectricFieldDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialElectricFields[i], {},
            trialAtoms[i], {});
        Interactions::computeEwaldFourierElectricFieldDifference(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.trialEik, system.forceField, system.simulationBox, trialElectricFields[i], {}, trialAtoms[i], {});
      }

      if (!system.forceField.omitInterPolarization)
      {
        electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));

        // step 1: making the fractional molecules integer changes the field they produce on the others
        for (std::size_t i = 0; i < groupSize; ++i)
        {
          accumulateInterMolecularFieldDelta(system, fractionalMolecules[i], oldFractionalMolecules[i],
                                             electricFieldNeighborDelta);
        }
        // step 2: field of the new trial molecules on the existing molecules
        for (std::size_t i = 0; i < groupSize; ++i)
        {
          accumulateInterMolecularFieldDelta(system, trialAtoms[i], {}, electricFieldNeighborDelta);
        }

        // mutual fields among the new trial molecules
        for (std::size_t i = 0; i < groupSize; ++i)
        {
          for (std::size_t j = i + 1; j < groupSize; ++j)
          {
            std::vector<double3> fieldOnIfromJ(trialAtoms[i].size());
            std::vector<double3> fieldOnJfromI(trialAtoms[j].size());
            [[maybe_unused]] std::optional<RunningEnergy> mutualEnergy =
                Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                    system.forceField, system.simulationBox, fieldOnIfromJ, fieldOnJfromI, std::span<double3>{},
                    trialAtoms[i], trialAtoms[j], {});
            for (std::size_t k = 0; k < trialElectricFields[i].size(); ++k)
              trialElectricFields[i][k] += fieldOnIfromJ[k];
            for (std::size_t k = 0; k < trialElectricFields[j].size(); ++k)
              trialElectricFields[j][k] += fieldOnJfromI[k];
          }
        }

        polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
            system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
            system.spanOfMoleculeAtoms());
      }

      for (std::size_t i = 0; i < groupSize; ++i)
      {
        polarizationDifference += Interactions::computePolarizationEnergyDifference(
            system.forceField, trialElectricFields[i], {}, trialAtoms[i], {});
      }
    }

    // Acceptance probability: one factor beta*f_c*V/(N_c + k) per group member, where k counts the
    // members of the same component placed so far (k = 1 for the first, 2 for the second, ...)
    double preFactor = 1.0;
    {
      std::vector<std::size_t> placedOfComponent(system.components.size(), 0);
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        const std::size_t componentId = members[i].componentId;
        const Component& memberComponent = system.components[componentId];
        const double fugacity =
            memberComponent.molFraction * memberComponent.fugacityCoefficient.value_or(1.0) * system.pressure;
        const std::size_t oldCount = system.numberOfIntegerMoleculesPerComponent[componentId];
        preFactor *= system.beta * fugacity * system.simulationBox.volume /
                     static_cast<double>(oldCount + 1 + placedOfComponent[componentId]);
        ++placedOfComponent[componentId];
      }
    }

    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    const double physicalPacc =
        preFactor * rosenbluthRatio *
        std::exp(-system.beta * (energyDifference.potentialEnergy() + trialExtraDifference.potentialEnergy() +
                                 polarizationDifference.potentialEnergy()));
    const double samplingPacc = physicalPacc * std::exp(biasTerm);

    if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN >= system.tmmc.maxMacrostate)
    {
      restoreFractionalGroup();
      return {std::nullopt, double3(0.0, 1.0 - physicalPacc, physicalPacc)};
    }

    const std::size_t newN = oldN == std::numeric_limits<std::size_t>::max() ? oldN : oldN + 1;
    const double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

    if (random.uniform() < biasTransitionMatrix * samplingPacc)
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

      lambda.setCurrentBin(newBin);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      // insert every trial molecule and swap it into its member's fractional slot (the former
      // fractional molecule, now integer, ends up at the end of the component's molecule list)
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        const std::size_t componentId = members[i].componentId;
        if (system.forceField.computePolarization)
        {
          system.insertMoleculePolarization(componentId, trialMolecules[i], trialAtoms[i], trialElectricFields[i]);
        }
        else
        {
          system.insertMolecule(componentId, trialMolecules[i], trialAtoms[i]);
        }
        const std::size_t lastMoleculeId = system.numberOfMoleculesPerComponent[componentId] - 1;
        std::span<Atom> lastMolecule = system.spanOfMolecule(componentId, lastMoleculeId);
        std::span<Atom> fractionalMolecule =
            system.spanOfMolecule(componentId, members[i].fractionalMoleculeIndex);
        std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
        if (system.forceField.computePolarization)
        {
          std::span<double3> fieldFractional =
              system.spanElectricFieldOld(componentId, members[i].fractionalMoleculeIndex);
          std::span<double3> fieldLast = system.spanElectricFieldOld(componentId, lastMoleculeId);
          std::swap_ranges(fieldFractional.begin(), fieldFractional.end(), fieldLast.begin());
        }
        std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentId, members[i].fractionalMoleculeIndex)],
                  system.moleculeData[system.moleculeIndexOfComponent(componentId, lastMoleculeId)]);
      }

      system.updateMoleculeAtomInformation();
      system.computeTailCorrectionCounts();

      central.mc_moves_statistics.addAccepted(move, 0);

      RunningEnergy totalDifference = energyDifference + trialExtraDifference + polarizationDifference;
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        totalDifference += trialGrowEnergies[i];
      }
      return {totalDifference, double3(0.0, 1.0 - physicalPacc, physicalPacc)};
    }

    restoreFractionalGroup();
    return {std::nullopt, double3(0.0, 1.0 - physicalPacc, physicalPacc)};
  }
  else if (selectedNewBin < 0)  // Deletion move
  {
    // Steps for deletion Lambda_new = -epsilon
    // ===================================================================
    // (1) The existing fractional group with lambda=lambda_old is removed (retraced with CBMC in
    //     the CB variant)
    // (2) Randomly selected integer molecules (one per member, drawn without replacement per
    //     component) become the new fractional group with lambda_new = 1 - epsilon

    central.mc_moves_statistics.addTrial(move, 1);

    // required integer molecules per component (a component occurring twice needs two)
    {
      std::vector<std::size_t> requiredOfComponent(system.components.size(), 0);
      for (const GroupMember& member : members) ++requiredOfComponent[member.componentId];
      for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
      {
        if (system.numberOfIntegerMoleculesPerComponent[componentId] < requiredOfComponent[componentId])
        {
          return {std::nullopt, double3(0.0, 1.0, 0.0)};
        }
      }
    }

    std::size_t newBin =
        static_cast<std::size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    // ordered selection without replacement: for every member draw a random integer molecule of its
    // component distinct from the ones already chosen for earlier members of the same component
    std::vector<std::size_t> selectedMolecules(groupSize);
    {
      std::vector<std::vector<std::size_t>> chosenPerComponent(system.components.size());
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        const std::size_t componentId = members[i].componentId;
        const std::size_t firstIntegerMolecule = system.numberOfFractionalMoleculesPerComponent[componentId];
        const std::size_t numberOfIntegers = system.numberOfIntegerMoleculesPerComponent[componentId];
        std::vector<std::size_t>& chosen = chosenPerComponent[componentId];

        std::size_t pick =
            static_cast<std::size_t>(random.uniform() * static_cast<double>(numberOfIntegers - chosen.size()));
        // map 'pick' to the actual molecule index, skipping already-chosen molecules
        std::size_t moleculeIndex = firstIntegerMolecule;
        for (std::size_t skipped = 0;; ++moleculeIndex)
        {
          if (std::find(chosen.begin(), chosen.end(), moleculeIndex) != chosen.end()) continue;
          if (skipped == pick) break;
          ++skipped;
        }
        selectedMolecules[i] = moleculeIndex;
        chosen.push_back(moleculeIndex);
      }
    }

    std::vector<std::span<Atom>> newFractionalMolecules(groupSize);
    std::vector<std::vector<Atom>> oldNewFractionalMolecules(groupSize);
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      newFractionalMolecules[i] = system.spanOfMolecule(members[i].componentId, selectedMolecules[i]);
      oldNewFractionalMolecules[i].assign(newFractionalMolecules[i].begin(), newFractionalMolecules[i].end());
    }

    auto restoreMolecules = [&]()
    {
      restoreFractionalGroup();
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        std::copy(oldNewFractionalMolecules[i].begin(), oldNewFractionalMolecules[i].end(),
                  newFractionalMolecules[i].begin());
      }
    };

    RunningEnergy energyDifference{};
    double runningNetCharge = system.netCharge;

    // (1a, CB variant) retrace the fractional group; the growth order of the reverse (insertion)
    // move is member 0 first, then 1, ..., each with the earlier members in the background, so
    // member i is retraced against the full background minus the members grown after it
    double rosenbluthRatio = 1.0;  // product of W_ideal/W_retrace (CB variant only; 1 for conventional)
    std::vector<RunningEnergy> retraceEnergies(groupSize);
    if (useCBMC)
    {
      time_begin = std::chrono::steady_clock::now();
      for (std::size_t i = groupSize; i-- > 0;)
      {
        std::vector<std::uint32_t> laterMoleculeIds;
        for (std::size_t j = i + 1; j < groupSize; ++j)
        {
          laterMoleculeIds.push_back(fractionalMolecules[j].front().moleculeId);
        }
        std::vector<Atom> background;
        background.reserve(system.spanOfMoleculeAtoms().size());
        for (const Atom& atom : system.spanOfMoleculeAtoms())
        {
          if (std::find(laterMoleculeIds.begin(), laterMoleculeIds.end(), atom.moleculeId) != laterMoleculeIds.end())
          {
            continue;
          }
          background.push_back(atom);
        }

        const CBMC::GrowContext retraceContext{system.hasExternalField, system.forceField, system.simulationBox,
                                               system.interpolationGrids, system.externalFieldInterpolationGrid,
                                               system.framework, system.spanOfFrameworkAtoms(), background,
                                               system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

        Component& memberComponent = system.components[members[i].componentId];
        ChainRetraceData retraceData =
            CBMC::retraceMoleculeSwapDeletion(random, retraceContext, memberComponent, fractionalMolecules[i]);

        if (system.forceField.useDualCutOff)
        {
          std::optional<RunningEnergy> correction =
              CBMC::computeDualCutOffCorrection(retraceContext, memberComponent, oldFractionalMolecules[i]);
          if (!correction.has_value())
          {
            restoreMolecules();
            return {std::nullopt, double3(0.0, 1.0, 0.0)};
          }
          retraceData.energies += correction.value();
          retraceData.RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
        }

        rosenbluthRatio *= memberComponent.idealGasRosenbluthWeight.value_or(1.0) / retraceData.RosenbluthWeight;
        retraceEnergies[i] = retraceData.energies;
      }
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
    }

    // (1b) Ewald and tail contributions of removing the fractional group, threaded in member order;
    // the 'external' group-tagged charge of member i consists of the not-yet-removed members j > i
    RunningEnergy removalExtraDifference{};
    time_begin = std::chrono::steady_clock::now();
    bool firstEwaldCall = true;
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      std::vector<std::span<const Atom>> laterMembers;
      for (std::size_t j = i + 1; j < groupSize; ++j) laterMembers.push_back(oldFractionalMolecules[j]);

      removalExtraDifference += Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          firstEwaldCall ? system.storedEik : system.trialEik, system.trialEik, system.forceField, system.simulationBox,
          {}, oldFractionalMolecules[i], runningNetCharge, groupChargeSumOfSpans(laterMembers));
      firstEwaldCall = false;
      runningNetCharge += scaledChargeDifference({}, oldFractionalMolecules[i]);
    }
    time_end = std::chrono::steady_clock::now();
    central.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      removalExtraDifference +=
          computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, {}, oldFractionalMolecules[i]);
    }
    time_end = std::chrono::steady_clock::now();
    central.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

    // (1c) switch the fractional group off; the conventional variant computes the real-space
    // (non-Ewald) removal energy sequentially in place (the CB variant has it in the retrace)
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      for (Atom& atom : fractionalMolecules[i]) atom.setScalingOff();

      if (!useCBMC)
      {
        time_begin = std::chrono::steady_clock::now();
        std::optional<RunningEnergy> nonEwaldDifference = computeNonEwaldEnergyDifference(
            system, fractionalMolecules[i], oldFractionalMolecules[i], system.spanOfMoleculeAtoms());
        time_end = std::chrono::steady_clock::now();
        central.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
        system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
        if (!nonEwaldDifference.has_value())
        {
          restoreMolecules();
          return {std::nullopt, double3(0.0, 1.0, 0.0)};
        }
        energyDifference += nonEwaldDifference.value();
      }
    }

    // (2) the selected integer molecules become fractional with lambda_new, sequentially in place;
    // the 'external' group-tagged charge of member i consists of the already-converted members j < i
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      for (Atom& atom : newFractionalMolecules[i])
      {
        atom.setScalingToFractional(newLambda, dUdlambdaGroupId);
      }

      if (system.insideBlockedPockets(system.components[members[i].componentId], newFractionalMolecules[i]))
      {
        restoreMolecules();
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> nonEwaldDifference = computeNonEwaldEnergyDifference(
          system, newFractionalMolecules[i], oldNewFractionalMolecules[i], system.spanOfMoleculeAtoms());
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
      if (!nonEwaldDifference.has_value())
      {
        restoreMolecules();
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }
      energyDifference += nonEwaldDifference.value();

      std::vector<std::span<const Atom>> earlierMembers;
      for (std::size_t j = 0; j < i; ++j) earlierMembers.push_back(newFractionalMolecules[j]);

      time_begin = std::chrono::steady_clock::now();
      energyDifference += Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.trialEik, system.forceField,
          system.simulationBox, newFractionalMolecules[i], oldNewFractionalMolecules[i], runningNetCharge,
          groupChargeSumOfSpans(earlierMembers));
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
      runningNetCharge += scaledChargeDifference(newFractionalMolecules[i], oldNewFractionalMolecules[i]);

      time_begin = std::chrono::steady_clock::now();
      energyDifference += computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts,
                                                      newFractionalMolecules[i], oldNewFractionalMolecules[i]);
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
    }

    central.mc_moves_statistics.addConstructed(move, 1);

    std::vector<double3> electricFieldNeighborDelta;
    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      // The removed group loses its (lambda-scaled) self polarization energy; the newly chosen
      // fractional group's polarization coupling is rescaled from 1 to scalingCoulomb(lambda_new)
      // (the fields they feel are unchanged).
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        std::span<double3> storedFieldFractional =
            system.spanElectricFieldOld(members[i].componentId, members[i].fractionalMoleculeIndex);
        polarizationDifference += Interactions::computePolarizationEnergyDifference(
            system.forceField, {}, storedFieldFractional, {}, oldFractionalMolecules[i]);

        std::span<double3> storedFieldNew =
            system.spanElectricFieldOld(members[i].componentId, selectedMolecules[i]);
        polarizationDifference += Interactions::computePolarizationEnergyDifference(
            system.forceField, storedFieldNew, storedFieldNew, newFractionalMolecules[i],
            oldNewFractionalMolecules[i]);
      }

      if (!system.forceField.omitInterPolarization)
      {
        electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
        for (std::size_t i = 0; i < groupSize; ++i)
        {
          accumulateInterMolecularFieldDelta(system, fractionalMolecules[i], oldFractionalMolecules[i],
                                             electricFieldNeighborDelta);
          accumulateInterMolecularFieldDelta(system, newFractionalMolecules[i], oldNewFractionalMolecules[i],
                                             electricFieldNeighborDelta);
        }

        // the removed molecules leave the system: their own field change must not contribute
        std::span<Atom> allAtoms = system.spanOfMoleculeAtoms();
        for (std::size_t i = 0; i < groupSize; ++i)
        {
          std::size_t offset = static_cast<std::size_t>(fractionalMolecules[i].data() - allAtoms.data());
          for (std::size_t k = 0; k < fractionalMolecules[i].size(); ++k)
          {
            electricFieldNeighborDelta[offset + k] = double3(0.0, 0.0, 0.0);
          }
        }

        polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
            system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
            system.spanOfMoleculeAtoms());
      }
    }

    // Acceptance probability: one factor (N_c - k)/(beta*f_c*V) per group member, where k counts
    // the members of the same component removed before it
    double preFactor = 1.0;
    {
      std::vector<std::size_t> removedOfComponent(system.components.size(), 0);
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        const std::size_t componentId = members[i].componentId;
        const Component& memberComponent = system.components[componentId];
        const double fugacity =
            memberComponent.molFraction * memberComponent.fugacityCoefficient.value_or(1.0) * system.pressure;
        const std::size_t oldCount = system.numberOfIntegerMoleculesPerComponent[componentId];
        preFactor *= static_cast<double>(oldCount - removedOfComponent[componentId]) /
                     (system.beta * fugacity * system.simulationBox.volume);
        ++removedOfComponent[componentId];
      }
    }

    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    const double physicalPacc =
        preFactor * rosenbluthRatio *
        std::exp(-system.beta * (energyDifference.potentialEnergy() + removalExtraDifference.potentialEnergy() +
                                 polarizationDifference.potentialEnergy()));
    const double samplingPacc = physicalPacc * std::exp(biasTerm);

    if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
    {
      restoreMolecules();
      return {std::nullopt, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
    }

    const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
    const double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

    if (random.uniform() < biasTransitionMatrix * samplingPacc)
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
      lambda.setCurrentBin(newBin);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      // per member: swap the newly chosen fractional molecule into the fractional slot and delete
      // the switched-off former fractional molecule. Deleting a molecule shifts the indices of the
      // later molecules of the same component, so the remaining selected indices are adjusted.
      std::vector<std::size_t> pendingSelected = selectedMolecules;
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        const std::size_t componentId = members[i].componentId;
        const std::size_t selectedMolecule = pendingSelected[i];

        std::span<Atom> fractionalMolecule =
            system.spanOfMolecule(componentId, members[i].fractionalMoleculeIndex);
        std::span<Atom> newFractionalMolecule = system.spanOfMolecule(componentId, selectedMolecule);
        std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
        if (system.forceField.computePolarization)
        {
          std::span<double3> storedFieldFractional =
              system.spanElectricFieldOld(componentId, members[i].fractionalMoleculeIndex);
          std::span<double3> storedFieldSelected = system.spanElectricFieldOld(componentId, selectedMolecule);
          std::swap_ranges(storedFieldFractional.begin(), storedFieldFractional.end(), storedFieldSelected.begin());
        }
        std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentId, selectedMolecule)],
                  system.moleculeData[system.moleculeIndexOfComponent(componentId, members[i].fractionalMoleculeIndex)]);
        system.deleteMolecule(componentId, selectedMolecule, newFractionalMolecule);

        for (std::size_t j = i + 1; j < groupSize; ++j)
        {
          if (members[j].componentId == componentId && pendingSelected[j] > selectedMolecule)
          {
            --pendingSelected[j];
          }
        }
      }
      system.computeTailCorrectionCounts();

      central.mc_moves_statistics.addAccepted(move, 1);

      RunningEnergy totalDifference = energyDifference + removalExtraDifference + polarizationDifference;
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        totalDifference -= retraceEnergies[i];
      }
      return {totalDifference, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
    }

    restoreMolecules();
    return {std::nullopt, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
  }
  else  // Lambda-change move
  {
    std::size_t newBin = static_cast<std::size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    central.mc_moves_statistics.addTrial(move, 2);

    RunningEnergy energyDifference{};
    double runningNetCharge = system.netCharge;
    bool firstEwaldCall = true;

    // rescale the fractional molecules of the group, sequentially in place
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      for (Atom& atom : fractionalMolecules[i]) atom.setScaling(newLambda);

      if (system.insideBlockedPockets(system.components[members[i].componentId], fractionalMolecules[i]))
      {
        restoreFractionalGroup();
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> nonEwaldDifference = computeNonEwaldEnergyDifference(
          system, fractionalMolecules[i], oldFractionalMolecules[i], system.spanOfMoleculeAtoms());
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::LambdaNonEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaNonEwald] += (time_end - time_begin);
      if (!nonEwaldDifference.has_value())
      {
        restoreFractionalGroup();
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }
      energyDifference += nonEwaldDifference.value();

      time_begin = std::chrono::steady_clock::now();
      energyDifference += Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          firstEwaldCall ? system.storedEik : system.trialEik, system.trialEik, system.forceField,
          system.simulationBox, fractionalMolecules[i], oldFractionalMolecules[i], runningNetCharge,
          groupChargeSumOfOthers(i));
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
      firstEwaldCall = false;
      runningNetCharge += scaledChargeDifference(fractionalMolecules[i], oldFractionalMolecules[i]);

      time_begin = std::chrono::steady_clock::now();
      energyDifference += computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts,
                                                      fractionalMolecules[i], oldFractionalMolecules[i]);
      time_end = std::chrono::steady_clock::now();
      central.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
    }

    central.mc_moves_statistics.addConstructed(move, 2);

    std::vector<double3> electricFieldNeighborDelta;
    if (system.forceField.computePolarization)
    {
      // Changing lambda rescales the polarization coupling of every fractional molecule of the
      // group (the fields they feel are unchanged because their positions do not move).
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        std::span<double3> storedField =
            system.spanElectricFieldOld(members[i].componentId, members[i].fractionalMoleculeIndex);
        energyDifference += Interactions::computePolarizationEnergyDifference(
            system.forceField, storedField, storedField, fractionalMolecules[i], oldFractionalMolecules[i]);
      }
    }
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        accumulateInterMolecularFieldDelta(system, fractionalMolecules[i], oldFractionalMolecules[i],
                                           electricFieldNeighborDelta);
      }
      energyDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }

    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
      central.mc_moves_statistics.addAccepted(move, 2);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      system.computeTailCorrectionCounts();

      lambda.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    }

    restoreFractionalGroup();
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}

}  // namespace

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::groupSwapMove_CFCMC(RandomNumber& random, System& system,
                                                                               std::size_t selectedComponent)
{
  return groupSwapMoveCFCMCImplementation(random, system, selectedComponent, false);
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::groupSwapMove_CFCMC_CBMC(RandomNumber& random,
                                                                                    System& system,
                                                                                    std::size_t selectedComponent)
{
  return groupSwapMoveCFCMCImplementation(random, system, selectedComponent, true);
}
