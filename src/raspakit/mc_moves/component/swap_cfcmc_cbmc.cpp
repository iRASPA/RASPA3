module;

module mc_moves_swap_cfcmc_cbmc;

import std;

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;
import scaling;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapMove_CFCMC_CBMC(RandomNumber& random, System& system,
                                                                               std::size_t selectedComponent,
                                                                               std::size_t selectedMolecule,
                                                                               bool insertionDisabled,
                                                                               bool deletionDisabled)
{
  // Initialize time points for performance measurement
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::SwapCBCFCMC;
  Component& component = system.components[selectedComponent];

  // Reference to the lambda histogram for the selected component
  PropertyLambdaProbabilityHistogram& lambda = component.lambdaGC;

  // Store old lambda values and bin index
  std::size_t oldBin = lambda.currentBin;
  double deltaLambda = lambda.delta;

  // Get maximum allowed change in lambda for this move
  double maxChange = component.mc_moves_statistics.getMaxChange(move, 2);

  // Select a new bin based on the maximum change
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);

  // Store the current number of integer molecules
  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  // Get index of the fractional molecule for the component
  std::size_t indexFractionalMolecule =
      system.indexOfFractionalMoleculeForMove(Move::Types::SwapCFCMC, selectedComponent);

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints))  // Insertion move
  {
    if (insertionDisabled)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Steps for insertion Lambda_new = 1 + epsilon
    // ===================================================================
    // (1) Unbiased: the fractional molecule with lambda=lambda_old is made integer (lambda=1) and deltaU is
    // computed
    // (2) Biased: a new fractional molecule is grown with lambda_new = epsilon

    // Determine cutoff distances based on whether dual cutoff is used.
    double cutOffFrameworkVDW =
        system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
    double cutOffMoleculeVDW =
        system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
    double cutOffCoulomb =
        system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
    // Calculate the new bin index and lambda value
    std::size_t newBin =
        static_cast<std::size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    // Update move statistics for insertion move
    component.mc_moves_statistics.addTrial(move, 0);

    // Get the fractional molecule for the component
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    // Make a copy of the old fractional molecule for possible restoration
    std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

    // Fractional particle becomes integer (lambda=1.0)
    // Before: lambda and fractional computed fully
    // After: 1.0 and integer computed using grids (or fully if without grids)
    for (Atom& atom : fractionalMolecule)
    {
      atom.setScalingToInteger();
    }

    // Check if the molecule is inside blocked pockets
    if ((system.insideBlockedPockets(component, fractionalMolecule)))
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionExternalField] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionExternalField] += (time_end - time_begin);
    if (!externalFieldDifference.has_value())
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionFramework] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionFramework] += (time_end - time_begin);
    if (!frameworkDifference.has_value())
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), fractionalMolecule,
        oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionMolecule] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionMolecule] += (time_end - time_begin);
    if (!moleculeDifference.has_value())
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy contribution
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldEnergyDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, fractionalMolecule, oldFractionalMolecule, system.netCharge);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);

    // Compute tail-correction energy contribution (Brick-CFCMC-style aggregated accounting).
    // Snapshot the committed effective type counts and thread them across the sequential sub-steps.
    time_begin = std::chrono::steady_clock::now();
    std::vector<double> tailEffectiveCounts = system.effectiveNumberOfPseudoAtomsVDW;
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupCounts =
        system.fractionalPseudoAtomCountsPerGroup;
    RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifferenceAggregated(system.forceField, system.simulationBox,
                                                                          tailEffectiveCounts, tailGroupCounts,
                                                                          fractionalMolecule, oldFractionalMolecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), fractionalMolecule,
                                                                   oldFractionalMolecule);
    Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, fractionalMolecule,
                                            oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    // Sum up all energy differences
    RunningEnergy energyDifference = externalFieldDifference.value() + frameworkDifference.value() +
                                     moleculeDifference.value() + EwaldEnergyDifference + tailEnergyDifference;

    // Grow molecule with newLambda (trial global molecule id, not a storage index)
    std::size_t newMolecule = system.numberOfMolecules();

    const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                        system.interpolationGrids, system.externalFieldInterpolationGrid,
                                        system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                        system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

    time_begin = std::chrono::steady_clock::now();
    std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
        random, growContext, component, selectedComponent, newMolecule, newLambda,
        system.components[selectedComponent].lambdaGC.dUdlambdaGroupId, true);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);

    if (!growData)
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    if (system.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
      // cut-offs, so that Rosenbluth weight and energies behave as if grown at the full cut-offs.
      std::optional<RunningEnergy> correctionNew =
          CBMC::computeDualCutOffCorrection(growContext, component, growData->atoms);
      if (!correctionNew.has_value())
      {
        // Reject move and restore the fractional molecule
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      growData->energies += correctionNew.value();
      growData->RosenbluthWeight *= std::exp(-system.beta * correctionNew->potentialEnergy());
    }

    // Check if the new molecule is inside blocked pockets
    if ((system.insideBlockedPockets(component, growData->atoms)))
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    component.mc_moves_statistics.addConstructed(move, 0);

    // Compute Ewald energy contribution for the new molecule
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.trialEik, system.forceField,
        system.simulationBox, std::span(growData->atoms.begin(), growData->atoms.end()), {}, system.netCharge);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);

    // Compute tail-correction energy contribution for the new molecule (threaded counts include the step-1 change).
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceGrow =
        Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
            system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts,
            std::span(growData->atoms.begin(), growData->atoms.end()), {}) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
            std::span(growData->atoms.begin(), growData->atoms.end()), {});
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    // Polarization: (step 1) making the fractional molecule integer only changes the field it produces on the other
    // molecules; (step 2) growing the new fractional molecule adds its own polarization energy and again changes the
    // field on the other molecules (its inter-molecular energy is already accounted for through CBMC).
    std::vector<double3> electricFieldNeighborDelta;
    std::vector<double3> grownElectricField(growData->atoms.size());
    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));

      std::vector<double3> fractionalFieldNew(fractionalMolecule.size());
      std::vector<double3> fractionalFieldOld(fractionalMolecule.size());
      [[maybe_unused]] std::optional<RunningEnergy> e1 =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, fractionalFieldNew,
              fractionalFieldOld, system.spanOfMoleculeAtoms(), fractionalMolecule, oldFractionalMolecule);

      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), grownElectricField,
                                                                    {}, growData->atoms, {});
      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.trialEik, system.forceField, system.simulationBox, grownElectricField, {}, growData->atoms, {});
      [[maybe_unused]] std::optional<RunningEnergy> e2 =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, grownElectricField,
              std::span<double3>{}, system.spanOfMoleculeAtoms(), growData->atoms, {});

      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, grownElectricField,
                                                                                 {}, growData->atoms, {}) +
                               Interactions::computePolarizationEnergyNeighborDifference(
                                   system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
                                   system.spanOfMoleculeAtoms());
    }
    else if (system.forceField.computePolarization)
    {
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), grownElectricField,
                                                                    {}, growData->atoms, {});
      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.trialEik, system.forceField, system.simulationBox, grownElectricField, {}, growData->atoms, {});
      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, grownElectricField,
                                                                                 {}, growData->atoms, {});
    }

    // Calculate correction factor for Ewald summation
    double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.potentialEnergy() +
                                                            tailEnergyDifferenceGrow.potentialEnergy() +
                                                            polarizationDifference.potentialEnergy()));

    // Compute acceptance probability
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * fugacity * system.simulationBox.volume /
                       static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    double physicalPacc = preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) *
                          std::exp(-system.beta * energyDifference.potentialEnergy());
    double samplingPacc = physicalPacc * std::exp(biasTerm);

    // Retrieve bias from transition matrix
    if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN >= system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - physicalPacc, physicalPacc)};
    }

    const std::size_t newN = oldN == std::numeric_limits<std::size_t>::max() ? oldN : oldN + 1;
    double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * samplingPacc)
    {
      // Accept the move and update Ewald sums
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

      component.lambdaGC.setCurrentBin(newBin);

      // Commit the field changes on the surrounding molecules before the new molecule (and its field) is appended.
      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      // Insert the new molecule into the system
      if (system.forceField.computePolarization)
      {
        system.insertMoleculePolarization(selectedComponent, growData->molecule, growData->atoms, grownElectricField);
      }
      else
      {
        system.insertMolecule(selectedComponent, growData->molecule, growData->atoms);
      }

      // Swap molecules to keep the fractional molecule at a fixed index
      std::size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
      std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
      fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
      if (system.forceField.computePolarization)
      {
        std::span<double3> fieldFractional = system.spanElectricFieldOld(selectedComponent, indexFractionalMolecule);
        std::span<double3> fieldLast = system.spanElectricFieldOld(selectedComponent, lastMoleculeId);
        std::swap_ranges(fieldFractional.begin(), fieldFractional.end(), fieldLast.begin());
      }
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)],
                system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeId)]);

      system.updateMoleculeAtomInformation();
      system.computeTailCorrectionCounts();

      component.mc_moves_statistics.addAccepted(move, 0);

      return {energyDifference + growData->energies + energyFourierDifference + tailEnergyDifferenceGrow +
                  polarizationDifference,
              double3(0.0, 1.0 - physicalPacc, physicalPacc)};
    };

    // Reject move and restore the fractional molecule
    std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

    return {std::nullopt, double3(0.0, 1.0 - physicalPacc, physicalPacc)};
  }
  else if (selectedNewBin < 0)  // Deletion move
  {
    if (deletionDisabled)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    // Steps for deletion Lambda_new = -epsilon
    // ===================================================================
    // (1) Biased: the existing fractional molecule is retraced using CBMC with lambda=lambda_o,
    //             fractional molecule is removed.
    // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.

    component.mc_moves_statistics.addTrial(move, 1);

    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
    {
      // Determine cutoff distances based on whether dual cutoff is used.
      double cutOffFrameworkVDW =
          system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
      double cutOffMoleculeVDW =
          system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
      double cutOffCoulomb =
          system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
      // Select a random integer molecule to be fractional
      selectedMolecule = system.randomIntegerMoleculeOfComponent(random, selectedComponent);

      // Get spans for fractional and new fractional molecules
      std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

      // Make copies of the old molecules for restoration if needed
      std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());
      std::vector<Atom> oldNewFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      const CBMC::GrowContext retraceContext{system.hasExternalField, system.forceField, system.simulationBox,
                                             system.interpolationGrids, system.externalFieldInterpolationGrid,
                                             system.framework, system.spanOfFrameworkAtoms(),
                                             system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW,
                                             cutOffMoleculeVDW, cutOffCoulomb};

      // Retrace the existing fractional molecule
      time_begin = std::chrono::steady_clock::now();
      ChainRetraceData retraceData =
          CBMC::retraceMoleculeSwapDeletion(random, retraceContext, component, fractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);

      if (system.forceField.useDualCutOff)
      {
        // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
        // cut-offs, so that Rosenbluth weight and energies behave as if retraced at the full cut-offs.
        std::optional<RunningEnergy> correctionOld =
            CBMC::computeDualCutOffCorrection(retraceContext, component, oldFractionalMolecule);
        if (!correctionOld.has_value())
        {
          return {std::nullopt, double3(0.0, 1.0, 0.0)};
        }

        retraceData.energies += correctionOld.value();
        retraceData.RosenbluthWeight *= std::exp(-system.beta * correctionOld->potentialEnergy());
      }

      // Compute Ewald energy difference for the retraced molecule
      time_begin = std::chrono::steady_clock::now();
      RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
          system.simulationBox, {}, fractionalMolecule, system.netCharge);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

      // Compute tail-correction energy difference for the retraced molecule (Brick-CFCMC-style aggregated accounting).
      // Snapshot the committed effective type counts and thread them across the sequential sub-steps.
      time_begin = std::chrono::steady_clock::now();
      std::vector<double> tailEffectiveCounts = system.effectiveNumberOfPseudoAtomsVDW;
      std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupCounts =
          system.fractionalPseudoAtomCountsPerGroup;
      RunningEnergy tailEnergyDifferenceRetrace =
          Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
              system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, {}, fractionalMolecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(
              system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {}, fractionalMolecule);
      Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, {}, fractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

      // Calculate correction factor for Ewald summation
      double correctionFactorEwald = std::exp(
          -system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifferenceRetrace.potentialEnergy()));

      // Deactivate scaling for the fractional molecule; it remains the (parked) fractional slot
      // Before: lambda and fractional computed fully
      // After: 0.0 and integer computed using grids (or fully if without grids)
      for (Atom& atom : fractionalMolecule)
      {
        atom.setScalingToInactiveFractional();
      }

      // Save the state of the new fractional molecule
      std::vector<Atom> savedFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // Calculate new bin and lambda value
      std::size_t newBin =
          static_cast<std::size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
      double newLambda = deltaLambda * static_cast<double>(newBin);

      // Update the new fractional molecule with the new lambda
      std::uint8_t groupId = system.components[selectedComponent].lambdaGC.dUdlambdaGroupId;
      for (Atom& atom : newFractionalMolecule)
      {
        atom.setScalingToFractional(newLambda, groupId);
      }

      // Check if the new fractional molecule is inside blocked pockets
      if ((system.insideBlockedPockets(component, newFractionalMolecule)))
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute external field energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
          newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionExternalField] += (time_end - time_begin);
      if (!externalFieldDifference.has_value())
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionFramework] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionFramework] += (time_end - time_begin);
      if (!frameworkDifference.has_value())
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute molecule-molecule energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newFractionalMolecule,
          savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionMolecule] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionMolecule] += (time_end - time_begin);
      if (!moleculeDifference.has_value())
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute Ewald energy contribution for the new fractional molecule
      time_begin = std::chrono::steady_clock::now();
      RunningEnergy EwaldEnergyDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.trialEik, system.forceField,
          system.simulationBox, newFractionalMolecule, savedFractionalMolecule, system.netCharge);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

      // Compute tail-correction energy contribution for the new fractional molecule (threaded counts include step 1).
      time_begin = std::chrono::steady_clock::now();
      RunningEnergy tailEnergyDifference =
          Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
              system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, newFractionalMolecule,
              savedFractionalMolecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                     system.spanOfFrameworkAtoms(), newFractionalMolecule,
                                                                     savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

      // Sum up all energy differences
      RunningEnergy energyDifference = externalFieldDifference.value() + frameworkDifference.value() +
                                       moleculeDifference.value() + EwaldEnergyDifference + tailEnergyDifference;

      component.mc_moves_statistics.addConstructed(move, 1);

      // Polarization: the removed fractional molecule loses its own polarization energy (evaluated with its stored
      // field), while the field on all remaining molecules changes because that molecule disappears and because the
      // newly chosen fractional molecule is scaled down.
      std::vector<double3> electricFieldNeighborDelta;
      RunningEnergy polarizationDifference;
      if (system.forceField.computePolarization)
      {
        std::span<double3> fractionalStoredField =
            system.spanElectricFieldOld(selectedComponent, indexFractionalMolecule);
        RunningEnergy removedSelf = Interactions::computePolarizationEnergyDifference(
            system.forceField, {}, fractionalStoredField, {}, fractionalMolecule);

        if (!system.forceField.omitInterPolarization)
        {
          electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));

          std::vector<double3> fractionalFieldNew(fractionalMolecule.size());
          std::vector<double3> fractionalFieldOld(fractionalMolecule.size());
          [[maybe_unused]] std::optional<RunningEnergy> e1 =
              Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                  system.forceField, system.simulationBox, electricFieldNeighborDelta, fractionalFieldNew,
                  fractionalFieldOld, system.spanOfMoleculeAtoms(), fractionalMolecule, oldFractionalMolecule);

          std::vector<double3> newFractionalFieldNew(newFractionalMolecule.size());
          std::vector<double3> newFractionalFieldOld(newFractionalMolecule.size());
          [[maybe_unused]] std::optional<RunningEnergy> e2 =
              Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                  system.forceField, system.simulationBox, electricFieldNeighborDelta, newFractionalFieldNew,
                  newFractionalFieldOld, system.spanOfMoleculeAtoms(), newFractionalMolecule, savedFractionalMolecule);

          std::span<Atom> allAtoms = system.spanOfMoleculeAtoms();
          std::size_t fractionalOffset = static_cast<std::size_t>(fractionalMolecule.data() - allAtoms.data());
          for (std::size_t k = 0; k < fractionalMolecule.size(); ++k)
          {
            electricFieldNeighborDelta[fractionalOffset + k] = double3(0.0, 0.0, 0.0);
          }

          polarizationDifference = removedSelf + Interactions::computePolarizationEnergyNeighborDifference(
                                                     system.forceField, system.spanOfMoleculeElectricField(),
                                                     electricFieldNeighborDelta, system.spanOfMoleculeAtoms());
        }
        else
        {
          polarizationDifference = removedSelf;
        }
      }

      // Compute acceptance probability
      double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
      double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
      double preFactor = correctionFactorEwald *
                         double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                         (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
      double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
      double physicalPacc =
          preFactor * (idealGasRosenbluthWeight / retraceData.RosenbluthWeight) *
          std::exp(-system.beta * (energyDifference.potentialEnergy() + polarizationDifference.potentialEnergy()));
      double samplingPacc = physicalPacc * std::exp(biasTerm);

      // Retrieve bias from transition matrix
      if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
      }

      const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
      double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

      // Apply acceptance/rejection rule
      if (random.uniform() < biasTransitionMatrix * samplingPacc)
      {
        // Accept the move and update Ewald sums
        Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
        component.lambdaGC.setCurrentBin(newBin);

        // Commit the field changes on the remaining molecules (the removed molecule's own entries were zeroed).
        if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
        {
          std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
          for (std::size_t i = 0; i < storedElectricField.size(); ++i)
          {
            storedElectricField[i] += electricFieldNeighborDelta[i];
          }
        }

        // Swap molecules (and their fields) to keep the fractional molecule at a fixed index
        std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
        if (system.forceField.computePolarization)
        {
          std::span<double3> fieldNewFractional = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
          std::span<double3> fieldFractional = system.spanElectricFieldOld(selectedComponent, indexFractionalMolecule);
          std::swap_ranges(fieldNewFractional.begin(), fieldNewFractional.end(), fieldFractional.begin());
        }
        std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)],
                  system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)]);

        // Delete the selected molecule
        system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);
        system.computeTailCorrectionCounts();

        component.mc_moves_statistics.addAccepted(move, 1);

        return {energyDifference + energyFourierDifference + tailEnergyDifferenceRetrace - retraceData.energies +
                    polarizationDifference,
                double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
      };

      // Restore the old molecules
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

      return {std::nullopt, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
    }
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  else  // Lambda-move
  {
    // Calculate new bin and lambda value
    std::size_t newBin = static_cast<std::size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    component.mc_moves_statistics.addTrial(move, 2);

    // Get the fractional molecule
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    // Create trial positions with new lambda scaling
    std::vector<Atom> trialPositions(molecule.begin(), molecule.end());
    std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
                   [&](Atom a)
                   {
                     a.setScaling(newLambda);
                     return a;
                   });

    // Check if the trial positions are inside blocked pockets
    if (system.insideBlockedPockets(component, trialPositions))
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy difference
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> externalFieldEnergyDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::ExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::ExternalField] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaExternalField] += (time_end - time_begin);
    }
    if (!externalFieldEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy difference
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::Framework] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::Framework] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaFramework] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaFramework] += (time_end - time_begin);
    }
    if (!frameworkEnergyDifference.has_value())
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy difference
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::Molecule] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::Molecule] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaMolecule] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaMolecule] += (time_end - time_begin);
    }
    if (!interEnergyDifference.has_value())
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy difference
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, trialPositions, molecule, system.netCharge);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::Ewald] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCBCFCMC][Move::Timing::Ewald] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
    }

    // Compute tail-correction energy difference (Brick-CFCMC-style aggregated accounting)
    time_begin = std::chrono::steady_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
            system.forceField, system.simulationBox, system.effectiveNumberOfPseudoAtomsVDW,
            system.fractionalPseudoAtomCountsPerGroup, trialPositions, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);

    // Sum up all energy differences
    RunningEnergy energyDifference = externalFieldEnergyDifference.value() + frameworkEnergyDifference.value() +
                                     interEnergyDifference.value() + EwaldFourierDifference + tailEnergyDifference;

    // Polarization: changing lambda only rescales the field the fractional molecule produces on the other molecules
    // (its own field, and therefore its own polarization energy, is unchanged because its position does not move).
    std::vector<double3> electricFieldNeighborDelta;
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));

      std::vector<double3> fractionalFieldNew(molecule.size());
      std::vector<double3> fractionalFieldOld(molecule.size());
      [[maybe_unused]] std::optional<RunningEnergy> e =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, fractionalFieldNew,
              fractionalFieldOld, system.spanOfMoleculeAtoms(), trialPositions, molecule);

      energyDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }

    component.mc_moves_statistics.addConstructed(move, 2);

    // Calculate bias term for acceptance probability
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    // Apply acceptance/rejection rule
    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      // Accept the move and update Ewald sums
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
      component.mc_moves_statistics.addAccepted(move, 2);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      // Update molecule positions with new scaling
      std::copy(trialPositions.begin(), trialPositions.end(), molecule.begin());
      system.computeTailCorrectionCounts();

      component.lambdaGC.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    };

    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}
