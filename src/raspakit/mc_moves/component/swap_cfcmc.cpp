module;

module mc_moves_swap_cfcmc;

import std;

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapMove_CFCMC(RandomNumber& random, System& system,
                                                                          std::size_t selectedComponent,
                                                                          std::size_t selectedMolecule,
                                                                          bool insertionDisabled, bool deletionDisabled)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::SwapCFCMC;
  Component& component = system.components[selectedComponent];

  // Retrieve lambda parameters and select a new lambda bin for the move
  PropertyLambdaProbabilityHistogram& lambda = component.lambdaGC;
  std::size_t oldBin = lambda.currentBin;
  double deltaLambda = lambda.delta;
  double maxChange = component.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);
  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

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
    // (1) Unbiased: the fractional molecule with lambda=lambda_old is made integer (lambda=1) and deltaU is computed
    // (2) Unbiased: a new fractional molecule is inserted with lambda_new = epsilon

    std::size_t newBin =
        static_cast<std::size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    // Update counts for insertion move
    component.mc_moves_statistics.addTrial(move, 0);

    std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    // Make copy of old fractional molecule for reference and restoring
    std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

    // Fractional particle becomes integer (lambda=1.0)
    for (Atom& atom : fractionalMolecule)
    {
      atom.setScalingToInteger();
    }

    if ((system.insideBlockedPockets(component, fractionalMolecule)))
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> externalFieldDifferenceStep1 = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionExternalField] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionExternalField] += (time_end - time_begin);
    if (!externalFieldDifferenceStep1.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceStep1 = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionFramework] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionFramework] += (time_end - time_begin);
    if (!frameworkDifferenceStep1.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceStep1 = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), fractionalMolecule,
        oldFractionalMolecule);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionMolecule] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionMolecule] += (time_end - time_begin);
    if (!moleculeDifferenceStep1.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy contribution
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldEnergyDifferenceStep1 = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
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
    RunningEnergy tailEnergyDifference1 =
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

    RunningEnergy energyDifferenceStep1 = externalFieldDifferenceStep1.value() + frameworkDifferenceStep1.value() +
                                          moleculeDifferenceStep1.value() + EwaldEnergyDifferenceStep1 +
                                          tailEnergyDifference1;

    // Generate trial molecule for insertion
    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        component.equilibratedMoleculeRandomInBox(random, selectedComponent, system.simulationBox);

    if ((system.insideBlockedPockets(component, trialMolecule.second)))
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Copy atoms from the old fractional molecule, including the groupIds
    std::size_t upcomingMoleculeId = system.numberOfMolecules();

    std::uint8_t groupId = system.components[selectedComponent].lambdaGC.dUdlambdaGroupId;
    std::for_each(std::begin(trialMolecule.second), std::end(trialMolecule.second),
                  [selectedComponent, upcomingMoleculeId, groupId, newLambda](Atom& atom)
                  {
                    atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeId);
                    atom.componentId = static_cast<std::uint8_t>(selectedComponent);
                    atom.setScalingToFractional(newLambda, groupId);
                  });

    if ((system.insideBlockedPockets(component, trialMolecule.second)))
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> externalFieldDifferenceStep2 = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        trialMolecule.second, {});
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionExternalField] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionExternalField] += (time_end - time_begin);
    if (!externalFieldDifferenceStep2.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceStep2 = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, {});
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionFramework] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionFramework] += (time_end - time_begin);
    if (!frameworkDifferenceStep2.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceStep2 = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionMolecule] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionMolecule] += (time_end - time_begin);

    if (!moleculeDifferenceStep2.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy contribution
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldEnergyDifferenceStep2 = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, trialMolecule.second, {}, system.netCharge);
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);

    // Compute tail-correction energy contribution (threaded effective counts already include the step-1 change).
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifference2 =
        Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
            system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, trialMolecule.second, {}) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, {});
    time_end = std::chrono::steady_clock::now();
    component.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    component.mc_moves_statistics.addConstructed(move, 0);

    RunningEnergy energyDifferenceStep2 = externalFieldDifferenceStep2.value() + frameworkDifferenceStep2.value() +
                                          moleculeDifferenceStep2.value() + EwaldEnergyDifferenceStep2 +
                                          tailEnergyDifference2;

    // Polarization: (step 1) making the fractional molecule integer only changes the field it produces on the other
    // molecules (its own field is position-independent, so its self-polarization is unchanged); (step 2) inserting the
    // new fractional molecule adds its own polarization energy and again changes the field on the other molecules.
    std::vector<double3> electricFieldNeighborDelta;
    std::vector<double3> trialElectricField(trialMolecule.second.size());
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
                                                                    system.spanOfFrameworkAtoms(), trialElectricField,
                                                                    {}, trialMolecule.second, {});
      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.totalEik, system.forceField, system.simulationBox, trialElectricField, {}, trialMolecule.second, {});
      [[maybe_unused]] std::optional<RunningEnergy> e2 =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, trialElectricField,
              std::span<double3>{}, system.spanOfMoleculeAtoms(), trialMolecule.second, {});

      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, trialElectricField,
                                                                                 {}, trialMolecule.second, {}) +
                               Interactions::computePolarizationEnergyNeighborDifference(
                                   system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
                                   system.spanOfMoleculeAtoms());
    }
    else if (system.forceField.computePolarization)
    {
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), trialElectricField,
                                                                    {}, trialMolecule.second, {});
      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.totalEik, system.forceField, system.simulationBox, trialElectricField, {}, trialMolecule.second, {});
      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, trialElectricField,
                                                                                 {}, trialMolecule.second, {});
    }

    // Calculate acceptance probability
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor = system.beta * fugacity * system.simulationBox.volume / static_cast<double>(1 + oldN);
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    double physicalPacc = preFactor * std::exp(-system.beta * (energyDifferenceStep1.potentialEnergy() +
                                                               energyDifferenceStep2.potentialEnergy() +
                                                               polarizationDifference.potentialEnergy()));
    double samplingPacc = physicalPacc * std::exp(biasTerm);

    if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN >= system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - physicalPacc, physicalPacc)};
    }

    const std::size_t newN = oldN == std::numeric_limits<std::size_t>::max() ? oldN : oldN + 1;
    double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * samplingPacc)
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

      component.lambdaGC.setCurrentBin(newBin);

      // Commit the field changes on the surrounding molecules before the new fractional molecule (and its field) is
      // appended.
      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      // Note: inserting invalidates iterators and spans (the vector could reallocate memory)
      if (system.forceField.computePolarization)
      {
        system.insertMoleculePolarization(selectedComponent, trialMolecule.first, trialMolecule.second,
                                          trialElectricField);
      }
      else
      {
        system.insertMolecule(selectedComponent, trialMolecule.first, trialMolecule.second);
      }

      // Swap molecules to maintain fractional molecule index
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

      return {energyDifferenceStep1 + energyDifferenceStep2 + polarizationDifference,
              double3(0.0, 1.0 - physicalPacc, physicalPacc)};
    };

    // Restore old lambda
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
    // (1) Unbiased: the existing fractional molecule with lambda=lambda_o is removed.
    // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.

    component.mc_moves_statistics.addTrial(move, 1);

    if (oldN > 0)
    {
      selectedMolecule = system.randomIntegerMoleculeOfComponent(random, selectedComponent);

      std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

      // Make copies of molecules for restoring if needed
      std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());
      std::vector<Atom> oldNewFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // Set scaling to zero for deletion; the molecule remains the (parked) fractional slot
      for (Atom& atom : fractionalMolecule)
      {
        atom.setScalingToInactiveFractional();
      }

      if ((system.insideBlockedPockets(component, fractionalMolecule)))
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute external field energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> externalFieldDifferenceStep1 = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
          fractionalMolecule, oldFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionExternalField] += (time_end - time_begin);
      if (!externalFieldDifferenceStep1.has_value())
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> frameworkDifferenceStep1 = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), fractionalMolecule, oldFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionFramework] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionFramework] += (time_end - time_begin);
      if (!frameworkDifferenceStep1.has_value())
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute molecule-molecule energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> moleculeDifferenceStep1 = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), fractionalMolecule,
          oldFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionMolecule] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionMolecule] += (time_end - time_begin);
      if (!moleculeDifferenceStep1.has_value())
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute Ewald energy contribution
      time_begin = std::chrono::steady_clock::now();
      RunningEnergy EwaldEnergyDifferenceStep1 = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
          system.simulationBox, fractionalMolecule, oldFractionalMolecule, system.netCharge);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

      // Compute tail-correction energy contribution (Brick-CFCMC-style aggregated accounting).
      // Snapshot the committed effective type counts and thread them across the sequential sub-steps.
      time_begin = std::chrono::steady_clock::now();
      std::vector<double> tailEffectiveCounts = system.effectiveNumberOfPseudoAtomsVDW;
      std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupCounts =
          system.fractionalPseudoAtomCountsPerGroup;
      RunningEnergy tailEnergyDifference1 =
          Interactions::computeInterMolecularTailEnergyDifferenceAggregated(system.forceField, system.simulationBox,
                                                                            tailEffectiveCounts, tailGroupCounts,
                                                                            fractionalMolecule, oldFractionalMolecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                     system.spanOfFrameworkAtoms(), fractionalMolecule,
                                                                     oldFractionalMolecule);
      Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, fractionalMolecule,
                                              oldFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

      RunningEnergy energyDifferenceStep1 = externalFieldDifferenceStep1.value() + frameworkDifferenceStep1.value() +
                                            moleculeDifferenceStep1.value() + EwaldEnergyDifferenceStep1 +
                                            tailEnergyDifference1;

      std::vector<Atom> savedFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.
      std::size_t newBin =
          static_cast<std::size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
      double newLambda = deltaLambda * static_cast<double>(newBin);

      // Update new fractional molecule with new lambda
      std::uint8_t groupId = system.components[selectedComponent].lambdaGC.dUdlambdaGroupId;
      for (Atom& atom : newFractionalMolecule)
      {
        atom.setScalingToFractional(newLambda, groupId);
      }

      // Compute external field energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> externalFieldDifferenceStep2 = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
          newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionExternalField] += (time_end - time_begin);
      if (!externalFieldDifferenceStep2.has_value())
      {
        // Restore molecules and reject
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> frameworkDifferenceStep2 = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionFramework] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionFramework] += (time_end - time_begin);
      if (!frameworkDifferenceStep2.has_value())
      {
        // Restore molecules and reject
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute molecule-molecule energy contribution
      time_begin = std::chrono::steady_clock::now();
      std::optional<RunningEnergy> moleculeDifferenceStep2 = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newFractionalMolecule,
          savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionMolecule] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionMolecule] += (time_end - time_begin);
      if (!moleculeDifferenceStep2.has_value())
      {
        // Restore molecules and reject
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute Ewald energy contribution
      time_begin = std::chrono::steady_clock::now();
      RunningEnergy EwaldEnergyDifferenceStep2 = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
          system.simulationBox, newFractionalMolecule, savedFractionalMolecule, system.netCharge);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

      // Compute tail-correction energy contribution (threaded effective counts already include the step-1 change).
      time_begin = std::chrono::steady_clock::now();
      RunningEnergy tailEnergyDifferenceStep2 =
          Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
              system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, newFractionalMolecule,
              savedFractionalMolecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                     system.spanOfFrameworkAtoms(),
                                                                     newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::steady_clock::now();
      component.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

      RunningEnergy energyDifferenceStep2 = externalFieldDifferenceStep2.value() + frameworkDifferenceStep2.value() +
                                            moleculeDifferenceStep2.value() + EwaldEnergyDifferenceStep2 +
                                            tailEnergyDifferenceStep2;

      component.mc_moves_statistics.addConstructed(move, 1);

      // Polarization: the removed fractional molecule loses its own polarization energy (evaluated with its stored
      // field), while the field on all remaining molecules changes both because that molecule disappears (step 1) and
      // because the newly chosen fractional molecule is scaled down (step 2).
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

          // The molecule being removed must not appear in the neighbor sum; its own polarization change is accounted
          // for by 'removedSelf'.
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

      // Calculate acceptance probability
      double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
      double preFactor = double(oldN) / (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
      double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
      double physicalPacc = preFactor * std::exp(-system.beta * (energyDifferenceStep1.potentialEnergy() +
                                                                 energyDifferenceStep2.potentialEnergy() +
                                                                 polarizationDifference.potentialEnergy()));
      double samplingPacc = physicalPacc * std::exp(biasTerm);

      if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
      }

      const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
      double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

      // Apply acceptance/rejection rule
      if (random.uniform() < biasTransitionMatrix * samplingPacc)
      {
        Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
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

        // Swap molecules (and their fields) to maintain fractional molecule index
        std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
        if (system.forceField.computePolarization)
        {
          std::span<double3> fieldNewFractional = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
          std::span<double3> fieldFractional = system.spanElectricFieldOld(selectedComponent, indexFractionalMolecule);
          std::swap_ranges(fieldNewFractional.begin(), fieldNewFractional.end(), fieldFractional.begin());
        }
        std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)],
                  system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, 0)]);

        system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);
        system.computeTailCorrectionCounts();

        component.mc_moves_statistics.addAccepted(move, 1);

        return {energyDifferenceStep1 + energyDifferenceStep2 + polarizationDifference,
                double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
      };

      // Restore the old and the newly chosen fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

      return {std::nullopt, double3(physicalPacc, 1.0 - physicalPacc, 0.0)};
    }
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  else  // Lambda-move
  {
    std::size_t newBin = static_cast<std::size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    component.mc_moves_statistics.addTrial(move, 2);

    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    // Update molecule with new lambda
    std::vector<Atom> trialPositions(molecule.begin(), molecule.end());
    std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
                   [&](Atom a)
                   {
                     a.setScaling(newLambda);
                     return a;
                   });

    if (system.insideBlockedPockets(component, trialPositions))
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> externalFieldEnergyDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::ExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::ExternalField] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaExternalField] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaExternalField] += (time_end - time_begin);
    }
    if (!externalFieldEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Framework] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Framework] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaFramework] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaFramework] += (time_end - time_begin);
    }
    if (!frameworkEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Molecule] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Molecule] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaMolecule] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaMolecule] += (time_end - time_begin);
    }
    if (!interEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute Ewald energy contribution
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy EwaldFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, trialPositions, molecule, system.netCharge);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Ewald] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Ewald] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
    }

    // Compute tail-correction energy contribution (Brick-CFCMC-style aggregated accounting).
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
            system.forceField, system.simulationBox, system.effectiveNumberOfPseudoAtomsVDW,
            system.fractionalPseudoAtomCountsPerGroup, trialPositions, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::steady_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Tail] += (time_end - time_begin);
      system.mc_moves_cputime[Move::Types::WidomCFCMC][Move::Timing::Tail] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
      system.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
    }

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

    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    // Apply acceptance/rejection rule
    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      component.mc_moves_statistics.addAccepted(move, 2);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      std::copy(trialPositions.begin(), trialPositions.end(), molecule.begin());
      system.computeTailCorrectionCounts();

      component.lambdaGC.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    };
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}
