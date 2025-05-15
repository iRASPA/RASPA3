module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

module mc_moves_swap_cfcmc_cbmc;

#ifndef USE_LEGACY_HEADERS
import <complex>;
import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
import <type_traits>;
#endif

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import energy_factor;
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
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapMove_CFCMC_CBMC(RandomNumber& random, System& system,
                                                                               size_t selectedComponent,
                                                                               size_t selectedMolecule,
                                                                               bool insertionDisabled,
                                                                               bool deletionDisabled)
{
  // Initialize time points for performance measurement
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapCBCFCMC;
  Component& component = system.components[selectedComponent];

  // Reference to the lambda histogram for the selected component
  PropertyLambdaProbabilityHistogram& lambda = component.lambdaGC;

  // Store old lambda values and bin index
  size_t oldBin = lambda.currentBin;
  double deltaLambda = lambda.delta;
  double oldLambda = component.lambdaGC.lambdaValue();

  // Get maximum allowed change in lambda for this move
  double maxChange = component.mc_moves_statistics.getMaxChange(move, 2);

  // Select a new bin based on the maximum change
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);

  // Store the current number of integer molecules
  size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  // Get index of the fractional molecule for the component
  size_t indexFractionalMolecule = system.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints))  // Insertion move
  {
    if (insertionDisabled)
    {
      return {std::nullopt, double3(0.0, 0.0, 0.0)};
    }

    // Steps for insertion Lambda_new = 1 + epsilon
    // ===================================================================
    // (1) Unbiased: the fractional molecule with lambda=lambda_old is made integer (lambda=1) and deltaU is
    // computed
    // (2) Biased: a new fractional molecule is grown with lambda_new = epsilon

    // Retrieve cutoff distances for interactions
    double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
    double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
    double cutOffCoulomb = system.forceField.cutOffCoulomb;
    // Get the grow type for the component
    Component::GrowType growType = component.growType;

    // Calculate the new bin index and lambda value
    size_t newBin = static_cast<size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    // Update move statistics for insertion move
    component.mc_moves_statistics.addTrial(move, 0);

    // Get the fractional molecule for the component
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    // Make a copy of the old fractional molecule for possible restoration
    std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

    // Fractional particle becomes integer (lambda=1.0)
    for (Atom& atom : fractionalMolecule)
    {
      atom.setScalingFullyOn();
    }

    // Check if the molecule is inside blocked pockets
    if ((system.insideBlockedPockets(component, fractionalMolecule)))
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-ExternalField"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-ExternalField"] += (time_end - time_begin);
    if (!externalFieldDifference.has_value())
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Framework"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Framework"] += (time_end - time_begin);
    if (!frameworkDifference.has_value())
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), fractionalMolecule,
        oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Molecule"] += (time_end - time_begin);
    if (!moleculeDifference.has_value())
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);

    // Compute tail-correction energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference = Interactions::computeInterMolecularTailEnergyDifference(
                                             system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(),
                                             fractionalMolecule, oldFractionalMolecule) +
                                         Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                             system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
                                             fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    // Sum up all energy differences
    RunningEnergy energyDifference = externalFieldDifference.value() + frameworkDifference.value() +
                                     moleculeDifference.value() + EwaldEnergyDifference + tailEnergyDifference;

    // Grow molecule with newLambda
    size_t newMolecule = system.numberOfMoleculesPerComponent[selectedComponent];

    time_begin = std::chrono::system_clock::now();
    std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
        random, component, system.hasExternalField, system.components, system.forceField, system.simulationBox,
        system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
        system.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, newMolecule,
        newLambda, static_cast<size_t>(oldFractionalMolecule.front().groupId), system.numberOfTrialDirections);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);

    if (!growData)
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Check if the new molecule is inside blocked pockets
    if ((system.insideBlockedPockets(component, growData->atom)))
    {
      // Reject move and restore the fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    component.mc_moves_statistics.addConstructed(move, 0);

    // Compute Ewald energy contribution for the new molecule
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, std::span(growData->atom.begin(), growData->atom.end()), {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);

    // Compute tail-correction energy contribution for the new molecule
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceGrow = Interactions::computeInterMolecularTailEnergyDifference(
                                                 system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(),
                                                 std::span(growData->atom.begin(), growData->atom.end()), {}) +
                                             Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                                 system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
                                                 std::span(growData->atom.begin(), growData->atom.end()), {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    // Calculate correction factor for Ewald summation
    double correctionFactorEwald = std::exp(
        -system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifferenceGrow.potentialEnergy()));

    // Compute acceptance probability
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * component.molFraction * fugacity *
                       system.simulationBox.volume /
                       static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    double Pacc = preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) *
                  exp(-system.beta * energyDifference.potentialEnergy() + biasTerm);

    // Retrieve bias from transition matrix
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

    if (system.tmmc.doTMMC)
    {
      size_t newN = oldN + 1;
      if (newN > system.tmmc.maxMacrostate)
      {
        return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
      }
    }

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      for (auto& atom : fractionalMolecule)
      {
        atom.groupId = uint8_t{0};
      }
      std::optional<RunningEnergy> energyDifferenceInterpolation =
          Interactions::computeFrameworkMoleculeEnergyDifferenceInterpolationExplicit(
              system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
              system.spanOfFrameworkAtoms(), fractionalMolecule);

      // Accept the move and update Ewald sums
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

      component.lambdaGC.setCurrentBin(newBin);

      // Insert the new molecule into the system
      system.insertMolecule(selectedComponent, growData->molecule, growData->atom);

      // Swap molecules to keep the fractional molecule at a fixed index
      size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
      std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
      fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
      std::swap(system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)],
                system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeId)]);

      component.mc_moves_statistics.addAccepted(move, 0);

      return {energyDifference + growData->energies + energyFourierDifference + tailEnergyDifferenceGrow -
                  energyDifferenceInterpolation.value(),
              double3(0.0, 1.0 - Pacc, Pacc)};
    };

    // Reject move and restore the fractional molecule
    std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  else if (selectedNewBin < 0)  // Deletion move
  {
    if (deletionDisabled)
    {
      return {std::nullopt, double3(0.0, 0.0, 0.0)};
    }
    // Steps for deletion Lambda_new = -epsilon
    // ===================================================================
    // (1) Biased: the existing fractional molecule is retraced using CBMC with lambda=lambda_o,
    //             fractional molecule is removed.
    // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.

    component.mc_moves_statistics.addTrial(move, 1);

    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
    {
      // Retrieve cutoff distances for interactions
      double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
      double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
      double cutOffCoulomb = system.forceField.cutOffCoulomb;

      // Select a random integer molecule to be fractional
      selectedMolecule = system.randomIntegerMoleculeOfComponent(random, selectedComponent);

      // Get spans for fractional and new fractional molecules
      std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

      // Make copies of the old molecules for restoration if needed
      std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());
      std::vector<Atom> oldNewFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // Retrace the existing fractional molecule
      time_begin = std::chrono::system_clock::now();
      ChainData retraceData = CBMC::retraceMoleculeSwapDeletion(
          random, component, system.hasExternalField, system.components, system.forceField, system.simulationBox,
          system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
          system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, indexFractionalMolecule,
          fractionalMolecule, oldLambda, system.numberOfTrialDirections);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);

      // Compute Ewald energy difference for the retraced molecule
      time_begin = std::chrono::system_clock::now();
      RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
          system.simulationBox, {}, fractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);

      // Compute tail-correction energy difference for the retraced molecule
      time_begin = std::chrono::system_clock::now();
      RunningEnergy tailEnergyDifferenceRetrace =
          Interactions::computeInterMolecularTailEnergyDifference(
              system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), {}, fractionalMolecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(
              system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {}, fractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

      // Calculate correction factor for Ewald summation
      double correctionFactorEwald = std::exp(
          -system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifferenceRetrace.potentialEnergy()));

      // Deactivate scaling for the fractional molecule
      for (Atom& atom : fractionalMolecule)
      {
        atom.setScalingFullyOff();
      }

      // Save the state of the new fractional molecule
      std::vector<Atom> savedFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // Calculate new bin and lambda value
      size_t newBin =
          static_cast<size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
      double newLambda = deltaLambda * static_cast<double>(newBin);

      // Update the new fractional molecule with the new lambda
      std::transform(newFractionalMolecule.begin(), newFractionalMolecule.end(), oldFractionalMolecule.begin(),
                     newFractionalMolecule.begin(), [newLambda](const Atom& a, const Atom& b)
                     { return Atom(a.position, a.charge, newLambda, a.moleculeId, a.type, a.componentId, b.groupId); });

      // Check if the new fractional molecule is inside blocked pockets
      if ((system.insideBlockedPockets(component, newFractionalMolecule)))
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute external field energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, newFractionalMolecule,
          savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-ExternalField"] += (time_end - time_begin);
      if (!externalFieldDifference.has_value())
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Framework"] += (time_end - time_begin);
      if (!frameworkDifference.has_value())
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute molecule-molecule energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newFractionalMolecule,
          savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Molecule"] += (time_end - time_begin);
      if (!moleculeDifference.has_value())
      {
        // Restore old molecules
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute Ewald energy contribution for the new fractional molecule
      time_begin = std::chrono::system_clock::now();
      RunningEnergy EwaldEnergyDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
          system.simulationBox, newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);

      // Compute tail-correction energy contribution for the new fractional molecule
      time_begin = std::chrono::system_clock::now();
      RunningEnergy tailEnergyDifference = Interactions::computeInterMolecularTailEnergyDifference(
                                               system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(),
                                               newFractionalMolecule, savedFractionalMolecule) +
                                           Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                               system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
                                               newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

      // Sum up all energy differences
      RunningEnergy energyDifference = externalFieldDifference.value() + frameworkDifference.value() +
                                       moleculeDifference.value() + EwaldEnergyDifference + tailEnergyDifference;

      component.mc_moves_statistics.addConstructed(move, 1);

      // Compute acceptance probability
      double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
      double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
      double preFactor = correctionFactorEwald *
                         double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                         (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
      double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
      double Pacc = preFactor * (idealGasRosenbluthWeight / retraceData.RosenbluthWeight) *
                    exp(-system.beta * energyDifference.potentialEnergy() + biasTerm);

      // Retrieve bias from transition matrix
      double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

      if (system.tmmc.doTMMC)
      {
        size_t newN = oldN - 1;
        if (newN < system.tmmc.minMacrostate)
        {
          return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
        }
      }

      // Apply acceptance/rejection rule
      if (random.uniform() < biasTransitionMatrix * Pacc)
      {
        for (auto& atom : fractionalMolecule)
        {
          atom.groupId = uint8_t{0};
        }
        for (auto& atom : savedFractionalMolecule)
        {
          atom.groupId = uint8_t{0};
        }

        std::optional<RunningEnergy> energyDifferenceInterpolation =
            Interactions::computeFrameworkMoleculeEnergyDifferenceInterpolationExplicit(
                system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                system.spanOfFrameworkAtoms(), savedFractionalMolecule);

        // Accept the move and update Ewald sums
        Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
        component.lambdaGC.setCurrentBin(newBin);

        // Swap molecules to keep the fractional molecule at a fixed index
        std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
        std::swap(
            system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)],
            system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)]);

        // Delete the selected molecule
        system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);

        component.mc_moves_statistics.addAccepted(move, 1);

        return {energyDifference + energyFourierDifference + tailEnergyDifferenceRetrace - retraceData.energies +
                    energyDifferenceInterpolation.value(),
                double3(Pacc, 1.0 - Pacc, 0.0)};
      };

      // Restore the old molecules
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
    }
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  else  // Lambda-move
  {
    // Calculate new bin and lambda value
    size_t newBin = static_cast<size_t>(selectedNewBin);
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
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> externalFieldEnergyDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["ExternalField"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-ExternalField"] += (time_end - time_begin);
    }
    if (!externalFieldEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy difference
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["Framework"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Framework"] += (time_end - time_begin);
    }
    if (!frameworkEnergyDifference.has_value())
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy difference
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["Molecule"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Molecule"] += (time_end - time_begin);
    }
    if (!interEnergyDifference.has_value())
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy difference
    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCBCFCMC]["Ewald"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
    }

    // Compute tail-correction energy difference
    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialPositions, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);

    // Sum up all energy differences
    RunningEnergy energyDifference = externalFieldEnergyDifference.value() + frameworkEnergyDifference.value() +
                                     interEnergyDifference.value() + EwaldFourierDifference + tailEnergyDifference;

    component.mc_moves_statistics.addConstructed(move, 2);

    // Calculate bias term for acceptance probability
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    // Apply acceptance/rejection rule
    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      // Accept the move and update Ewald sums
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      component.mc_moves_statistics.addAccepted(move, 2);

      // Update molecule positions with new scaling
      std::copy(trialPositions.begin(), trialPositions.end(), molecule.begin());

      component.lambdaGC.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    };

    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}
