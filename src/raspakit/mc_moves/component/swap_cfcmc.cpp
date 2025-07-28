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

module mc_moves_swap_cfcmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

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
import scaling;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapMove_CFCMC(RandomNumber& random, System& system,
                                                                          std::size_t selectedComponent,
                                                                          std::size_t selectedMolecule,
                                                                          bool insertionDisabled, bool deletionDisabled)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapCFCMC;
  Component& component = system.components[selectedComponent];

  // Retrieve lambda parameters and select a new lambda bin for the move
  PropertyLambdaProbabilityHistogram& lambda = component.lambdaGC;
  std::size_t oldBin = lambda.currentBin;
  double deltaLambda = lambda.delta;
  double maxChange = component.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);
  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  std::size_t indexFractionalMolecule = system.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints))  // Insertion move
  {
    if (insertionDisabled)
    {
      return {std::nullopt, double3(0.0, 0.0, 0.0)};
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
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> externalFieldDifferenceStep1 = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-ExternalField"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-ExternalField"] += (time_end - time_begin);
    if (!externalFieldDifferenceStep1.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceStep1 = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Framework"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Framework"] += (time_end - time_begin);
    if (!frameworkDifferenceStep1.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceStep1 = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), fractionalMolecule,
        oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Molecule"] += (time_end - time_begin);
    if (!moleculeDifferenceStep1.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceStep1 = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);

    // Compute tail-correction energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference1 = Interactions::computeInterMolecularTailEnergyDifference(
                                              system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(),
                                              fractionalMolecule, oldFractionalMolecule) +
                                          Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                              system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
                                              fractionalMolecule, oldFractionalMolecule);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    RunningEnergy energyDifferenceStep1 = externalFieldDifferenceStep1.value() + frameworkDifferenceStep1.value() +
                                          moleculeDifferenceStep1.value() + EwaldEnergyDifferenceStep1 +
                                          tailEnergyDifference1;

    // Generate trial molecule for insertion
    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        component.equilibratedMoleculeRandomInBox(random, system.simulationBox);

    if ((system.insideBlockedPockets(component, trialMolecule.second)))
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Copy atoms from the old fractional molecule, including the groupIds
    std::size_t upcomingMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent];
    bool groupId = system.components[selectedComponent].lambdaGC.computeDUdlambda;
    std::for_each(std::begin(trialMolecule.second), std::end(trialMolecule.second),
                  [selectedComponent, upcomingMoleculeId, groupId, newLambda](Atom& atom)
                  {
                    atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeId);
                    atom.componentId = static_cast<std::uint8_t>(selectedComponent);
                    atom.groupId = groupId;
                    atom.isFractional = true;
                    atom.setScaling(newLambda);
                  });

    if ((system.insideBlockedPockets(component, trialMolecule.second)))
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute external field energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> externalFieldDifferenceStep2 = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-ExternalField"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-ExternalField"] += (time_end - time_begin);
    if (!externalFieldDifferenceStep2.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceStep2 = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Framework"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Framework"] += (time_end - time_begin);
    if (!frameworkDifferenceStep2.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceStep2 = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Molecule"] += (time_end - time_begin);

    if (!moleculeDifferenceStep2.has_value())
    {
      // Reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute Ewald energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceStep2 = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, trialMolecule.second, {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);

    // Compute tail-correction energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference2 =
        Interactions::computeInterMolecularTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {}) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, {});
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    component.mc_moves_statistics.addConstructed(move, 0);

    RunningEnergy energyDifferenceStep2 = externalFieldDifferenceStep2.value() + frameworkDifferenceStep2.value() +
                                          moleculeDifferenceStep2.value() + EwaldEnergyDifferenceStep2 +
                                          tailEnergyDifference2;

    // Calculate acceptance probability
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor =
        system.beta * component.molFraction * fugacity * system.simulationBox.volume / static_cast<double>(1 + oldN);
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    double Pacc =
        preFactor *
        std::exp(-system.beta * (energyDifferenceStep1.potentialEnergy() + energyDifferenceStep2.potentialEnergy()) +
                 biasTerm);

    double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

    if (system.tmmc.doTMMC)
    {
      std::size_t newN = oldN + 1;
      if (newN > system.tmmc.maxMacrostate)
      {
        return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
      }
    }

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

      component.lambdaGC.setCurrentBin(newBin);

      // Note: inserting invalidates iterators and spans (the vector could reallocate memory)
      system.insertMolecule(selectedComponent, trialMolecule.first, trialMolecule.second);

      // Swap molecules to maintain fractional molecule index
      std::size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
      std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
      fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
      std::swap(system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)],
                system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeId)]);

      component.mc_moves_statistics.addAccepted(move, 0);

      return {energyDifferenceStep1 + energyDifferenceStep2, double3(0.0, 1.0 - Pacc, Pacc)};
    };

    // Restore old lambda
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

      // Set scaling to zero for deletion
      for (Atom& atom : fractionalMolecule)
      {
        atom.setScalingFullyOff();
        atom.groupId = false;
        atom.isFractional = true;
      }

      if ((system.insideBlockedPockets(component, fractionalMolecule)))
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute external field energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> externalFieldDifferenceStep1 = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, fractionalMolecule, oldFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-ExternalField"] += (time_end - time_begin);
      if (!externalFieldDifferenceStep1.has_value())
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> frameworkDifferenceStep1 = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), fractionalMolecule, oldFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Framework"] += (time_end - time_begin);
      if (!frameworkDifferenceStep1.has_value())
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute molecule-molecule energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> moleculeDifferenceStep1 = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), fractionalMolecule,
          oldFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Molecule"] += (time_end - time_begin);
      if (!moleculeDifferenceStep1.has_value())
      {
        // Reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute Ewald energy contribution
      time_begin = std::chrono::system_clock::now();
      RunningEnergy EwaldEnergyDifferenceStep1 = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
          system.simulationBox, fractionalMolecule, oldFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);

      // Compute tail-correction energy contribution
      time_begin = std::chrono::system_clock::now();
      RunningEnergy tailEnergyDifference1 = Interactions::computeInterMolecularTailEnergyDifference(
                                                system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(),
                                                fractionalMolecule, oldFractionalMolecule) +
                                            Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                                system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
                                                fractionalMolecule, oldFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

      RunningEnergy energyDifferenceStep1 = externalFieldDifferenceStep1.value() + frameworkDifferenceStep1.value() +
                                            moleculeDifferenceStep1.value() + EwaldEnergyDifferenceStep1 +
                                            tailEnergyDifference1;

      std::vector<Atom> savedFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.
      std::size_t newBin =
          static_cast<std::size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
      double newLambda = deltaLambda * static_cast<double>(newBin);

      // Update new fractional molecule with new lambda
      bool groupId = system.components[selectedComponent].lambdaGC.computeDUdlambda;
      for (Atom& atom : newFractionalMolecule)
      {
        atom.scalingVDW = Scaling::scalingVDW(newLambda);
        atom.scalingCoulomb = Scaling::scalingCoulomb(newLambda);
        atom.groupId = groupId;
        atom.isFractional = true;
      }

      // Compute external field energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> externalFieldDifferenceStep2 = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, newFractionalMolecule,
          savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-ExternalField"] += (time_end - time_begin);
      if (!externalFieldDifferenceStep2.has_value())
      {
        // Restore molecules and reject
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> frameworkDifferenceStep2 = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Framework"] += (time_end - time_begin);
      if (!frameworkDifferenceStep2.has_value())
      {
        // Restore molecules and reject
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute molecule-molecule energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> moleculeDifferenceStep2 = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newFractionalMolecule,
          savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Molecule"] += (time_end - time_begin);
      if (!moleculeDifferenceStep2.has_value())
      {
        // Restore molecules and reject
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      // Compute Ewald energy contribution
      time_begin = std::chrono::system_clock::now();
      RunningEnergy EwaldEnergyDifferenceStep2 = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
          system.simulationBox, newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);

      // Compute tail-correction energy contribution
      time_begin = std::chrono::system_clock::now();
      RunningEnergy tailEnergyDifferenceStep2 =
          Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfMoleculeAtoms(), newFractionalMolecule,
                                                                  savedFractionalMolecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                     system.spanOfFrameworkAtoms(),
                                                                     newFractionalMolecule, savedFractionalMolecule);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

      RunningEnergy energyDifferenceStep2 = externalFieldDifferenceStep2.value() + frameworkDifferenceStep2.value() +
                                            moleculeDifferenceStep2.value() + EwaldEnergyDifferenceStep2 +
                                            tailEnergyDifferenceStep2;

      component.mc_moves_statistics.addConstructed(move, 1);

      // Calculate acceptance probability
      double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
      double preFactor = double(oldN) / (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
      double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
      double Pacc =
          preFactor *
          std::exp(-system.beta * (energyDifferenceStep1.potentialEnergy() + energyDifferenceStep2.potentialEnergy()) +
                   biasTerm);

      double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

      if (system.tmmc.doTMMC)
      {
        std::size_t newN = oldN - 1;
        if (newN < system.tmmc.minMacrostate)
        {
          return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
        }
      }

      // Apply acceptance/rejection rule
      if (random.uniform() < biasTransitionMatrix * Pacc)
      {
        Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
        component.lambdaGC.setCurrentBin(newBin);

        // Swap molecules to maintain fractional molecule index
        std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
        std::swap(system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)],
                  system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, 0)]);

        system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);

        component.mc_moves_statistics.addAccepted(move, 1);

        return {energyDifferenceStep1 + energyDifferenceStep2, double3(Pacc, 1.0 - Pacc, 0.0)};
      };

      // Restore the old and the newly chosen fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
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
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> externalFieldEnergyDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCFCMC]["ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCFCMC]["ExternalField"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-ExternalField"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-ExternalField"] += (time_end - time_begin);
    }
    if (!externalFieldEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCFCMC]["Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCFCMC]["Framework"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Framework"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Framework"] += (time_end - time_begin);
    }
    if (!frameworkEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCFCMC]["Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCFCMC]["Molecule"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Molecule"] += (time_end - time_begin);
    }
    if (!interEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute Ewald energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCFCMC]["Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCFCMC]["Ewald"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
    }

    // Compute tail-correction energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialPositions, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialPositions, molecule);
    time_end = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      component.mc_moves_cputime[MoveTypes::WidomCFCMC]["Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[MoveTypes::WidomCFCMC]["Tail"] += (time_end - time_begin);
    }
    else
    {
      component.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);
    }

    RunningEnergy energyDifference = externalFieldEnergyDifference.value() + frameworkEnergyDifference.value() +
                                     interEnergyDifference.value() + EwaldFourierDifference + tailEnergyDifference;

    component.mc_moves_statistics.addConstructed(move, 2);

    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    // Apply acceptance/rejection rule
    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      component.mc_moves_statistics.addAccepted(move, 2);

      std::copy(trialPositions.begin(), trialPositions.end(), molecule.begin());

      component.lambdaGC.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    };
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}
