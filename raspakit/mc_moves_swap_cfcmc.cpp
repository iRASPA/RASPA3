module;

module mc_moves;

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


import component;
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
import move_statistics;
import mc_moves_probabilities_particles;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapMove_CFCMC(System& system, size_t selectedComponent, size_t selectedMolecule,
  bool insertionDisabled, bool deletionDisabled)
{
  PropertyLambdaProbabilityHistogram& lambda = system.components[selectedComponent].lambdaGC;
  size_t oldBin = lambda.currentBin;
  double deltaLambda = lambda.delta;
  double maxChange = system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.maxChange[2];
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(maxChange);
  size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  size_t indexFractionalMolecule = system.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfBins)) // Insertion move
  {
    if (insertionDisabled)
    {
      return {std::nullopt, double3(0.0, 0.0, 0.0)};
    }

    // Steps for insertion Lambda_new = 1 + epsilon
    // ===================================================================
    // (1) Unbiased: the fractional molecule with lambda=lambda_old is made integer (lambda=1) and deltaU is computed
    // (2) Unbiased: a new fractional molecule is inserted with lambda_new = epsilon

    size_t newBin = static_cast<size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfBins));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.counts[0] += 1;
    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalCounts[0] += 1;

    
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    // make copy of old fractional molecule for reference and restoring
    std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

    // fractional particle becomes integer (lambda=1.0)
    for (Atom& atom : fractionalMolecule)
    {
      atom.setScalingToInteger();
    }

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceStep1 = system.computeFrameworkMoleculeEnergyDifference(fractionalMolecule, oldFractionalMolecule);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (t2 - t1);
    system.mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (t2 - t1);

    if (!frameworkDifferenceStep1.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceStep1 = system.computeInterMolecularEnergyDifference(fractionalMolecule, oldFractionalMolecule);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (u2 - u1);
    system.mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (u2 - u1);

    if (!moleculeDifferenceStep1.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceStep1 = system.energyDifferenceEwaldFourier(system.storedEik, fractionalMolecule, oldFractionalMolecule);
    std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapLambdaInsertionMoveCFCMCEwald += (v2 - v1);
    system.mc_moves_cputime.swapLambdaInsertionMoveCFCMCEwald += (v2 - v1);

    RunningEnergy energyDifferenceStep1 = frameworkDifferenceStep1.value() + moleculeDifferenceStep1.value() + EwaldEnergyDifferenceStep1;

    // copy atoms from the old-fractional molecule, including the groupdIds
    std::vector<Atom> newatoms = system.components[selectedComponent].copyAtomsRandomlyRotatedAt(system.simulationBox.randomPosition(), oldFractionalMolecule, newLambda, system.numberOfMoleculesPerComponent[selectedComponent]);

    std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceStep2 = system.computeFrameworkMoleculeEnergyDifference(newatoms, {});
    std::chrono::system_clock::time_point t4 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (t4 - t3);
    system.mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (t4 - t3);

    if (!frameworkDifferenceStep2.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    std::chrono::system_clock::time_point u3 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceStep2 = system.computeInterMolecularEnergyDifference(newatoms, {});
    std::chrono::system_clock::time_point u4 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (u4 - u3);
    system.mc_moves_cputime.swapLambdaInsertionMoveCFCMCNonEwald += (u4 - u3);

    if (!moleculeDifferenceStep1.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    std::chrono::system_clock::time_point v3 = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceStep2 = system.energyDifferenceEwaldFourier(system.totalEik, newatoms, {});
    std::chrono::system_clock::time_point v4 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapLambdaInsertionMoveCFCMCEwald += (v4 - v3);
    system.mc_moves_cputime.swapLambdaInsertionMoveCFCMCEwald += (v4 - v3);
   
    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.constructed[0] += 1;
    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalConstructed[0] += 1;

    RunningEnergy energyDifferenceStep2= frameworkDifferenceStep2.value() + moleculeDifferenceStep2.value() + EwaldEnergyDifferenceStep2;

    //EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
    //                                    system.computeTailCorrectionVDWOldEnergy();

    double preFactor = system.beta * system.components[selectedComponent].molFraction * system.pressure * system.simulationBox.volume /
      static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    double Pacc = preFactor * exp(-system.beta * (energyDifferenceStep1.total() + energyDifferenceStep2.total()) + biasTerm);

    double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

    if(system.tmmc.doTMMC)
    {
      size_t newN = oldN + 1;
      if(newN > system.tmmc.maxMacrostate)
      {
        return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
      }
    }

    if (RandomNumber::Uniform() < biasTransitionMatrix  * Pacc)
    {
      system.acceptEwaldMove();

      system.components[selectedComponent].lambdaGC.setCurrentBin(newBin);

      // Note: inserting invalidates iterators and spans (the vector could reallocate memory)
      system.insertMolecule(selectedComponent, newatoms);

      // swap first and last molecule (selectedMolecule) so that molecule 'indexFractionalMolecule' is always the fractional molecule 
      size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
      std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
      fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());

      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.accepted[0] += 1;
      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalAccepted[0] += 1;

      return {energyDifferenceStep1 + energyDifferenceStep2, double3(0.0, 1.0 - Pacc, Pacc)};
    };

    // Restore old lamba
    std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  else if (selectedNewBin < 0) // Deletion move
  {
    if (deletionDisabled)
    {
      return {std::nullopt, double3(0.0, 0.0, 0.0)};
    }
    // Steps for deletion Lambda_new = -epsilon
    // ===================================================================
    // (1) Unbiased: the existing fractional molecule with lambda=lambda_o is removed.
    // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.

    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
    {
      selectedMolecule = system.randomIntegerMoleculeOfComponent(selectedComponent);

      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.counts[1] += 1;
      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalCounts[1] += 1;

      std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
      std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

      // make copy of old fractional molecule for reference and restoring
      std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());
      std::vector<Atom> oldNewFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      for (Atom& atom : fractionalMolecule)
      {
        atom.setScalingFullyOff();
        atom.groupId = std::byte{ 0 };
      }

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<RunningEnergy> frameworkDifferenceStep1 = system.computeFrameworkMoleculeEnergyDifference(fractionalMolecule, oldFractionalMolecule);
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      system.components[selectedComponent].mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (t2 - t1);
      system.mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (t2 - t1);

      if (!frameworkDifferenceStep1.has_value())
      {
        // reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
      std::optional<RunningEnergy> moleculeDifferenceStep1 = system.computeInterMolecularEnergyDifference(fractionalMolecule, oldFractionalMolecule);
      std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
      system.components[selectedComponent].mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (u2 - u1);
      system.mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (u2 - u1);

      if (!moleculeDifferenceStep1.has_value())
      {
        // reject, set fractional molecule back to old state
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
      RunningEnergy EwaldEnergyDifferenceStep1 = system.energyDifferenceEwaldFourier(system.storedEik, fractionalMolecule, oldFractionalMolecule);
      std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
      system.components[selectedComponent].mc_moves_cputime.swapLambdaDeletionMoveCFCMCEwald += (v2 - v1);

      RunningEnergy energyDifferenceStep1 = frameworkDifferenceStep1.value() + moleculeDifferenceStep1.value() + EwaldEnergyDifferenceStep1;

      std::vector<Atom> savedFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());

      // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.
      size_t newBin = static_cast<size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfBins));
      double newLambda = deltaLambda * static_cast<double>(newBin);

      // get the groupIds from the fractional molecule, set new Lambda
      std::transform(newFractionalMolecule.begin(), newFractionalMolecule.end(), oldFractionalMolecule.begin(), newFractionalMolecule.begin(),
        [newLambda](const Atom& a, const Atom& b) { return Atom(a.position, a.charge, newLambda, a.moleculeId, a.type, a.componentId, b.groupId); });

      std::chrono::system_clock::time_point v3 = std::chrono::system_clock::now();
      std::optional<RunningEnergy> frameworkDifferenceStep2 = system.computeFrameworkMoleculeEnergyDifference(newFractionalMolecule, savedFractionalMolecule);
      std::chrono::system_clock::time_point v4 = std::chrono::system_clock::now();
      system.components[selectedComponent].mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (v4 - v3);
      system.mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (v4 - v3);
      if (!frameworkDifferenceStep2.has_value())
      {
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      std::chrono::system_clock::time_point w3 = std::chrono::system_clock::now();
      std::optional<RunningEnergy> moleculeDifferenceStep2 = system.computeInterMolecularEnergyDifference(newFractionalMolecule, savedFractionalMolecule);
      std::chrono::system_clock::time_point w4 = std::chrono::system_clock::now();
      system.components[selectedComponent].mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (w4 - w3);
      system.mc_moves_cputime.swapLambdaDeletionMoveCFCMCNonEwald += (w4 - w3);
      if (!moleculeDifferenceStep2.has_value())
      {
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      }

      std::chrono::system_clock::time_point y3 = std::chrono::system_clock::now();
      RunningEnergy EwaldEnergyDifferenceStep2 = system.energyDifferenceEwaldFourier(system.totalEik, newFractionalMolecule, savedFractionalMolecule);
      std::chrono::system_clock::time_point y4 = std::chrono::system_clock::now();
      system.components[selectedComponent].mc_moves_cputime.swapLambdaDeletionMoveCFCMCEwald += (y4 - y3);
      system.mc_moves_cputime.swapLambdaDeletionMoveCFCMCEwald += (y4 - y3);

      //EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWRemoveEnergy(selectedComponent) - 
      //                                    system.computeTailCorrectionVDWOldEnergy();
      RunningEnergy tailEnergyDifferenceStep2;
      RunningEnergy energyDifferenceStep2 = frameworkDifferenceStep2.value() + moleculeDifferenceStep2.value() + EwaldEnergyDifferenceStep2 + tailEnergyDifferenceStep2;

      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.constructed[1] += 1;
      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalConstructed[1] += 1;

      double preFactor = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
        (system.beta * system.components[selectedComponent].molFraction * system.pressure * system.simulationBox.volume);
      double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
      double Pacc = preFactor *exp(-system.beta * (energyDifferenceStep1.total() + energyDifferenceStep2.total()) + biasTerm);

      double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

      if(system.tmmc.doTMMC)
      {
        size_t newN = oldN - 1;
        if(newN < system.tmmc.minMacrostate)
        {
          return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
        }
      }

      if (RandomNumber::Uniform() < biasTransitionMatrix * Pacc)
      {
        system.acceptEwaldMove();
        system.components[selectedComponent].lambdaGC.setCurrentBin(newBin);

        // Swap first and last molecule (selectedMolecule) so that molecule 'indexFractionalMolecule' is always the fractional molecule 
        std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());

        system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);

        system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.accepted[1] += 1;
        system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalAccepted[1] += 1;

        return {energyDifferenceStep1 + energyDifferenceStep2, double3(Pacc, 1.0 - Pacc, 0.0)};
      };

      // Restore the old- and the newly chosen fractional molecule
      std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
      std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());

      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
    }
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  else // Lambda-move
  {
    size_t newBin = static_cast<size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.counts[2] += 1;
    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalCounts[2] += 1;

    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);

    std::vector<Atom> trialPositions(molecule.begin(), molecule.end());
    std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
      [&](Atom a) { a.setScaling(newLambda); return a; });

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = system.computeFrameworkMoleculeEnergyDifference(trialPositions, molecule);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      system.components[selectedComponent].mc_moves_cputime.WidomMoveCFCMCNonEwald += (t2 - t1);
      system.mc_moves_cputime.WidomMoveCFCMCNonEwald += (t2 - t1);
    }
    else 
    {
      system.components[selectedComponent].mc_moves_cputime.swapLambdaChangeMoveCFCMCNonEwald += (t2 - t1);
      system.mc_moves_cputime.swapLambdaChangeMoveCFCMCNonEwald += (t2 - t1);
    }
   
    if (!frameworkEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = system.computeInterMolecularEnergyDifference(trialPositions, molecule);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      system.components[selectedComponent].mc_moves_cputime.WidomMoveCFCMCNonEwald += (u2 - u1);
      system.mc_moves_cputime.WidomMoveCFCMCNonEwald += (u2 - u1);
    }
    else 
    {
      system.components[selectedComponent].mc_moves_cputime.swapLambdaChangeMoveCFCMCNonEwald += (u2 - u1);
      system.mc_moves_cputime.swapLambdaChangeMoveCFCMCNonEwald += (u2 - u1);
    }
   

    if (!interEnergyDifference.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, trialPositions, molecule);
    std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
    if (insertionDisabled || deletionDisabled)
    {
      system.components[selectedComponent].mc_moves_cputime.WidomMoveCFCMCEwald += (v2 - v1);
      system.mc_moves_cputime.WidomMoveCFCMCEwald += (v2 - v1);
    }
    else 
    {
      system.components[selectedComponent].mc_moves_cputime.swapLambdaChangeMoveCFCMCEwald += (v2 - v1);
      system.mc_moves_cputime.swapLambdaChangeMoveCFCMCEwald += (v2 - v1);
    }
   

    RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() + EwaldFourierDifference;

    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.constructed[2] += 1;
    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalConstructed[2] += 1;

    double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    if (RandomNumber::Uniform() < std::exp(-system.beta * energyDifference.total() + biasTerm))
    {
      system.acceptEwaldMove();
      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.accepted[2] += 1;
      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC.totalAccepted[2] += 1;

      std::copy(trialPositions.begin(), trialPositions.end(), molecule.begin());

      system.components[selectedComponent].lambdaGC.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    };
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}

