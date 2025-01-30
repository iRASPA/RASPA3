module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_swap;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
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
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapMove(RandomNumber &random, System &selectedSystem, System &selectedSecondSystem,
                                 size_t selectedComponent, size_t &fractionalMoleculeSystem)
{
    if (random.uniform() < 0.5)
    {
        return MC_Moves::insertionMove(random, selectedSystem, selectedComponent);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
      return MC_Moves::deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);
    }
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMove(RandomNumber& random, System& system,
                                                                         size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove.totalCounts += 1;
  // Initialize selected molecule and update swap insertion move counts.

  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      system.components[selectedComponent].equilibratedMoleculeRandomInBox(random, system.simulationBox);
  // Generate a trial molecule with a random position inside the simulation box.

  if (system.insideBlockedPockets(system.components[selectedComponent], trialMolecule.second))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  // Check if the trial molecule is inside blocked pockets; reject if true.

  std::for_each(std::begin(trialMolecule.second), std::end(trialMolecule.second),
                [selectedComponent, selectedMolecule](Atom& atom)
                {
                  atom.moleculeId = static_cast<uint32_t>(selectedMolecule);
                  atom.componentId = static_cast<uint8_t>(selectedComponent);
                  atom.groupId = static_cast<uint8_t>(0);
                  atom.setScaling(1.0);
                });
  // Assign molecule ID, component ID, group ID, and set scaling factors for each atom.

  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove.constructed += 1;
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove.totalConstructed += 1;
  // Update constructed counts for swap insertion moves.

  // compute external field energy contribution
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, {});
  if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // compute framework-molecule energy contribution
  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, {});
  if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // compute molecule-molecule energy contribution
  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
  if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, trialMolecule.second, {});
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveEwald += (time_end - time_begin);
  system.mc_moves_cputime.swapInsertionMoveEwald += (time_end - time_begin);
  // Compute Ewald Fourier energy difference and update CPU time statistics.

  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), trialMolecule.second, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, {});
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveTail += (time_end - time_begin);
  system.mc_moves_cputime.swapInsertionMoveTail += (time_end - time_begin);
  // Compute tail energy difference and update CPU time statistics.

  // get the total difference in energy
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   energyFourierDifference + tailEnergyDifference;

  double fugacity = system.components[selectedComponent].fugacityCoefficient.value_or(1.0) * system.pressure;
  double preFactor = system.beta * system.components[selectedComponent].molFraction * fugacity *
                     system.simulationBox.volume /
                     double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
  double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
  size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);
  // Calculate acceptance probability and bias from the transition matrix.

  if (system.tmmc.doTMMC)
  {
    size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }
  // Check if the new macrostate exceeds the maximum allowed; reject if true.

  // apply acceptance/rejection rule
  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    system.components[selectedComponent].mc_moves_statistics.swapInsertionMove.accepted += 1;
    system.components[selectedComponent].mc_moves_statistics.swapInsertionMove.totalAccepted += 1;

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMolecule(selectedComponent, trialMolecule.first, trialMolecule.second);

    return {energyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
  };

  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}


std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMove(RandomNumber& random, System& system,
                                                                        size_t selectedComponent,
                                                                        size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  // Increment swap deletion move counts for the selected component
  system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.totalCounts += 1;

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    // Increment constructed swap deletion move counts
    system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.constructed += 1;
    system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.totalConstructed += 1;

    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // Compute external field energy contribution
    std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, {}, molecule);
    if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy contribution
    std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {}, molecule);
    if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute molecule-molecule energy contribution
    std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), {}, molecule);
    if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute Ewald Fourier energy difference
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, {}, molecule);
    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for Ewald calculations
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveEwald += (time_end - time_begin);
    system.mc_moves_cputime.swapDeletionMoveEwald += (time_end - time_begin);

    // Compute tail correction energy difference
    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for tail corrections
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveTail += (time_end - time_begin);
    system.mc_moves_cputime.swapDeletionMoveTail += (time_end - time_begin);

    // Get the total difference in energy
    RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                     energyFourierDifference + tailEnergyDifference;

    // Calculate the acceptance probability
    double fugacity = system.components[selectedComponent].fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor =
        double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
        (system.beta * system.components[selectedComponent].molFraction * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
    size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    // Check if TMMC is enabled and if new state is below minimum macrostate
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
      // Move accepted; update acceptance statistics
      system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.accepted += 1;
      system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.totalAccepted += 1;

      // Accept Ewald move and delete molecule from system
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {-energyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  // No molecules to delete; return default values
  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
