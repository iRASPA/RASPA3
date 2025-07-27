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
#include <vector>
#endif

module mc_moves_deletion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import simd_quatd;
import component;
import atom;
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
import interactions_polarization;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMove(RandomNumber& random, System& system,
                                                                        std::size_t selectedComponent,
                                                                        std::size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::Swap;
  Component& component = system.components[selectedComponent];

  // Increment swap deletion move counts for the selected component
  component.mc_moves_statistics.addTrial(move, 1);

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // Copy the current electric field if polarization is computed
    std::vector<double3> electricFieldMoleculeOld(molecule.size());

    // Compute external field energy contribution
    std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, {}, molecule);
    if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy contribution
    std::optional<RunningEnergy> frameworkMolecule;
    if (system.forceField.computePolarization)
    {
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), {}, electricFieldMoleculeOld, {}, molecule);
    }
    else
    {
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), {}, molecule);
    }
    if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute molecule-molecule energy contribution
    std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), {}, molecule);
    if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute Ewald Fourier energy difference
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference;
    if (system.forceField.computePolarization)
    {
      energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.totalEik, system.forceField, system.simulationBox, {}, electricFieldMoleculeOld, {}, molecule);
    }
    else
    {
      energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
          system.simulationBox, {}, molecule);
    }
    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for Ewald calculations
    component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

    // Compute tail correction energy difference
    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for tail corrections
    component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      // Compute polarization energy difference
      polarizationDifference = Interactions::computePolarizationEnergyDifference(
          system.forceField, {}, electricFieldMoleculeOld, {}, molecule);
    }

    // Get the total difference in energy
    RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                     energyFourierDifference + tailEnergyDifference + polarizationDifference;

    // Increment constructed swap deletion move counts
    component.mc_moves_statistics.addConstructed(move, 1);

    // Calculate the acceptance probability
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
    std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    // Check if TMMC is enabled and if new state is below minimum macrostate
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
      // Move accepted; update acceptance statistics
      component.mc_moves_statistics.addAccepted(move, 1);

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
