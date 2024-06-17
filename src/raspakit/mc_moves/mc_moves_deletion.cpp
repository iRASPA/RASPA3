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

module mc_moves_deletion;

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
import move_statistics;
import mc_moves_probabilities_particles;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMove(RandomNumber& random, System& system,
                                                                        size_t selectedComponent,
                                                                        size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.totalCounts += 1;

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.constructed += 1;
    system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.totalConstructed += 1;

    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // compute external field energy contribution
    std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, {}, molecule);
    if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // compute framework-molecule energy contribution
    std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {}, molecule);
    if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // compute molecule-molecule energy contribution
    std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), {}, molecule);
    if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, {}, molecule);
    time_end = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveEwald += (time_end - time_begin);
    system.mc_moves_cputime.swapDeletionMoveEwald += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveTail += (time_end - time_begin);
    system.mc_moves_cputime.swapDeletionMoveTail += (time_end - time_begin);

    // get the total difference in energy
    RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                     energyFourierDifference + tailEnergyDifference;

    double fugacity = system.components[selectedComponent].fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor =
        double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
        (system.beta * system.components[selectedComponent].molFraction * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
    size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    if (system.tmmc.doTMMC)
    {
      size_t newN = oldN - 1;
      if (newN < system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    // apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.accepted += 1;
      system.components[selectedComponent].mc_moves_statistics.swapDeletionMove.totalAccepted += 1;

      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {-energyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
