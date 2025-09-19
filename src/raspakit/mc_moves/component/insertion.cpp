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

module mc_moves_insertion;

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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMove(RandomNumber& random, System& system,
                                                                         std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::Swap;
  Component& component = system.components[selectedComponent];

  // Initialize selected molecule and update swap insertion move counts.
  std::size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  component.mc_moves_statistics.addTrial(move, 0);

  // Generate a trial molecule with a random position inside the simulation box.
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      component.equilibratedMoleculeRandomInBox(random, system.simulationBox);

  std::vector<double3> electricFieldMoleculeNew(trialMolecule.second.size());

  // Check if the trial molecule is inside blocked pockets; reject if true.
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Assign molecule ID, component ID, group ID, and set scaling factors for each atom.
  std::for_each(std::begin(trialMolecule.second), std::end(trialMolecule.second),
                [selectedComponent, selectedMolecule](Atom& atom)
                {
                  atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
                  atom.componentId = static_cast<std::uint8_t>(selectedComponent);
                  atom.groupId = static_cast<std::uint8_t>(0);
                  atom.setScaling(1.0);
                });

  // Update constructed counts for swap insertion moves.
  component.mc_moves_statistics.addConstructed(move, 0);

  // compute external field energy contribution
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, {});
  if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // compute framework-molecule energy contribution
  std::optional<RunningEnergy> frameworkMolecule;
  if (system.forceField.computePolarization)
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, {}, trialMolecule.second, {});
  }
  else
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, {});
  }
  if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // compute molecule-molecule energy contribution
  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
  if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // Compute Ewald Fourier energy difference and update CPU time statistics.
  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference;
  if (system.forceField.computePolarization)
  {
    energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, {}, trialMolecule.second,
        {});
  }
  else
  {
    energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, trialMolecule.second, {});
  }
  time_end = std::chrono::system_clock::now();

  component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Compute tail energy difference and update CPU time statistics.
  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), trialMolecule.second, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, {});
  time_end = std::chrono::system_clock::now();

  component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    // Compute polarization energy difference
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, {}, trialMolecule.second, {});
  }

  // get the total difference in energy
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   energyFourierDifference + tailEnergyDifference + polarizationDifference;

  double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
  double preFactor = system.beta * fugacity * system.simulationBox.volume /
                     double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
  double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  // Calculate acceptance probability and bias from the transition matrix.
  if (system.tmmc.doTMMC)
  {
    std::size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }
  // Check if the new macrostate exceeds the maximum allowed; reject if true.

  // apply acceptance/rejection rule
  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    component.mc_moves_statistics.addAccepted(move, 0);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMoleculePolarization(selectedComponent, trialMolecule.first, trialMolecule.second,
                                      electricFieldMoleculeNew);

    return {energyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
  };

  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
