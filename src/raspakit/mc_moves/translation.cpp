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

module mc_moves_translation;

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
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import move_statistics;
import mc_moves_probabilities_particles;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;

std::optional<RunningEnergy> MC_Moves::translationMove(RandomNumber &random, System &system, size_t selectedComponent,
                                                       size_t selectedMolecule,
                                                       const std::vector<Component> &components, Molecule &molecule,
                                                       std::span<Atom> molecule_atoms)
{
  double3 displacement{};
  std::chrono::system_clock::time_point time_begin, time_end;

  // Get the maximum displacement allowed for the selected component
  double3 maxDisplacement = system.components[selectedComponent].mc_moves_statistics.translationMove.maxChange;
  // Randomly select a direction (0 for x, 1 for y, 2 for z)
  size_t selectedDirection = size_t(3.0 * random.uniform());
  // Compute the displacement along the selected direction
  displacement[selectedDirection] = maxDisplacement[selectedDirection] * 2.0 * (random.uniform() - 0.5);
  // Update move statistics for the selected direction
  system.components[selectedComponent].mc_moves_statistics.translationMove.counts[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.translationMove.totalCounts[selectedDirection] += 1;

  // Copy the current electric field if polarization is computed
  if (system.forceField.computePolarization)
  {
    std::copy(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin());
  }

  // Construct the trial molecule with the new displacement
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      components[selectedComponent].translate(molecule, molecule_atoms, displacement);

  // Check if the trial molecule is inside blocked pockets
  if (system.insideBlockedPockets(system.components[selectedComponent], trialMolecule.second))
  {
    return std::nullopt;
  }

  // Compute external field energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveExternalFieldMolecule += (time_end - time_begin);
  system.mc_moves_cputime.translationMoveExternalFieldMolecule += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> frameworkMolecule;
  if (system.forceField.computePolarization)
  {
    // Get the new electric field for the molecule
    std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(selectedComponent, selectedMolecule);
    frameworkMolecule = Interactions::computeFrameworkMoleculeElectricFieldDifference(
        system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), electricFieldMoleculeNew,
        trialMolecule.second, molecule_atoms);
  }
  else
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveFrameworkMolecule += (time_end - time_begin);
  system.mc_moves_cputime.translationMoveFrameworkMolecule += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> interMolecule;
  if (system.forceField.computePolarization)
  {
    // Get the new electric fields for all molecules
    std::span<double3> electricFieldNew = system.spanOfMoleculeElectricFieldNew();
    std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(selectedComponent, selectedMolecule);
    interMolecule = Interactions::computeInterMolecularElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNew, electricFieldMoleculeNew,
        system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveMoleculeMolecule += (time_end - time_begin);
  system.mc_moves_cputime.translationMoveMoleculeMolecule += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  time_begin = std::chrono::system_clock::now();
  RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveEwald += (time_end - time_begin);
  system.mc_moves_cputime.translationMoveEwald += (time_end - time_begin);

  RunningEnergy polarization;
  if (system.forceField.computePolarization)
  {
    // Compute polarization energy difference
    polarization = Interactions::computePolarizationEnergyDifference(
        system.forceField, system.spanOfMoleculeElectricFieldNew(), system.spanOfMoleculeElectricField(),
        system.spanOfMoleculeAtoms());
  }

  // Calculate the total energy difference
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + polarization;

  // Update move construction statistics
  system.components[selectedComponent].mc_moves_statistics.translationMove.constructed[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.translationMove.totalConstructed[selectedDirection] += 1;

  // Apply acceptance/rejection rule based on Metropolis criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    // Update acceptance statistics
    system.components[selectedComponent].mc_moves_statistics.translationMove.accepted[selectedDirection] += 1;
    system.components[selectedComponent].mc_moves_statistics.translationMove.totalAccepted[selectedDirection] += 1;

    // Accept the Ewald move
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    // Update the molecule atoms with the trial positions
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    // Update the electric field if polarization is computed
    if (system.forceField.computePolarization)
    {
      std::copy(system.electricFieldNew.begin(), system.electricFieldNew.end(), system.electricField.begin());
    }

    return energyDifference;
  };
  return std::nullopt;
}
