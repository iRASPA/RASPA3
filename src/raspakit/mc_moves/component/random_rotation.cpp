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

module mc_moves_random_rotation;

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
import atom;
import molecule;
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

std::optional<RunningEnergy> MC_Moves::randomRotationMove(RandomNumber &random, System &system,
                                                          size_t selectedComponent,
                                                          const std::vector<Component> &components, Molecule &molecule,
                                                          std::span<Atom> molecule_atoms)
{
  double3 angle{};
  std::chrono::system_clock::time_point time_begin, time_end;

  std::array<double3, 3> axes{double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0)};
  double3 maxAngle = system.components[selectedComponent].mc_moves_statistics.randomRotationMove.maxChange;
  size_t selectedDirection = size_t(3.0 * random.uniform());
  angle[selectedDirection] = maxAngle[selectedDirection] * 2.0 * (random.uniform() - 0.5);

  // Update move statistics for the selected direction
  system.components[selectedComponent].mc_moves_statistics.randomRotationMove.counts[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.randomRotationMove.totalCounts[selectedDirection] += 1;

  // Construct the trial positions by rotating the molecule
  double rotationAngle = angle[selectedDirection];
  double3 rotationAxis = double3(axes[selectedDirection]);
  simd_quatd q = simd_quatd::fromAxisAngle(rotationAngle, rotationAxis);
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      components[selectedComponent].rotate(molecule, molecule_atoms, q);

  // Check if the trial molecule is inside any blocked pockets
  if (system.insideBlockedPockets(system.components[selectedComponent], trialMolecule.second))
  {
    return std::nullopt;
  }

  // Compute external field energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.randomRotationMoveExternalFieldMolecule +=
      (time_end - time_begin);
  system.mc_moves_cputime.randomRotationMoveExternalFieldMolecule += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.randomRotationMoveFrameworkMolecule += (time_end - time_begin);
  system.mc_moves_cputime.randomRotationMoveFrameworkMolecule += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.randomRotationMoveMoleculeMolecule += (time_end - time_begin);
  system.mc_moves_cputime.randomRotationMoveMoleculeMolecule += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  time_begin = std::chrono::system_clock::now();
  RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.randomRotationMoveEwald += (time_end - time_begin);
  system.mc_moves_cputime.randomRotationMoveEwald += (time_end - time_begin);

  // Get the total difference in energy
  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

  // Update constructed move statistics
  system.components[selectedComponent].mc_moves_statistics.randomRotationMove.constructed[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.randomRotationMove.totalConstructed[selectedDirection] += 1;

  // Apply acceptance/rejection rule based on Metropolis criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    // Move accepted; update statistics
    system.components[selectedComponent].mc_moves_statistics.randomRotationMove.accepted[selectedDirection] += 1;
    system.components[selectedComponent].mc_moves_statistics.randomRotationMove.totalAccepted[selectedDirection] += 1;

    // Accept Ewald move and update molecule atoms
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    return energyDifference;
  };
  // Move rejected
  return std::nullopt;
}
