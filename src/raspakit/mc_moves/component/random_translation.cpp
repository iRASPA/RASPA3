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

module mc_moves_random_translation;

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
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::randomTranslationMove(RandomNumber &random, System &system,
                                                             size_t selectedComponent,
                                                             const std::vector<Component> &components,
                                                             Molecule &molecule, std::span<Atom> molecule_atoms)
{
  double3 s, displacement{};
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::RandomTranslation;
  Component &component = system.components[selectedComponent];

  // Select a random direction (0: x, 1: y, 2: z)
  size_t selectedDirection = size_t(3.0 * random.uniform());
  // Generate a random displacement along the selected direction
  s[selectedDirection] = random.uniform();
  displacement = system.simulationBox.cell * s;

  // Update move counts for the selected direction
  component.mc_moves_statistics.addTrial(move, selectedDirection);

  // Construct the trial positions
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      components[selectedComponent].translate(molecule, molecule_atoms, displacement);

  // Reject move if trial positions are inside blocked pockets
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  // Compute external field energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["ExternalField-Molecule"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["ExternalField-Molecule"] += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["Molecule-Molecule"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Molecule-Molecule"] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  time_begin = std::chrono::system_clock::now();
  RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Get the total difference in energy
  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

  // Update constructed move counts for the selected direction
  component.mc_moves_statistics.addConstructed(move, selectedDirection);

  // Apply Metropolis acceptance criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    // Move accepted, update accepted move counts
    component.mc_moves_statistics.addAccepted(move, selectedDirection);

    // Accept Ewald move and update system state
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    return energyDifference;
  };
  // Move rejected
  return std::nullopt;
}
