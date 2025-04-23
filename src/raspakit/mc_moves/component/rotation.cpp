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
#include <vector>
#endif

module mc_moves_rotation;

#ifndef USE_LEGACY_HEADERS
import <complex>;
import <vector>;
import <array>;
import <optional>;
import <span>;
import <optional>;
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
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::rotationMove(RandomNumber &random, System &system, size_t selectedComponent,
                                                    const std::vector<Component> &components, Molecule &molecule,
                                                    std::span<Atom> molecule_atoms)
{
  double3 angle{};
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::Rotation;
  Component &component = system.components[selectedComponent];

  std::array<double3, 3> axes{double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0)};
  size_t selectedDirection = size_t(3.0 * random.uniform());

  double maxAngle = component.mc_moves_statistics.getMaxChange(move, selectedDirection);

  angle[selectedDirection] = maxAngle * 2.0 * (random.uniform() - 0.5);

  component.mc_moves_statistics.addTrial(move, selectedDirection);

  // construct the trial positions
  double rotationAngle = angle[selectedDirection];
  double3 rotationAxis = double3(axes[selectedDirection]);
  simd_quatd q = simd_quatd::fromAxisAngle(rotationAngle, rotationAxis);
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      components[selectedComponent].rotate(molecule, molecule_atoms, q);

  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  // compute external field energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["ExternalField-Molecule"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["ExternalField-Molecule"] += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // compute framework-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // compute molecule-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["Molecule-Molecule"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Molecule-Molecule"] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // compute Ewald energy contribution
  time_begin = std::chrono::system_clock::now();
  RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // get the total difference in energy
  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

  component.mc_moves_statistics.addConstructed(move, selectedDirection);

  // apply acceptance/rejection rule
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, selectedDirection);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    return energyDifference;
  };
  return std::nullopt;
}
