module;

module mc_moves_rotation;

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
import move_statistics;
import mc_moves_probabilities_particles;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;

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


std::optional<RunningEnergy> 
MC_Moves::rotationMove(RandomNumber &random, System& system, size_t selectedComponent, std::span<Atom> molecule)
{
  double3 angle{};
  std::chrono::system_clock::time_point time_begin, time_end;

  std::array<double3,3> axes{double3(1.0,0.0,0.0), double3(0.0,1.0,0.0) ,double3(0.0,0.0,1.0) };
  double3 maxAngle = system.components[selectedComponent].mc_moves_statistics.rotationMove.maxChange;
  size_t selectedDirection = size_t(3.0 * random.uniform());
  angle[selectedDirection] = maxAngle[selectedDirection] * 2.0 * (random.uniform() - 0.5);

  system.components[selectedComponent].mc_moves_statistics.rotationMove.counts[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.rotationMove.totalCounts[selectedDirection] += 1;

  // construct the trial positions
  size_t startingBead = system.components[selectedComponent].startingBead;
  std::vector<Atom> trialMolecule(molecule.size());
  double rotationAngle = angle[selectedDirection];
  double3 rotationAxis = double3(axes[selectedDirection]);
  double3x3 rotationMatrix = double3x3(simd_quatd::fromAxisAngle(rotationAngle, rotationAxis));
  std::transform(molecule.begin(), molecule.end(), trialMolecule.begin(),
          [&](Atom a) { a.position = rotationMatrix * (a.position - molecule[startingBead].position) 
                        + molecule[startingBead].position; return a; });

  // compute external field energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule =
    Interactions::computeExternalFieldEnergyDifference(system.hasExternalField, system.forceField, system.simulationBox,
                                                       trialMolecule, molecule);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.rotationMoveExternalFieldMolecule += (time_end - time_begin);
  system.mc_moves_cputime.rotationMoveExternalFieldMolecule += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // compute framework-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> frameworkMolecule = 
    Interactions::computeFrameworkMoleculeEnergyDifference(system.forceField, system.simulationBox,                     
                                                           system.spanOfFrameworkAtoms(), trialMolecule, molecule);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.rotationMoveFrameworkMolecule += (time_end - time_begin);
  system.mc_moves_cputime.rotationMoveFrameworkMolecule += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // compute molecule-molecule energy contribution
  time_begin = std::chrono::system_clock::now();
  std::optional<RunningEnergy> interMolecule = 
    Interactions::computeInterMolecularEnergyDifference(system.forceField, system.simulationBox,                     
                                                        system.spanOfMoleculeAtoms(), trialMolecule, molecule);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.rotationMoveMoleculeMolecule += (time_end - time_begin);
  system.mc_moves_cputime.rotationMoveMoleculeMolecule += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // compute Ewald energy contribution
  time_begin = std::chrono::system_clock::now();
  RunningEnergy ewaldFourierEnergy = 
    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                               system.storedEik, system.totalEik,
                                               system.forceField, system.simulationBox,
                                               trialMolecule, molecule);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.rotationMoveEwald += (time_end - time_begin);
  system.mc_moves_cputime.rotationMoveEwald += (time_end - time_begin);

  // get the total difference in energy
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() +
                                   interMolecule.value() + ewaldFourierEnergy;

  system.components[selectedComponent].mc_moves_statistics.rotationMove.constructed[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.rotationMove.totalConstructed[selectedDirection] += 1;

  // apply acceptance/rejection rule
  if (random.uniform() < std::exp(-system.beta * energyDifference.total()))
  {
    system.components[selectedComponent].mc_moves_statistics.rotationMove.accepted[selectedDirection] += 1;
    system.components[selectedComponent].mc_moves_statistics.rotationMove.totalAccepted[selectedDirection] += 1;

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    std::copy(trialMolecule.cbegin(), trialMolecule.cend(), molecule.begin());

    return energyDifference;
  };
  return std::nullopt;
}
