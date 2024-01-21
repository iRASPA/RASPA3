module;

module mc_moves_translation;

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
MC_Moves::translationMove(RandomNumber &random, System &system, size_t selectedComponent, 
                          std::span<Atom> molecule)
{
  double3 displacement{};
  double3 maxDisplacement = system.components[selectedComponent].mc_moves_statistics.translationMove.maxChange;
  size_t selectedDirection = size_t(3.0 * random.uniform());
  displacement[selectedDirection] = maxDisplacement[selectedDirection] * 2.0 * (random.uniform() - 0.5);
  system.components[selectedComponent].mc_moves_statistics.translationMove.counts[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.translationMove.totalCounts[selectedDirection] += 1;

  std::vector<Atom> trialPositions(molecule.size());
  std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
      [&](Atom a) { a.position += displacement; return a; });
  std::span<Atom> newMolecule{trialPositions.begin(), trialPositions.end()};

  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  std::optional<RunningEnergy> frameworkMolecule = 
    Interactions::computeFrameworkMoleculeEnergyDifference(system.forceField, system.simulationBox,
                                                           system.spanOfFrameworkAtoms(), newMolecule, molecule);
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveNonEwald += (t2 - t1);
  system.mc_moves_cputime.translationMoveNonEwald += (t2 - t1);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
  std::optional<RunningEnergy> interMolecule = 
    Interactions::computeInterMolecularEnergyDifference(system.forceField, system.simulationBox,                     
                                                        system.spanOfMoleculeAtoms(), newMolecule, molecule);
  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveNonEwald += (u2 - u1);
  system.mc_moves_cputime.translationMoveNonEwald += (u2 - u1);
  if (!interMolecule.has_value()) return std::nullopt;

  std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
  RunningEnergy ewaldFourierEnergy = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
  std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.translationMoveEwald += (v2 - v1);
  system.mc_moves_cputime.translationMoveEwald += (v2 - v1);

  RunningEnergy energyDifference = frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

  system.components[selectedComponent].mc_moves_statistics.translationMove.constructed[selectedDirection] += 1;
  system.components[selectedComponent].mc_moves_statistics.translationMove.totalConstructed[selectedDirection] += 1;

  if (random.uniform() < std::exp(-system.beta * energyDifference.total()))
  {
    system.components[selectedComponent].mc_moves_statistics.translationMove.accepted[selectedDirection] += 1;
    system.components[selectedComponent].mc_moves_statistics.translationMove.totalAccepted[selectedDirection] += 1;

    system.acceptEwaldMove();
    std::copy(trialPositions.cbegin(), trialPositions.cend(), molecule.begin());

    return energyDifference;
  };
  return std::nullopt;
}
