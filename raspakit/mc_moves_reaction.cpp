module;

module mc_moves;

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

// mc_moves_reaction.cpp

std::optional<RunningEnergy> MC_Moves::reactionMove([[maybe_unused]] RandomNumber &random, System& system, [[maybe_unused]] const std::vector<size_t> reactantStoichiometry, [[maybe_unused]] const std::vector<size_t> productStoichiometry) const
{
  size_t selectedComponent = 0;
  size_t selectedMolecule = 0;

  double cutOffVDW = system.forceField.cutOffVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;

  Component::GrowType growType = system.components[selectedComponent].growType;

  [[maybe_unused]] std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  std::vector<Atom> atoms = system.components[selectedComponent].newAtoms(1.0, system.numberOfMoleculesPerComponent[selectedComponent]);
  std::optional<ChainData> growData = system.growMoleculeSwapInsertion(random, growType, cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, 1.0, atoms);
  if (!growData) return std::nullopt;


  return std::nullopt;
}
