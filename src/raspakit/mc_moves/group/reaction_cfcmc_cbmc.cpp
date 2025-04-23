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

module mc_moves_reaction_cfcmc_cbmc;

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
import cbmc_chain_data;
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
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

std::optional<RunningEnergy> MC_Moves::reactionMove_CFCMC_CBMC(
    [[maybe_unused]] RandomNumber& random, System& system,
    [[maybe_unused]] const std::vector<size_t> reactantStoichiometry,
    [[maybe_unused]] const std::vector<size_t> productStoichiometry)
{
  size_t selectedComponent = 0;
  size_t selectedMolecule = 0;

  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = system.components[selectedComponent].growType;

  [[maybe_unused]] std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

  // Fix groupId
  std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
      random, system.components[selectedComponent], system.hasExternalField, system.components, system.forceField,
      system.simulationBox, system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), system.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
      selectedComponent, selectedMolecule, 1.0, 1uz, system.numberOfTrialDirections);

  if (!growData) return std::nullopt;

  return std::nullopt;
}
