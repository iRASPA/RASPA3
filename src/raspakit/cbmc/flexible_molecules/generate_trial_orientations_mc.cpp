module;
      
#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric> 
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif
    
module cbmc_generate_trialorientations_mc;
    
#ifndef USE_LEGACY_HEADERS
import std;
#endif
  
import randomnumbers;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
import double3x3;
import simulationbox; 
import energy_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;   
import cbmc_multiple_first_bead;
import cbmc_interactions;
import forcefield;  
import energy_factor;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;
import bond_potential;


std::vector<std::pair<Atom, double>> CBMC::generateTrialOrientationsMonteCarloScheme(RandomNumber &random, double beta,
            const std::vector<Atom> atoms, std::size_t previousBead, std::size_t currentBead, std::vector<std::size_t> nextBeads,
            std::size_t numberOfTrialDirections)
{
  std::vector<std::pair<Atom, double>> trial_positions(numberOfTrialDirections);

  double3 last_bond_vector = (atoms[previousBead].position - atoms[currentBead].position).normalized();

  return trial_positions;
}
