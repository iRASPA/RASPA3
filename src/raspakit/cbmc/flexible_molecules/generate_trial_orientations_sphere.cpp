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
    
module cbmc_generate_trialorientations_sphere;
    
#ifndef USE_LEGACY_HEADERS
import std;
#endif
  
import randomnumbers;
import component;
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


std::vector<std::vector<Atom>>
  CBMC::generateTrialOrientationsSimpleSphere(RandomNumber &random, double beta, const Atom &currentBead,
                                              std::size_t numberOfTrialDirections, const BondPotential &bond)
{
  std::vector<std::vector<Atom>> trial_positions(numberOfTrialDirections);

  for(std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    double bond_length = bond.generateBondLength(random, beta);
    double3 unit_vector = random.randomVectorOnUnitSphere();

    Atom trial_atom = currentBead;
    trial_atom.position = currentBead.position + bond_length * unit_vector;

    trial_positions[i] = { trial_atom };
  }
  return trial_positions;
}
