module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module cbmc_generate_trialorientations_sphere;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import atom;
import simulationbox;
import cbmc_chain_data;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import framework;
import component;
import interpolation_energy_grid;
import bond_potential;

export namespace CBMC
{
  std::vector<std::vector<Atom>> 
    generateTrialOrientationsSimpleSphere(RandomNumber &random, double beta, const Atom &first_bead,
                                          std::size_t numberOfTrialDirections, const BondPotential &bond);
}
