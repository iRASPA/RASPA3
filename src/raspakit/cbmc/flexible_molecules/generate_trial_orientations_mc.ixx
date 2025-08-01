module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module cbmc_generate_trialorientations_mc;

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
  std::vector<std::pair<Atom, double>> 
    generateTrialOrientationsMonteCarloScheme(RandomNumber &random, double beta,
        const std::vector<Atom> atoms, std::size_t previousBead, std::size_t currentBead, std::vector<std::size_t> nextBeads,
            std::size_t numberOfTrialDirections);
}
