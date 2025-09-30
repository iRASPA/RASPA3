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
import intra_molecular_potentials;

export namespace CBMC
{
enum class MoveType : std::size_t
{
  BondLengthChange = 0,
  BendAngleChange = 1,
  ConeChange = 2
};

std::vector<Atom> generateTrialOrientationsMonteCarloScheme(
    RandomNumber &random,  std::size_t numberOfTrialMovesPerOpenBead, double beta, 
    Component &component, const std::vector<Atom> chain_atoms,
    std::size_t previousBead, std::size_t currentBead, std::vector<std::size_t> nextBeads,
    const Potentials::IntraMolecularPotentials &intraMolecularInteractions);
}  // namespace CBMC
