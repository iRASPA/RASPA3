module;

export module mc_moves_gibbs_swap_cbcfcmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

/**
 * \brief Performs a Gibbs swap move in a CFCMC/CBMC simulation.
 *
 * Same sub-move structure as GibbsSwapMove_CFCMC, but integer insertions and deletions use CBMC grow and retrace.
 */
export namespace MC_Moves
{
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMove_CBCFCMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent,
    [[maybe_unused]] std::size_t& fractionalMoleculeSystem);
}
