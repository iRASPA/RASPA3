module;

export module mc_moves_gibbs_conventional_cbcfcmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

/**
 * \brief Gibbs conventional CB/CFCMC particle transfer (Shi and Maginn).
 *
 * Same sub-move structure as GibbsConventionalCFCMCMove, but integer insertions and deletions use
 * CBMC grow and retrace.
 */
export namespace MC_Moves
{
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsConventionalCBCFCMCMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent);
}
