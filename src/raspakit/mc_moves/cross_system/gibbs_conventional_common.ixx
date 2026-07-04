module;

export module mc_moves_gibbs_conventional_common;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

/**
 * \brief Shared implementation of the Gibbs conventional CFCMC particle transfer (Shi and Maginn).
 *
 * Each Gibbs box holds one active fractional molecule per swappable component, coupled such that
 * lambda_boxA + lambda_boxB = 1. A coupled lambda change in one box is mirrored in the other.
 * Integer insertions and deletions use random placement (CFCMC) or CBMC grow/retrace (CBCFCMC).
 */
export namespace MC_Moves::GibbsConventionalCommon
{
std::optional<std::pair<RunningEnergy, RunningEnergy>> gibbsConventionalMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent, bool useCBMC);
}
