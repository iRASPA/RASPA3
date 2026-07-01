module;

export module mc_moves_gibbs_conventional_cfcmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

/**
 * \brief Gibbs conventional CFCMC particle transfer (Shi and Maginn).
 *
 * Each Gibbs box holds one active fractional molecule per swappable component, coupled such that
 * lambda_boxA + lambda_boxB = 1. A coupled lambda change in one box is mirrored in the other.
 */
export namespace MC_Moves
{
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsConventionalCFCMCMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent);

std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsConventionalCFCMCCBMCMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent);
}
