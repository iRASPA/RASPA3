module;

export module mc_moves_pair_swap_cfcmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Ion-pair swap move in Continuous Fractional Component Monte Carlo (CFCMC).
 *
 * The pair of linked components (selectedComponent and its pairComponentId) each hold one fractional
 * molecule; both fractional molecules are coupled to the same lambda (the lambdaGC histogram of the
 * lower-index component). Because the coulombic scaling of both fractional molecules is identical,
 * the system remains charge neutral at every lambda.
 *
 * Depending on the randomly selected new lambda bin the move is:
 * - **Insertion move** (lambda crosses 1): the fractional pair is made integer, and a new fractional
 *   pair is inserted at random positions with lambda_new = epsilon.
 * - **Deletion move** (lambda crosses 0): the fractional pair is removed, and a randomly selected
 *   integer molecule of each component becomes the new fractional pair with lambda_new = 1 - epsilon.
 * - **Lambda-change move**: the lambda of both fractional molecules is changed simultaneously.
 *
 * Only the lower-index component of the pair performs the move (to avoid double counting).
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \return A pair containing the energy difference if accepted, and the acceptance probabilities for
 *         deletion, lambda-change, and insertion moves respectively.
 */
std::pair<std::optional<RunningEnergy>, double3> pairSwapMove_CFCMC(RandomNumber& random, System& system,
                                                                    std::size_t selectedComponent);
}  // namespace MC_Moves
