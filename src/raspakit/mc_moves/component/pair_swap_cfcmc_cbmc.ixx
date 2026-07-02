module;

export module mc_moves_pair_swap_cfcmc_cbmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Ion-pair swap move in Continuous Fractional Component Monte Carlo with CBMC (CB/CFCMC).
 *
 * Identical to the ion-pair CFCMC swap move, except that the new fractional pair of the insertion
 * move is grown with configurational-bias Monte Carlo, and the fractional pair of the deletion move
 * is retraced with CBMC. Both fractional molecules of the pair are coupled to the same lambda (the
 * lambdaPairSwapCB histogram of the lower-index component), keeping the system charge neutral at every
 * lambda.
 *
 * Only the lower-index component of the pair performs the move (to avoid double counting).
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \return A pair containing the energy difference if accepted, and the acceptance probabilities for
 *         deletion, lambda-change, and insertion moves respectively.
 */
std::pair<std::optional<RunningEnergy>, double3> pairSwapMove_CFCMC_CBMC(RandomNumber& random, System& system,
                                                                         std::size_t selectedComponent);
}  // namespace MC_Moves
