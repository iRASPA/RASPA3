module;

export module mc_moves_translation_smart_mc_all;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Attempts a force-biased (translation smart Monte Carlo) move of all molecules simultaneously.
 *
 * This is the collective variant of \ref translationSmartMCMove: every molecule is displaced in a
 * single trial move, each along the force acting on it,
 * \f[
 *   \Delta r_i = \frac{\beta \sigma^2}{2}\, F_{\text{old},i} + \sigma\, \xi_i,
 * \f]
 * with independent standard-normal vectors \f$\xi_i\f$ and a shared (auto-optimized) step size
 * \f$\sigma\f$. The whole trial configuration is accepted or rejected with a single Metropolis-Hastings
 * criterion that includes the product of the per-molecule forward/reverse proposal ratios, which
 * requires the force on every molecule before and after the collective displacement.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 *
 * \return An optional RunningEnergy containing the new total running energy if the move is accepted;
 *         std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> translationSmartMCMoveAll(RandomNumber &random, System &system);
}  // namespace MC_Moves
