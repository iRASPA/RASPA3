module;

export module mc_moves_rotation_smart_mc_all;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Attempts a torque-biased (rotation smart Monte Carlo) move of all molecules simultaneously.
 *
 * This is the collective variant of \ref rotationSmartMCMove: every multi-atomic molecule is rotated
 * in a single trial move, each about its own center of mass along the torque acting on it,
 * \f[
 *   \Delta\boldsymbol{\phi}_i = \frac{\beta \sigma^2}{2}\, \boldsymbol{\tau}_{\text{old},i} + \sigma\, \boldsymbol{\xi}_i,
 * \f]
 * with independent standard-normal vectors \f$\boldsymbol{\xi}_i\f$ and a shared (auto-optimized)
 * angular step size \f$\sigma\f$ (radians). Monoatomic molecules are left unchanged. The whole trial
 * configuration is accepted or rejected with a single Metropolis-Hastings criterion that includes
 * the product of the per-molecule forward/reverse proposal ratios.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 *
 * \return An optional RunningEnergy containing the new total running energy if the move is accepted;
 *         std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> rotationSmartMCMoveAll(RandomNumber &random, System &system);
}  // namespace MC_Moves
