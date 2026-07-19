module;

export module mc_moves_translation_rotation_smart_mc;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Attempts a combined force- and torque-biased (smart Monte Carlo) move of a single molecule.
 *
 * This move combines \ref translationSmartMCMove and \ref rotationSmartMCMove into a single trial
 * move: the center-of-mass displacement is biased along the force and the rotation vector is biased
 * along the torque acting on the selected molecule,
 * \f[
 *   \Delta r = \frac{\beta \sigma_t^2}{2}\, F_\text{old} + \sigma_t\, \xi_t, \qquad
 *   \Delta\boldsymbol{\phi} = \frac{\beta \sigma_r^2}{2}\, \boldsymbol{\tau}_\text{old} + \sigma_r\, \xi_r,
 * \f]
 * with \f$\xi_t, \xi_r\f$ independent standard normal vectors and \f$\sigma_t\f$ (Angstrom) and
 * \f$\sigma_r\f$ (radians) the two (auto-optimized) step sizes. This is equivalent to one Euler step
 * of overdamped rigid-body Brownian dynamics, accepted with a Metropolis-Hastings criterion.
 *
 * The translation is applied first and the rotation second; because the rotation pivots about a
 * material point of the molecule (the center of mass for rigid molecules, the starting bead for
 * flexible ones) the two operations commute and the reverse proposal \f$(-\Delta r, -\Delta\phi)\f$
 * maps the trial state exactly back onto the old state. Since the translational and rotational
 * noises are independent, the Hastings correction is simply the sum of the translational and
 * rotational log-bias terms, each requiring the force respectively torque both before and after the
 * trial move. Both are obtained from a single per-atom gradient evaluation per endpoint, making
 * this move cheaper than performing the two separate smart MC moves.
 *
 * Monoatomic molecules (no rotational degrees of freedom) are rejected immediately.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule within the component.
 *
 * \return An optional RunningEnergy containing the energy difference if the move is accepted;
 *         std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> translationRotationSmartMCMove(RandomNumber &random, System &system,
                                                            std::size_t selectedComponent,
                                                            std::size_t selectedMolecule);
}  // namespace MC_Moves
