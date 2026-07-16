module;

export module mc_moves_rotation_smart_mc;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Attempts a torque-biased (rotation smart Monte Carlo) move of a single molecule.
 *
 * This move is the rotational analogue of \ref translationSmartMCMove: instead of drawing an
 * isotropic random rotation, the trial rotation vector is biased along the torque acting on the
 * selected molecule,
 * \f[
 *   \Delta\boldsymbol{\phi} = \frac{\beta \sigma^2}{2}\, \boldsymbol{\tau}_\text{old} + \sigma\, \boldsymbol{\xi},
 * \f]
 * with \f$\boldsymbol{\xi}\f$ a standard normal vector and \f$\sigma\f$ the (auto-optimized) angular
 * step size (radians). The displacement is converted to a unit quaternion via
 * \c simd_quatd::fromAxisAngle and applied with \c Component::rotate, which left-multiplies
 * the molecule's orientation quaternion and rebuilds atomic positions from the body-frame
 * template. Detailed balance is restored through a Metropolis-Hastings acceptance rule that
 * includes the ratio of the (asymmetric) forward and reverse proposal probabilities, which
 * requires the torque on the molecule both before and after the trial move.
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
std::optional<RunningEnergy> rotationSmartMCMove(RandomNumber &random, System &system,
                                                 std::size_t selectedComponent, std::size_t selectedMolecule);
}  // namespace MC_Moves
