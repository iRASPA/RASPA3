module;

export module mc_moves_force_biased_translation;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Attempts a force-biased (smart Monte Carlo) translation of a single molecule.
 *
 * This move is the force-biased analogue of the ordinary translation move: instead of drawing an
 * isotropic random displacement, the trial center-of-mass displacement is biased along the force
 * acting on the selected molecule,
 * \f[
 *   \Delta r = \frac{\beta \sigma^2}{2}\, F_\text{old} + \sigma\, \xi,
 * \f]
 * with \f$\xi\f$ a standard normal vector and \f$\sigma\f$ the (auto-optimized) step size. Detailed
 * balance is restored through a Metropolis-Hastings acceptance rule that includes the ratio of the
 * (asymmetric) forward and reverse proposal probabilities, which requires the force on the molecule
 * both before and after the trial move.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule within the component.
 *
 * \return An optional RunningEnergy containing the energy difference if the move is accepted;
 *         std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> forceBiasTranslationMove(RandomNumber &random, System &system,
                                                      std::size_t selectedComponent, std::size_t selectedMolecule);
}  // namespace MC_Moves
