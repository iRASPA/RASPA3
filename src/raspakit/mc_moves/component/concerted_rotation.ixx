module;

export module mc_moves_concerted_rotation;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a concerted-rotation (ConRot) move on a flexible chain molecule.
 *
 * Dodd, Boone and Theodorou, Mol. Phys. 78, 961 (1993). A window of eight (nine) consecutive
 * backbone atoms a0 ... a7 (a8) is selected uniformly from the precomputed valid windows of the
 * component. The driver torsion about the a0-a1 bond is changed by a random angle, which displaces
 * a2; the trimer a3, a4, a5 is then re-bridged between the displaced a2 and the unchanged a6, a7
 * such that all bond lengths and all backbone bend angles of the window keep their values (see
 * mc_moves_concerted_rotation_geometry). Everything from a6 on, everything before a2, and every
 * atom outside the window stays where it is; side groups of a2 ... a5 are carried rigidly. The
 * move thus changes seven consecutive torsion angles at once while leaving both chain ends in
 * place, which is what makes it effective for long chains in dense phases where pivot moves fail.
 *
 * The closure has a discrete set of solutions; one of the n_new solutions is chosen at random. For
 * the reverse move the driver is rotated back and the closure solved again, which must reproduce
 * the old trimer among its n_old solutions (if it does not, a numerical failure, the move is
 * rejected). The move is a deterministic map between one-dimensional families of the torsion
 * angles; its acceptance rule carries the ratio of the closure Jacobians J = |det dP/dphi| of the
 * old and the new state (P: the fixed downstream placement, phi: the torsions about a1a2 ... a6a7)
 * and the ratio of the solution counts:
 *
 *   acc = min(1, exp(-beta dU) * (n_new / n_old) * (J_old / J_new)).
 *
 * The Jacobian is the one of the map in internal coordinates; the Cartesian measure factor
 * (products of bond lengths squared and sines of bend angles) is invariant under the move, so the
 * rule is exact for flexible (harmonic) as well as for FIXED bonds and bends.
 *
 * The driver angle is sampled with mixed step sizes: a fraction
 * component.concertedRotationRandomizationFraction of the attempts randomizes it completely
 * (uniform in [-pi, pi], statistics channel 1), the remainder perturbs it within an adaptive
 * window (channel 0).
 *
 * Valid windows: a simple path a1 ... a7 of the bond graph with a further neighbour a0 of a1,
 * where the side groups of a2 ... a5 do not connect back to the window (no ring through the
 * window), rigid fragments lie entirely in one rigidly moved group or entirely outside, and no
 * FIXED bend at a1 or a6 involves a side group together with a2 or a5. Molecules without a valid
 * window and simulations with polarization are not supported and lead to rejection.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> concertedRotationMove(RandomNumber &random, System &system,
                                                   std::size_t selectedComponent, std::size_t selectedMolecule);
}  // namespace MC_Moves
