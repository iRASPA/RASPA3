module;

export module mc_moves_double_bridging;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a (symmetric) double-bridging move between two chains of the same component.
 *
 * Double bridging (Karayiannis, Mavrantzas and Theodorou, Phys. Rev. Lett. 88, 105503 (2002);
 * J. Chem. Phys. 117, 5465 (2002)) is a connectivity-altering move for dense chain systems: a
 * trimer is excised from each of two neighbouring chains at the same backbone site, and the two
 * chain tails are exchanged by bridging the head of each chain onto the tail of the other with a
 * new trimer. Using the same site on both chains keeps every chain length unchanged, so the move
 * samples the same monodisperse ensemble as the other moves while changing the large-scale
 * conformations of two chains at once, which decorrelates the end-to-end vectors orders of
 * magnitude faster than local moves.
 *
 * The molecule must be acyclic with a backbone of at least seven units (see
 * Component::BridgingTopology). Side groups of the tails travel with their backbone atoms, side
 * groups of the re-bridged trimers are carried rigidly with the local frames of their backbone
 * atoms. Since the component's canonical topology is shared by all its molecules, the exchange is
 * a copy of positions into the slots of the other chain.
 *
 * Proposal. The site is chosen uniformly among the valid ones; the partner chain uniformly among
 * the chains whose anchors are within bridging reach of the selected chain (a symmetric, minimum-
 * image distance criterion on the fixed anchor atoms only, so that the candidate counts of the new
 * state can be evaluated exactly). Each of the two closures picks one of its discrete solutions
 * uniformly. The acceptance rule is the Metropolis factor times the closure factors of both
 * bridges (solution-count ratio and closure-Jacobian ratio, see mc_moves_bridging_common) and the
 * ratio of the partner-selection probabilities summed over the two ways to initiate the move.
 *
 * Only whole (non-fractional) molecules take part. Polarization is not supported.
 *
 * \param random Random number generator instance.
 * \param system The simulation system.
 * \param selectedComponent Index of the component.
 * \param selectedMolecule Index of the first chain within the component (must be an integer molecule).
 *
 * \return Energy difference if the move is accepted; \c std::nullopt otherwise.
 */
std::optional<RunningEnergy> doubleBridgingMove(RandomNumber& random, System& system, std::size_t selectedComponent,
                                                std::size_t selectedMolecule);
}  // namespace MC_Moves
