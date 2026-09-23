module;

export module mc_moves_reptation;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a reptation (slithering-snake) move for a periodic chain molecule.
 *
 * The molecule must declare its repeat units ('RepeatUnits' in the molecule JSON): an
 * ordered partition of the atoms into monomer blocks (backbone plus side group) that is
 * validated at parse time to be shift-periodic (identical types, charges, connectivity,
 * potential terms and rigid fragments under the one-unit shift). The move removes the
 * repeat unit at one chain end and grows a new unit at the opposite end with CBMC,
 * translating the chain by one monomer along its own contour.
 *
 * Because of the shift periodicity the resulting molecule is chemically identical to
 * the original under a relabeling: the surviving units shift one block toward the
 * vacated end and the freshly grown unit occupies the last block, so the component's
 * canonical topology (connectivity, potentials, fragments) is reused unchanged. The
 * direction (head-to-tail or tail-to-head) is chosen with equal probability; the
 * acceptance rule is the standard configurational-bias ratio of the grow and retrace
 * Rosenbluth weights, evaluated in log space.
 *
 * The move reuses the partial-reinsertion CBMC machinery: the departing unit is
 * retraced in the current configuration, the arriving unit is grown attached to the
 * shifted chain, and both use cached growth plans. Rigid fragments and rings inside a
 * repeat unit are grown with the standard rigid-fragment/ring-closure growth.
 *
 * Because the acceptance pairs the grow weight of one chain end's plan against the
 * retrace weight of the other end's plan, the two plans must be congruent (validated
 * at parse time): the approximate CBMC trial-conformation samplers deviate slightly
 * per step structure, and these deviations only cancel across the two plans when the
 * repeat unit grows with the same step sequence from either chain end. Direction-
 * asymmetric units (e.g. a branch point that forms a multi-bead step from one side
 * only) are rejected by the 'RepeatUnits' validation.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> reptationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                           std::size_t selectedMolecule);
}  // namespace MC_Moves
