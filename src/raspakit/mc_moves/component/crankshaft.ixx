module;

export module mc_moves_crankshaft;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a crankshaft move for a given flexible molecule.
 *
 * Selects a crankshaft unit uniformly from the precomputed valid units of the
 * component: a pair of anchor atoms whose removal isolates a small connected
 * segment attached only to those two anchors. The segment is rigidly rotated
 * about the axis through the anchors by a random angle. Because every rotated
 * atom keeps its distance to both (fixed) anchors, all bond lengths are
 * preserved; only the bend angles and torsions at the junctions, the
 * intramolecular non-bonded energy, and the external (framework,
 * intermolecular, Ewald) energy change. The proposal is symmetric, so the move
 * is accepted with the Metropolis criterion on the total energy difference.
 *
 * The angle is sampled with mixed step sizes: a fraction
 * component.crankshaftRandomizationFraction of the attempts randomizes the
 * angle completely (uniform in [-pi, pi], statistics channel 1), the remainder
 * perturbs it within an adaptive window (channel 0) tuned by the standard move
 * statistics machinery (Vitalis & Pappu, Methods 46 (2009)).
 *
 * The move relaxes the chain interior locally without propagating the change
 * to the chain ends (complementary to the pivot move, whose lever arm grows
 * with chain length), and it can rotate segments of flexible rings, which the
 * pivot move cannot sample. Segment sizes are capped by
 * component.crankshaftMaxSegmentSize; rigid fragments must lie entirely inside
 * the rotated segment or entirely outside it. Molecules without any valid unit
 * (fully rigid molecules, diatomics, too-small chains) and simulations with
 * polarization are not supported and lead to rejection.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> crankshaftMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                            std::size_t selectedMolecule);
}  // namespace MC_Moves
