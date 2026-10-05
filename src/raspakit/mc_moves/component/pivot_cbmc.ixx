module;

export module mc_moves_pivot_cbmc;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a configurational-bias pivot move for a flexible molecule.
 *
 * The geometry is that of the plain pivot move (mc_moves_pivot): a pivot bond is selected
 * uniformly from the precomputed valid bonds of the component and the smaller of the two chain
 * parts hanging off it is rotated rigidly about the bond axis, which preserves all bond lengths
 * and bend angles. Instead of a single random angle, k trial angles are drawn and one of them is
 * selected with probability proportional to its Boltzmann factor (the Rosenbluth selection of
 * configurational-bias Monte Carlo, applied to the single torsional degree of freedom of the
 * pivot; Frenkel and Smit, "Understanding Molecular Simulation", orientational bias). In dense
 * phases, where almost every random pivot angle sweeps the rotated part through neighbouring
 * molecules, the biased selection finds the few angles that fit and raises the acceptance of
 * large rotations by up to an order of magnitude at k times the cost of a plain pivot.
 *
 * Proposal. The k trial angles phi_i are drawn uniformly from [-maxAngle, maxAngle] (channel 0,
 * adaptive window) or from [-pi, pi] (channel 1, full randomization; a fraction
 * component.pivotCBMCRandomizationFraction of the attempts). For each trial the "bias energy"
 * difference dU_i with respect to the current configuration is evaluated: external field,
 * framework, intermolecular (van der Waals and real-space Coulomb), cross-link and the complete
 * intramolecular energy (the torsions through the pivot bond and the intramolecular non-bonded
 * terms are the ones that change). Trial i is selected with probability w_i / W_new,
 * w_i = exp(-beta dU_i), W_new = sum_i w_i; an overlap or a blocked pocket gives w_i = 0. The
 * reverse move draws its k - 1 other trial angles about the NEW configuration, so that the old
 * configuration is one of its trials: W_old = 1 + sum_{j=1}^{k-1} exp(-beta dU_j) with the
 * trial angles phi_selected + psi_j, psi_j uniform in the same window. The move is accepted with
 *
 *   acc = min(1, (W_new / W_old) * exp(-beta dU_Fourier)),
 *
 * where dU_Fourier is the Ewald Fourier-space difference of the selected trial, which is too
 * costly to include in the bias and enters as a correction (as in all CBMC moves of this code).
 * The uniform angle proposals are symmetric and cancel. With k = 1 the move reduces to the plain
 * pivot.
 *
 * Molecules without a valid pivot bond (fully rigid molecules, pure rings, diatomics) and
 * simulations with polarization are not supported and lead to rejection.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> pivotCBMCMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                           std::size_t selectedMolecule);
}  // namespace MC_Moves
