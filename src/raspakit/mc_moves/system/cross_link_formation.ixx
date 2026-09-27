module;

export module mc_moves_cross_link_formation;

import std;

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
/**
 * \brief Formation/scission move for cross-links: creates a link between two free reactive sites or
 *        removes an existing one.
 *
 * Formation (chosen with probability 1/2): a free site i is picked uniformly among the N_F free sites
 * and a partner j uniformly among its n_i candidates (different molecule, matching bond type, free,
 * within the capture radius, not yet linked to i). The pair can be generated from either end, so the
 * proposal probability is (1/N_F)(1/n_i + 1/n_j) and the reverse (scission) picks the new link among
 * N_b + 1 links. Acceptance:
 *
 *     P_acc = min[1, N_F / ((N_b + 1) (1/n_i + 1/n_j)) exp(-beta dU)]
 *
 * Scission: a link is picked uniformly among the N_b links; it is rejected when its sites are outside
 * the capture radius (the reverse move could not propose it). With N_F', n_i', n_j' the free-site and
 * candidate counts in the state without the link:
 *
 *     P_acc = min[1, N_b (1/n_i' + 1/n_j') / N_F' exp(-beta dU)]
 *
 * dU includes the bond, the junction bends, the non-bonded exclusion corrections and the formation
 * energy of the bond type, so the formation energy sets the association equilibrium.
 *
 * \return The energy difference when accepted, std::nullopt otherwise.
 */
std::optional<RunningEnergy> crossLinkFormationScissionMove(RandomNumber& random, System& system);
}  // namespace MC_Moves
