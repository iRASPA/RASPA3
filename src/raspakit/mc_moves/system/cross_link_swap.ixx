module;

export module mc_moves_cross_link_swap;

import std;

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
/**
 * \brief Bond-swap move for cross-links: relocates one end of an existing link to another free site.
 *
 * The move conserves the number of links (Smallenburg-Sciortino bond swap). A link (i,j) and one of
 * its ends (the pivot i) are chosen uniformly; the other end is moved from j to a partner k chosen
 * uniformly among the free candidates of i (different molecule, matching bond type, within the
 * capture radius, not yet linked to i). With the link (i,j) counted as absent the candidate set of i
 * is the same in the old and the new state, so the proposal is symmetric and the move is accepted
 * with the Metropolis rule on the energy change of replacing link (i,j) by (i,k). The move is
 * rejected when j lies outside the capture radius (the reverse move could not propose it).
 *
 * \return The energy difference when accepted, std::nullopt otherwise.
 */
std::optional<RunningEnergy> crossLinkSwapMove(RandomNumber& random, System& system);
}  // namespace MC_Moves
