module;

export module mc_moves_cross_link_exchange;

import std;

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
/**
 * \brief Bond-exchange move for cross-links: two links trade partners.
 *
 *   (a,b) + (c,d)  ->  (a,c) + (b,d)
 *
 * No site is freed or occupied, so the number of links AND the link count of every site are
 * conserved: this is the metathesis / transesterification-type exchange of vitrimers and the only
 * topology move that still rewires a network whose sites are all saturated (where the bond swap has
 * no free partner to relocate to). A link (a,b) and its pivot end a are chosen uniformly; the exchange
 * partner c is chosen uniformly among the LINKED candidates of a (different molecule, matching bond
 * type, within the capture radius, not linked to a), and the link (c,d) that c gives up uniformly
 * among the links of c. The reverse route from the new state (link (a,c), pivot a, candidate b, its
 * link (b,d)) sees a candidate set of the same size (b takes the place of c), so the proposal ratio
 * reduces to n_links(c) / n_links(b) and the move is accepted with
 *
 *   min(1, exp(-beta [u(a,c) + u(b,d) - u(a,b) - u(c,d)]) * n_links(c) / n_links(b)),
 *
 * where n_links is the number of links a site carries (1 for valence-1 sites: a plain Metropolis rule).
 * The move is rejected when b lies outside the capture radius of a (the reverse route could not
 * propose it), when the new link (b,d) is not admissible (same molecule, no bond type, already linked,
 * or b == d), and when fewer than two links exist.
 *
 * \return The energy difference when accepted, std::nullopt otherwise.
 */
std::optional<RunningEnergy> crossLinkExchangeMove(RandomNumber& random, System& system);
}  // namespace MC_Moves
