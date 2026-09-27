module;

export module mc_moves_cross_link_common;

import std;

import double3;
import running_energy;
import system;
import cross_links;

/**
 * Shared bookkeeping of the cross-link topology moves (CrossLinkSwap, CrossLinkFormationScission and
 * CrossLinkExchange).
 *
 * A "site" is a reactive atom of one molecule (see Component::reactiveSites); it is "free" when it
 * carries fewer links than its valence. Two sites are "partners" when they lie on different molecules,
 * a bond type matches their site types, they are not already linked to each other, and their separation
 * is within the capture radius of that bond type. The moves propose partners uniformly among the
 * candidates of a chosen site (free candidates for a new link, linked candidates for an exchange), so
 * the candidate counts are the proposal probabilities that enter the acceptance rules.
 */
export namespace MC_Moves::CrossLinkCommon
{
struct Candidate
{
  CrossLinkSite site;
  std::size_t bondTypeId;
};

/// The reactive site type of 'site'.
const ReactiveSite &reactiveSiteOf(const System &system, const CrossLinkSite &site);

/// The position of 'site' in the current configuration.
double3 positionOf(const System &system, const CrossLinkSite &site);

/// True when 'site' carries fewer links than its valence.
bool isFree(const System &system, const CrossLinkSite &site);

/// All sites with free capacity (integer molecules only).
std::vector<CrossLinkSite> freeSites(const System &system);

/// The free partner candidates of 'pivot' (see above); the current links of the table count as occupied.
std::vector<Candidate> partnerCandidates(const System &system, const CrossLinkSite &pivot);

/// The LINKED sites that could be a partner of 'pivot' (different molecule, matching bond type, within the
/// capture radius, carrying at least one link, not linked to 'pivot'): the exchange partners of the
/// bond-exchange move. Their number does not depend on which of them is linked to what.
std::vector<Candidate> linkedCandidates(const System &system, const CrossLinkSite &pivot);

/// The ids of the links that 'site' carries.
std::vector<std::size_t> linksOfSite(const System &system, const CrossLinkSite &site);

/// The bond type of two sites, if any.
std::optional<std::size_t> bondTypeFor(const System &system, const CrossLinkSite &a, const CrossLinkSite &b);

/// True when the sites are within the capture radius of bond type 'bondTypeId'.
bool withinCaptureRadius(const System &system, const CrossLinkSite &a, const CrossLinkSite &b,
                         std::size_t bondTypeId);

/// The energy of one (hypothetical or existing) link in the current configuration.
RunningEnergy linkEnergy(const System &system, const CrossLink &link);
}  // namespace MC_Moves::CrossLinkCommon
