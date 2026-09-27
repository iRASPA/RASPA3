module;

module mc_moves_cross_link_exchange;

import std;

import running_energy;
import randomnumbers;
import system;
import cross_links;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_cross_link_common;

std::optional<RunningEnergy> MC_Moves::crossLinkExchangeMove(RandomNumber& random, System& system)
{
  const Move::Types move = Move::Types::CrossLinkExchange;
  system.mc_moves_statistics.addTrial(move);

  CrossLinkTable& table = system.crossLinks;
  if (table.links.size() < 2) return std::nullopt;

  // The first link (a,b) and its pivot end a.
  const std::size_t firstId = static_cast<std::size_t>(random.uniform() * static_cast<double>(table.links.size()));
  const CrossLink first = table.links[firstId];
  const bool pivotIsA = random.uniform() < 0.5;
  const CrossLinkSite a = pivotIsA ? first.a : first.b;
  const CrossLinkSite b = pivotIsA ? first.b : first.a;

  // The reverse route proposes b as the exchange partner of a: b has to be within the capture radius.
  if (!CrossLinkCommon::withinCaptureRadius(system, a, b, first.bondTypeId)) return std::nullopt;

  // The exchange partner c: a linked site near a that is not linked to a. (This set has the same size
  // in the new state, where c is linked to a and b is not.)
  const std::vector<CrossLinkCommon::Candidate> candidates = CrossLinkCommon::linkedCandidates(system, a);
  if (candidates.empty()) return std::nullopt;
  const CrossLinkCommon::Candidate& chosen =
      candidates[static_cast<std::size_t>(random.uniform() * static_cast<double>(candidates.size()))];
  const CrossLinkSite c = chosen.site;

  // The link (c,d) that c gives up, uniformly among the links of c.
  const std::vector<std::size_t> linksOfC = CrossLinkCommon::linksOfSite(system, c);
  const std::size_t secondId = linksOfC[static_cast<std::size_t>(random.uniform() * static_cast<double>(linksOfC.size()))];
  const CrossLink second = table.links[secondId];
  const CrossLinkSite d = (second.a == c) ? second.b : second.a;

  // The new link (b,d) has to be admissible: different molecules, distinct sites, a bond type, and not
  // already present (b and d each keep their link count, so the valences are automatically respected).
  if (b == d || b.sameMolecule(d)) return std::nullopt;
  if (table.isLinked(b, d)) return std::nullopt;
  const std::optional<std::size_t> bondTypeBD = CrossLinkCommon::bondTypeFor(system, b, d);
  if (!bondTypeBD.has_value()) return std::nullopt;

  const CrossLink newFirst{a, c, chosen.bondTypeId};
  const CrossLink newSecond{b, d, bondTypeBD.value()};

  const RunningEnergy energyDifference =
      MC_Moves::timed(system, move, Move::Timing::MoleculeMolecule,
                      [&]
                      {
                        return CrossLinkCommon::linkEnergy(system, newFirst) +
                               CrossLinkCommon::linkEnergy(system, newSecond) -
                               CrossLinkCommon::linkEnergy(system, first) - CrossLinkCommon::linkEnergy(system, second);
                      });

  system.mc_moves_statistics.addConstructed(move);

  // Proposal ratio: the reverse route draws the link of b among n_links(b) links where the forward route
  // drew among the n_links(c) links of c (both counts are unchanged by the exchange).
  const double proposalRatio =
      static_cast<double>(linksOfC.size()) / static_cast<double>(table.linkCount(b));

  if (random.uniform() < proposalRatio * std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    system.mc_moves_statistics.addAccepted(move);
    // Remove the higher id first: removing a link moves the last link into its slot.
    table.removeLink(std::max(firstId, secondId));
    table.removeLink(std::min(firstId, secondId));
    table.addLink(newFirst);
    table.addLink(newSecond);
    return energyDifference;
  }

  return std::nullopt;
}
