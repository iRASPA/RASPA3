module;

module mc_moves_cross_link_formation;

import std;

import running_energy;
import randomnumbers;
import system;
import cross_links;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_cross_link_common;

namespace
{
constexpr std::size_t formationDirection = 0;
constexpr std::size_t scissionDirection = 1;

// Number of partner candidates of 'site'; zero when the site is not free.
std::size_t candidateCount(const System& system, const CrossLinkSite& site)
{
  return MC_Moves::CrossLinkCommon::partnerCandidates(system, site).size();
}
}  // namespace

std::optional<RunningEnergy> MC_Moves::crossLinkFormationScissionMove(RandomNumber& random, System& system)
{
  const Move::Types move = Move::Types::CrossLinkFormationScission;
  CrossLinkTable& table = system.crossLinks;
  if (!table.enabled()) return std::nullopt;

  if (random.uniform() < 0.5)
  {
    // Formation.
    system.mc_moves_statistics.addTrial(move, formationDirection);

    const std::vector<CrossLinkSite> free = CrossLinkCommon::freeSites(system);
    if (free.empty()) return std::nullopt;
    const double numberOfFreeSites = static_cast<double>(free.size());

    const CrossLinkSite i = free[static_cast<std::size_t>(random.uniform() * numberOfFreeSites)];
    const std::vector<CrossLinkCommon::Candidate> candidates = CrossLinkCommon::partnerCandidates(system, i);
    if (candidates.empty()) return std::nullopt;
    const CrossLinkCommon::Candidate& chosen =
        candidates[static_cast<std::size_t>(random.uniform() * static_cast<double>(candidates.size()))];
    const CrossLinkSite j = chosen.site;

    // The pair could also have been generated from j; the candidates of j include i.
    const double n_i = static_cast<double>(candidates.size());
    const double n_j = static_cast<double>(candidateCount(system, j));

    const CrossLink newLink{i, j, chosen.bondTypeId};
    const RunningEnergy energyDifference = MC_Moves::timed(system, move, Move::Timing::MoleculeMolecule,
                                                           [&] { return CrossLinkCommon::linkEnergy(system, newLink); });

    system.mc_moves_statistics.addConstructed(move, formationDirection);

    const double numberOfLinksAfter = static_cast<double>(table.links.size() + 1);
    const double proposalRatio = numberOfFreeSites / (numberOfLinksAfter * (1.0 / n_i + 1.0 / n_j));

    if (random.uniform() < proposalRatio * std::exp(-system.beta * energyDifference.potentialEnergy()))
    {
      system.mc_moves_statistics.addAccepted(move, formationDirection);
      table.addLink(newLink);
      return energyDifference;
    }
    return std::nullopt;
  }

  // Scission.
  system.mc_moves_statistics.addTrial(move, scissionDirection);
  if (table.links.empty()) return std::nullopt;
  const double numberOfLinksBefore = static_cast<double>(table.links.size());

  const std::size_t linkId = static_cast<std::size_t>(random.uniform() * numberOfLinksBefore);
  const CrossLink oldLink = table.links[linkId];

  // The reverse (formation) move can only propose the pair when it is within the capture radius.
  if (!CrossLinkCommon::withinCaptureRadius(system, oldLink.a, oldLink.b, oldLink.bondTypeId)) return std::nullopt;

  const RunningEnergy energyDifference = MC_Moves::timed(
      system, move, Move::Timing::MoleculeMolecule, [&] { return -CrossLinkCommon::linkEnergy(system, oldLink); });

  // Proposal counts in the state without the link.
  table.removeLink(linkId);
  const double numberOfFreeSites = static_cast<double>(CrossLinkCommon::freeSites(system).size());
  const double n_i = static_cast<double>(candidateCount(system, oldLink.a));
  const double n_j = static_cast<double>(candidateCount(system, oldLink.b));

  system.mc_moves_statistics.addConstructed(move, scissionDirection);

  const double proposalRatio = numberOfLinksBefore * (1.0 / n_i + 1.0 / n_j) / numberOfFreeSites;

  if (random.uniform() < proposalRatio * std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    system.mc_moves_statistics.addAccepted(move, scissionDirection);
    return energyDifference;
  }

  table.addLink(oldLink);
  return std::nullopt;
}
