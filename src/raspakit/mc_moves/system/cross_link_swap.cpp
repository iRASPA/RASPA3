module;

module mc_moves_cross_link_swap;

import std;

import running_energy;
import randomnumbers;
import system;
import cross_links;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_cross_link_common;

std::optional<RunningEnergy> MC_Moves::crossLinkSwapMove(RandomNumber& random, System& system)
{
  const Move::Types move = Move::Types::CrossLinkSwap;
  system.mc_moves_statistics.addTrial(move);

  CrossLinkTable& table = system.crossLinks;
  if (table.links.empty()) return std::nullopt;

  // Choose the link and the end that stays (the pivot).
  const std::size_t linkId = static_cast<std::size_t>(random.uniform() * static_cast<double>(table.links.size()));
  const CrossLink oldLink = table.links[linkId];
  const bool pivotIsA = random.uniform() < 0.5;
  const CrossLinkSite pivot = pivotIsA ? oldLink.a : oldLink.b;
  const CrossLinkSite oldPartner = pivotIsA ? oldLink.b : oldLink.a;

  // The reverse move proposes the old partner among the candidates of the pivot: it has to be
  // within the capture radius.
  if (!CrossLinkCommon::withinCaptureRadius(system, pivot, oldPartner, oldLink.bondTypeId)) return std::nullopt;

  // Candidates of the pivot with the old link counted as absent; the old partner is then free again
  // and is part of this (state-independent) set. Exclude it from the actual choice.
  table.removeLink(linkId);
  std::vector<CrossLinkCommon::Candidate> candidates = CrossLinkCommon::partnerCandidates(system, pivot);
  std::erase_if(candidates, [&](const CrossLinkCommon::Candidate& c) { return c.site == oldPartner; });
  if (candidates.empty())
  {
    table.addLink(oldLink);
    return std::nullopt;
  }
  const CrossLinkCommon::Candidate& chosen =
      candidates[static_cast<std::size_t>(random.uniform() * static_cast<double>(candidates.size()))];

  CrossLink newLink{pivot, chosen.site, chosen.bondTypeId};

  const RunningEnergy energyDifference =
      MC_Moves::timed(system, move, Move::Timing::MoleculeMolecule,
                      [&] { return CrossLinkCommon::linkEnergy(system, newLink) - CrossLinkCommon::linkEnergy(system, oldLink); });

  system.mc_moves_statistics.addConstructed(move);

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    system.mc_moves_statistics.addAccepted(move);
    table.addLink(newLink);
    return energyDifference;
  }

  table.addLink(oldLink);
  return std::nullopt;
}
