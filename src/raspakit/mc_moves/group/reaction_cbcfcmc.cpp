module;

module mc_moves_reaction_cbcfcmc;

import std;

import randomnumbers;
import reaction;
import mc_moves_move_types;
import running_energy;
import system;
import mc_moves_reaction_common;

std::optional<RunningEnergy> MC_Moves::reactionMove_CBCFCMC(RandomNumber& random, System& system)
{
  const Move::Types move = Move::Types::ReactionCBCFCMC;

  if (system.reactions.list.empty())
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addTrial(move);

  std::vector<Reaction*> serialReactions;
  serialReactions.reserve(system.reactions.list.size());
  for (Reaction& candidate : system.reactions.list)
  {
    if (candidate.isSerialRxCFC())
    {
      serialReactions.push_back(&candidate);
    }
  }
  if (serialReactions.empty())
  {
    return std::nullopt;
  }

  Reaction& reaction = *serialReactions[random.uniform_integer(0uz, serialReactions.size() - 1uz)];

  return ReactionCommon::serialReactionMove(random, system, reaction, move, true);
}
