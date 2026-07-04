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

  Reaction& reaction = *serialReactions[static_cast<std::size_t>(
      random.uniform_integer(0, static_cast<int>(serialReactions.size()) - 1))];

  const std::vector<std::size_t>& activeIdsExist =
      reaction.fractionalSideIsReactants ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  if (ReactionCommon::totalStoichiometry(activeIdsExist) == 0)
  {
    return std::nullopt;
  }

  ReactionCommon::SerialMoveKind subMove = ReactionCommon::SerialMoveKind::LambdaChange;
  if (reaction.lambdaSwitchPoint > 1.0)
  {
    subMove = ReactionCommon::SerialMoveKind::LambdaChange;
  }
  else if (reaction.currentLambda < reaction.lambdaSwitchPoint)
  {
    subMove = random.uniform() < 0.5 ? ReactionCommon::SerialMoveKind::LambdaChange
                                      : ReactionCommon::SerialMoveKind::FractionalReaction;
  }
  else
  {
    subMove = random.uniform() < 0.5 ? ReactionCommon::SerialMoveKind::LambdaChange
                                      : ReactionCommon::SerialMoveKind::WholeMoleculeReaction;
  }

  switch (subMove)
  {
    case ReactionCommon::SerialMoveKind::LambdaChange:
      return ReactionCommon::serialLambdaChangeMove(random, system, reaction, move);
    case ReactionCommon::SerialMoveKind::FractionalReaction:
      return ReactionCommon::serialFractionalReactionMove(random, system, reaction, move, true);
    case ReactionCommon::SerialMoveKind::WholeMoleculeReaction:
      return ReactionCommon::serialWholeMoleculeReactionMove(random, system, reaction, move, true);
  }

  return std::nullopt;
}
