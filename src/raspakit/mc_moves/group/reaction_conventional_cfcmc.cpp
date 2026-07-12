module;

module mc_moves_reaction_conventional_cfcmc;

import std;

import randomnumbers;
import mc_moves_move_types;
import running_energy;
import system;
import mc_moves_reaction_common;

// Parallel Rx/CFC reaction move following Rosch and Maginn, J. Chem. Theory Comput. 7 (2011) 269-279:
// lambda changes within [0, 1] (eq. 26) and in-place boundary-crossing reactions (eq. 27). The shared
// conventional implementation uses one equilibrated random molecular configuration per insertion
// and a complete Metropolis energy difference, without Rosenbluth weights.
std::optional<RunningEnergy> MC_Moves::reactionMove_ConventionalCFCMC(RandomNumber& random, System& system)
{
  return ReactionCommon::parallelReactionMove(random, system, Move::Types::ReactionConventionalCFCMC, false);
}
