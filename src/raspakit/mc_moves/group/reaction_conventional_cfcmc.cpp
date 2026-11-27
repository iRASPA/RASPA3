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
  // the one-dimensional N chain of TMMC cannot hold this move, and the reactions change the identity of molecules
  // without the induced-dipole terms (the input reader refuses both mixes up front; a system built in code is held
  // here)
  if (system.tmmc.doTMMC || system.forceField.computePolarization) return std::nullopt;

  return ReactionCommon::parallelReactionMove(random, system, Move::Types::ReactionConventionalCFCMC, false);
}
