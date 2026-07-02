module;

module mc_moves_reaction_conventional_cfcmc_cbmc;

import std;

import randomnumbers;
import mc_moves_move_types;
import running_energy;
import system;
import mc_moves_reaction_common;

// Parallel Rx/CFC reaction move following Rosch and Maginn, J. Chem. Theory Comput. 7 (2011) 269-279:
// lambda changes within [0, 1] (eq. 26) and in-place boundary-crossing reactions (eq. 27). The shared
// implementation uses CBMC growth referenced to the ideal-gas Rosenbluth weights and the staged
// coupling schedule of scaling.ixx.
std::optional<RunningEnergy> MC_Moves::reactionMove_ConventionalCFCMCCBMC(RandomNumber& random, System& system)
{
  return ReactionCommon::parallelReactionMove(random, system, Move::Types::ReactionConventionalCFCMCCBMC);
}
