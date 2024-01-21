export module mc_moves_reaction_cfcmc_cbmc;

import <cstddef>;
import <optional>;
import <span>;
import <tuple>;

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
  std::optional<RunningEnergy>
  reactionMove_CFCMC_CBMC([[maybe_unused]] RandomNumber &random, System& system,
                          [[maybe_unused]] const std::vector<size_t> reactantStoichiometry,
                          [[maybe_unused]] const std::vector<size_t> productStoichiometry);
}  
