module;

export module mc_moves_reaction;

import std;

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
std::optional<RunningEnergy> reactionMove([[maybe_unused]] RandomNumber& random, System& system,
                                          [[maybe_unused]] const std::vector<std::size_t> reactantStoichiometry,
                                          [[maybe_unused]] const std::vector<std::size_t> productStoichiometry);
}
