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
std::optional<RunningEnergy> reactionMove_CBMC(RandomNumber& random, System& system);
}
