export module mc_moves_volume;

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
  std::optional<RunningEnergy> volumeMove(RandomNumber &random, System &system);
}
