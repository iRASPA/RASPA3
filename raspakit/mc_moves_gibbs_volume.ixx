export module mc_moves_gibbs_volume;

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
  std::optional<std::pair<RunningEnergy, RunningEnergy>>
  GibbsVolumeMove(RandomNumber &random, System &systemA, System &systemB);
}
