export module mc_moves_insertion_cbmc;

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
  std::pair<std::optional<RunningEnergy>, double3>
  insertionMoveCBMC(RandomNumber &random, System& system, size_t selectedComponent);
}

