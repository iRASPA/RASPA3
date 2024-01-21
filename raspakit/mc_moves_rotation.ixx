export module mc_moves_rotation;

import <cstddef>;
import <optional>;
import <span>;

import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
  std::optional<RunningEnergy>
  rotationMove(RandomNumber &random, System &system, size_t selectedComponent, std::span<Atom> molecule);
}
