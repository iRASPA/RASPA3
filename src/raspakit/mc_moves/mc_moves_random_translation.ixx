export module mc_moves_random_translation;

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
  randomTranslationMove(RandomNumber &random, System &system, size_t selectedComponent, std::span<Atom> molecule);
}
