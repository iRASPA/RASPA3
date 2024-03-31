export module mc_moves_identity_change;

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
  identityChangeMove([[maybe_unused]] RandomNumber &random, [[maybe_unused]] System& system,
                     [[maybe_unused]] size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule,
                     [[maybe_unused]] std::span<Atom> atoms);
}
