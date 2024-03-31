export module mc_moves_reinsertion;

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
  reinsertionMove(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule, 
                  std::span<Atom> molecule);
}

