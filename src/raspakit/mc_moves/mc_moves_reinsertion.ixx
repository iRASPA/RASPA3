module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#endif

export module mc_moves_reinsertion;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <optional>;
import <span>;
#endif

import randomnumbers;
import running_energy;
import atom;
import molecule;
import system;

export namespace MC_Moves
{
std::optional<RunningEnergy> reinsertionMove(RandomNumber &random, System &system, size_t selectedComponent,
                                             size_t selectedMolecule, Molecule &molecule,
                                             std::span<Atom> molecule_atoms);
}
