module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#endif

export module mc_moves_random_rotation;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <optional>;
import <span>;
#endif

import randomnumbers;
import running_energy;
import atom;
import molecule;
import component;
import system;

export namespace MC_Moves
{
std::optional<RunningEnergy> randomRotationMove(RandomNumber &random, System &system, size_t selectedComponent,
                                                const std::vector<Component> &components, Molecule &molecule,
                                                std::span<Atom> molecule_atoms);
}
