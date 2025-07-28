module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module mc_moves_random_translation;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import running_energy;
import atom;
import molecule;
import component;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a random translation move for a given molecule.
 *
 * Attempts to translate a molecule in a random direction and evaluates
 * the energy difference resulting from the move. The move is accepted or
 * rejected based on the Metropolis criterion.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param components List of all components in the system.
 * @param molecule The molecule to be moved.
 * @param molecule_atoms Span of atoms belonging to the molecule.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> randomTranslationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                   std::size_t selectedMolecule,
                                                   const std::vector<Component> &components, Molecule &molecule,
                                                   std::span<Atom> molecule_atoms);
}  // namespace MC_Moves
