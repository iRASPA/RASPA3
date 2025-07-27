module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module mc_moves_translation;

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
 * \brief Attempts to translate a molecule within the system.
 *
 * This function performs a translation move on a selected molecule by applying
 * a random displacement along a randomly chosen axis. It computes the energy
 * difference resulting from the move and decides whether to accept or reject
 * the move based on the Metropolis criterion.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule within the component.
 * \param components Vector of components in the system.
 * \param molecule Reference to the molecule to be moved.
 * \param molecule_atoms Span of atoms belonging to the molecule.
 *
 * \return An optional RunningEnergy object containing the energy difference if
 *         the move is accepted; returns std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> translationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                             std::size_t selectedMolecule, const std::vector<Component> &components,
                                             Molecule &molecule, std::span<Atom> molecule_atoms);
}  // namespace MC_Moves
