module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module mc_moves_random_rotation;

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
 * \brief Performs a random rotation move on a molecule within the system.
 *
 * Attempts to rotate a selected molecule around one of the principal axes by a random angle.
 * Calculates the energy differences resulting from the rotation, including contributions from
 * external fields, framework interactions, intermolecular interactions, and Ewald summation.
 * The move is accepted or rejected based on the Metropolis acceptance criterion.
 *
 * \param random A reference to the random number generator.
 * \param system The simulation system containing all molecules and settings.
 * \param selectedComponent The index of the selected component to rotate.
 * \param components A vector of all components in the system.
 * \param molecule The molecule to be rotated.
 * \param molecule_atoms A span of atoms belonging to the molecule.
 * \return An optional RunningEnergy object containing the energy difference if the move is accepted;
 *         returns std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> randomRotationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                std::size_t selectedMolecule, const std::vector<Component> &components,
                                                Molecule &molecule, std::span<Atom> molecule_atoms);
}  // namespace MC_Moves
