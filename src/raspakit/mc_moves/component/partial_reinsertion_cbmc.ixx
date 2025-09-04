module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#endif

export module mc_moves_partial_reinsertion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import running_energy;
import atom;
import molecule;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a CBMC reinsertion move for a selected molecule.
 *
 * Attempts to reinsert a molecule using Configurational Bias Monte Carlo (CBMC) method.
 * Updates the system if the move is accepted and returns the energy difference.
 *
 * \param random Random number generator instance.
 * \param system The simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule within the component.
 * \param molecule Reference to the molecule to be reinserted.
 * \param molecule_atoms Span of atoms belonging to the molecule.
 *
 * \return An optional RunningEnergy containing the energy difference if the move is accepted; std::nullopt otherwise.
 */
std::optional<RunningEnergy> partialReinsertionMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                    std::size_t selectedMolecule, Molecule &molecule,
                                                    std::span<Atom> molecule_atoms);
}  // namespace MC_Moves
