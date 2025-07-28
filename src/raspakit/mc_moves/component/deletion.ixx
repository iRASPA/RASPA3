module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_deletion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a deletion move in the Monte Carlo simulation.
 *
 * Attempts to delete a molecule from the system, updating energies and statistics accordingly.
 *
 * \param random             Reference to the random number generator.
 * \param system             Reference to the simulation system.
 * \param selectedComponent  Index of the component from which to delete the molecule.
 * \param selectedMolecule   Index of the molecule to delete.
 *
 * \return A pair containing an optional RunningEnergy (energy difference if move accepted) and a double3
 *         representing acceptance probabilities.
 */
std::pair<std::optional<RunningEnergy>, double3> deletionMove(RandomNumber& random, System& system,
                                                              std::size_t selectedComponent,
                                                              std::size_t selectedMolecule);
}  // namespace MC_Moves
