module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_deletion_cbmc;

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
 * \brief Performs a CBMC deletion move for a specified molecule.
 *
 * Attempts to delete a molecule of the selected component using the Configurational Bias Monte Carlo (CBMC) algorithm.
 * It calculates the acceptance probability based on energy differences and updates the system accordingly.
 *
 * \param random            The random number generator to use for stochastic processes.
 * \param system            The simulation system containing all components and molecules.
 * \param selectedComponent The index of the component from which the molecule will be deleted.
 * \param selectedMolecule  The index of the molecule to delete within the selected component.
 *
 * \return A pair consisting of an optional RunningEnergy (energy difference if move is accepted) and a double3
 * containing acceptance probabilities.
 */
std::pair<std::optional<RunningEnergy>, double3> deletionMoveCBMC(RandomNumber& random, System& system,
                                                                  std::size_t selectedComponent,
                                                                  std::size_t selectedMolecule);
}  // namespace MC_Moves
