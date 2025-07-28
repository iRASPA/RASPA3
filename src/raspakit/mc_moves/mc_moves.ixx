module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module mc_moves;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import energy_status;
import running_energy;
import mc_moves_translation;
export namespace MC_Moves
{
/**
 * \brief Performs a random Monte Carlo move on the selected system.
 *
 * This function selects and executes a random Monte Carlo move for the specified component in the selected system,
 * based on the configured move probabilities.
 *
 * The function updates the system's state and running energies as necessary after the move.
 *
 * \param random Reference to the random number generator.
 * \param selectedSystem The system on which to perform the move.
 * \param selectedSecondSystem A secondary system, used in moves involving two systems (e.g., Gibbs ensemble moves).
 * \param selectedComponent The index of the component on which the move is to be performed.
 * \param fractionalMoleculeSystem Reference to the system index holding the fractional molecule (used in CFCMC moves).
 */
void performRandomMove(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                       std::size_t selectedComponent, std::size_t& fractionalMoleculeSystem);

/**
 * \brief Performs a random Monte Carlo move during production runs, with statistics tracking.
 *
 * This function selects and executes a random Monte Carlo move for the specified component in the selected system,
 * similar to performRandomMove. Additionally, it records statistics such as move counts, acceptance rates, and CPU
 * time, which are important for analyzing the simulation's performance during production runs.
 *
 * \param random Reference to the random number generator.
 * \param selectedSystem The system on which to perform the move.
 * \param selectedSecondSystem A secondary system, used in moves involving two systems (e.g., Gibbs ensemble moves).
 * \param selectedComponent The index of the component on which the move is to be performed.
 * \param fractionalMoleculeSystem Reference to the system index holding the fractional molecule (used in CFCMC moves).
 * \param currentBlock The current block number, used for statistics aggregation.
 */
void performRandomMoveProduction(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                                 std::size_t selectedComponent, std::size_t& fractionalMoleculeSystem,
                                 std::size_t currentBlock);
};  // namespace MC_Moves
