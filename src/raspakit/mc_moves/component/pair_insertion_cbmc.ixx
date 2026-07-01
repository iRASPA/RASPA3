module;

export module mc_moves_pair_insertion_cbmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Distance-biased ion-pair CBMC insertion (Orkoulas & Panagiotopoulos).
 *
 * Inserts a pair of molecules from linked components. The second ion is placed at a distance r sampled
 * uniformly in [0, R_max] from the first ion's starting bead, with acceptance corrected by 3 r^2 / R_max^2.
 */
std::pair<std::optional<RunningEnergy>, double3> pairInsertionMoveCBMC(RandomNumber& random, System& system,
                                                                       std::size_t selectedComponent);

/**
 * \brief Conventional ion-pair CBMC insertion without distance bias.
 *
 * Inserts a linked pair of molecules. The second ion is placed at a distance r sampled uniformly in the
 * sphere of radius R_max around the first ion's starting bead.
 */
std::pair<std::optional<RunningEnergy>, double3> pairInsertionMove(RandomNumber& random, System& system,
                                                                   std::size_t selectedComponent);
}  // namespace MC_Moves
