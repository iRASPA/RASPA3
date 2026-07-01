module;

export module mc_moves_pair_deletion_cbmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Distance-biased ion-pair CBMC deletion (Orkoulas & Panagiotopoulos).
 *
 * Deletes a linked pair of molecules. The distance bias factor R_max^2 / (3 r^2) is included in the acceptance.
 */
std::pair<std::optional<RunningEnergy>, double3> pairDeletionMoveCBMC(RandomNumber& random, System& system,
                                                                       std::size_t selectedComponent,
                                                                       std::size_t selectedMolecule);

/**
 * \brief Conventional ion-pair CBMC deletion without distance bias.
 *
 * Deletes a linked pair of molecules without the Orkoulas distance bias correction.
 */
std::pair<std::optional<RunningEnergy>, double3> pairDeletionMove(RandomNumber& random, System& system,
                                                                  std::size_t selectedComponent,
                                                                  std::size_t selectedMolecule);
}  // namespace MC_Moves
