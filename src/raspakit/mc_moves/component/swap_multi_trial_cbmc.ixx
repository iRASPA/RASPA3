module;

export module mc_moves_swap_multi_trial_cbmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Multiple-try CBMC insertion: M trial molecules are grown, one is selected with probability
 *        proportional to its (long-range corrected) Rosenbluth weight, and the move is accepted with the
 *        mean weight of the trial set.
 *
 * \return The energy difference when accepted, and the (deletion, stay, insertion) acceptance triple
 *         for the transition matrix.
 */
std::pair<std::optional<RunningEnergy>, double3> insertionMoveMultiTrialCBMC(RandomNumber& random, System& system,
                                                                             std::size_t selectedComponent);

/**
 * \brief Multiple-try CBMC deletion: the selected molecule is retraced and M-1 reference molecules are grown
 *        into the system without it; the move is accepted with the mean weight of that reference set.
 */
std::pair<std::optional<RunningEnergy>, double3> deletionMoveMultiTrialCBMC(RandomNumber& random, System& system,
                                                                            std::size_t selectedComponent,
                                                                            std::size_t selectedMolecule);
}  // namespace MC_Moves
