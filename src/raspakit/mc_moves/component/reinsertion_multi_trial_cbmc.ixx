module;

export module mc_moves_reinsertion_multi_trial_cbmc;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Multiple-try CBMC reinsertion: the molecule is regrown M times at random positions, one trial is
 *        selected with probability proportional to its (long-range corrected) Rosenbluth weight, and the
 *        old configuration is retraced together with M-1 reference regrowths.
 *
 * \return The energy difference when accepted; std::nullopt otherwise.
 */
std::optional<RunningEnergy> reinsertionMoveMultiTrialCBMC(RandomNumber& random, System& system,
                                                           std::size_t selectedComponent, std::size_t selectedMolecule);
}  // namespace MC_Moves
