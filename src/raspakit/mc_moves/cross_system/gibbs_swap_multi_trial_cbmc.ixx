module;

export module mc_moves_gibbs_swap_multi_trial_cbmc;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Multiple-try CBMC Gibbs swap: a molecule of the component is transferred from system B to system A.
 *        M trial molecules are grown in A and one is selected by its (long-range corrected) Rosenbluth weight;
 *        the molecule leaving B is retraced together with M-1 reference molecules grown into B without it.
 *
 * \return The energy differences of systems A and B when accepted; std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMoveMultiTrialCBMC(RandomNumber& random,
                                                                                   System& systemA, System& systemB,
                                                                                   std::size_t selectedComponent);
}  // namespace MC_Moves
