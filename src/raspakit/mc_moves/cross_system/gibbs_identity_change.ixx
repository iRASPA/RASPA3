module;

export module mc_moves_gibbs_identity_change;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a Gibbs ensemble CBMC identity-change move between two systems.
 *
 * In box I, a molecule of component A is replaced by a molecule of component B.
 * In box II, a molecule of component B is replaced by a molecule of component A.
 *
 * \param random Random number generator.
 * \param systemI The first Gibbs box.
 * \param systemII The second Gibbs box.
 * \param componentA The index of the component whose identity is changed.
 * \return A pair of RunningEnergy differences for systems I and II if the move is accepted;
 *         std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsIdentityChangeMove_CBMC(RandomNumber& random,
                                                                                    System& systemI, System& systemII,
                                                                                    std::size_t componentA);
}  // namespace MC_Moves
