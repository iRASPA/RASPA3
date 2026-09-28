module;

export module mc_moves_parallel_tempering_swap;

import std;

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a parallel tempering swap between two systems.
 *
 * Attempts to swap the configurations of two systems in a parallel tempering Monte Carlo simulation.
 * The swap is accepted or rejected based on the Metropolis criterion, taking into account differences
 * in temperature, pressure, potential energy, and, for matching CFCMC replicas, the replica-local
 * lambda (and TMMC) biases at the incoming lambda coordinates. Rigid, flexible and semi-flexible
 * components all travel with the configuration (positions, rigid-body state and cross-links), the
 * topology being the same component definition in every replica.
 *
 * Reference: "Hyper-parallel tempering Monte Carlo: Application to the Lennard-Jones fluid and the
 * restricted primitive model", G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999.
 *
 * \param random   Random number generator used for acceptance probability.
 * \param systemA  First system involved in the swap.
 * \param systemB  Second system involved in the swap.
 * \return         An optional pair of RunningEnergy objects if the swap is accepted; std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> ParallelTemperingSwap(RandomNumber &random, System &systemA,
                                                                             System &systemB);

/// Log acceptance ratio for X_A ↔ X_B. Empty when the pair is incompatible
/// (no RNG should be consumed). Exposed so tests can check the ratio.
std::optional<double> ParallelTemperingLogAcceptance(const System &systemA, const System &systemB);

/**
 * \brief Replica-exchange swap between two molecular-dynamics replicas.
 *
 * Performs the configuration swap of ParallelTemperingSwap (same acceptance rule, based on the
 * potential energies) and then brings the molecular-dynamics state of both replicas in line with
 * their (unchanged) temperatures: the momenta that travelled with the configuration are rescaled
 * by sqrt(T_new / T_old) (Sugita & Okamoto, Chem. Phys. Lett. 314, 141-151, 1999), so the kinetic
 * energy remains canonical at the replica temperature and the momentum part of the acceptance
 * rule cancels; the gradients, kinetic energies and thermostat bookkeeping are recomputed and the
 * conserved-energy reference of both replicas is reset (the extended-system energy jumps at a
 * swap by construction). The thermostat chain is a property of the heat bath and stays with the
 * replica.
 *
 * \param random   Random number generator used for the acceptance test.
 * \param systemA  First replica.
 * \param systemB  Second replica.
 * \return         The running energies of both replicas after the swap if accepted; std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> ParallelTemperingSwapMolecularDynamics(RandomNumber &random,
                                                                                              System &systemA,
                                                                                              System &systemB);
}  // namespace MC_Moves
