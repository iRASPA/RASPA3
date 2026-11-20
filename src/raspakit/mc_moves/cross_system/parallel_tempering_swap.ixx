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

/// Why replica-exchange swaps can not be done for replicas of 'system' (nothing when they can). The parallel
/// tempering drivers refuse such a system up front; the swap itself refuses it as well. The cases:
///  - a framework with mobile atoms: the swap exchanges the adsorbates only, and a replica would keep its own
///    framework configuration with the partner's energies;
///  - a temperature-dependent pair potential (Feynman-Hibbs): the replicas then have different Hamiltonians,
///    and the acceptance rule needs the energy of each configuration under the partner's potential;
///  - solute tempering with polarization: the solute Hamiltonian change has no polarization term;
///  - solute tempering with a flexible molecule in an open ensemble: the intramolecular configuration integral
///    of the scaled Hamiltonians is not known. (Across temperatures it is computed, see
///    'prepareIdealGasReferences'; the swap refuses a pair of temperatures that was not prepared.)
std::optional<std::string> parallelTemperingUnsupportedReason(const System& system, bool soluteTempering);

/// The ideal-gas reference of a flexible molecule at one temperature.
struct IdealGasReference
{
  double temperature{0.0};
  /// ln of the intramolecular configuration integral q(beta) = int dx exp(-beta U_intra(x)), relative to the lowest
  /// of the temperatures it was computed with.
  double logPartitionFunction{0.0};
  /// ln of the ideal-gas Rosenbluth weight <W>(beta) of the configurational-bias growth of the isolated molecule.
  double logRosenbluthWeight{0.0};
  /// The disagreement of the two directions of the free-energy perturbation to the temperature below, summed
  /// up to this temperature: the scale of the error of 'logPartitionFunction'.
  double logPartitionFunctionError{0.0};
};

/// The ideal-gas references of component 'componentId' of 'system' at 'temperatures' (returned in their order).
///
/// An open ensemble with a flexible molecule samples pi ~ (beta f / q(beta))^N exp(-beta U), with q the
/// intramolecular configuration integral of the molecule: the growth draws the bonded coordinates from normalized
/// Boltzmann distributions, whose normalization Z_b(beta) is in the sampled distribution next to the ideal-gas
/// Rosenbluth weight <W>(beta) = q(beta) / Z_b(beta) the acceptance divides by. Two temperatures are compared (an
/// exchange of configurations with different molecule counts, a reweighting) through q(beta') / q(beta), which the
/// weights alone do not give: a diatomic with one harmonic bond has <W> = 1 at every temperature and q ~ sqrt(T).
///
/// Both follow from grows of the isolated molecule, for which generation density times weight is
/// exp(-beta U_intra) / Z_b: <W> is the mean weight, and q(beta') / q(beta) = <W exp(-(beta' - beta) U_intra)> / <W>
/// (taken in both directions between neighbouring temperatures, with temperatures put in between where the step
/// is wide). The generator is seeded with fixed numbers: the same input gives the same references.
std::vector<IdealGasReference> idealGasReferences(const System& system, std::size_t componentId,
                                                  std::span<const double> temperatures);

/// Gives the swappable flexible components of the replicas of a ladder what belongs to the temperature of each:
/// the ideal-gas Rosenbluth weight (the entry of 'Component::idealGasRosenbluthWeightLadder' for the temperature of
/// the replica in 'ladderTemperatures', or the computed one when no weight was given) and, when the replicas are
/// not all at one temperature, 'Component::logInsertionNormalization' (see 'idealGasReferences').
///
/// 'ladderTemperatures' is the ladder of a driver that replicated one system (empty for systems that were declared
/// one by one, each with its own weight). Returns, per component, the ln of the configuration integral at
/// 'otherTemperatures' on the scale of the normalizations of the replicas (zeros for a component that needs none):
/// the reference of a state point that is reweighted to but not sampled.
std::vector<std::vector<double>> prepareIdealGasReferences(std::span<System> systems,
                                                           std::span<const double> ladderTemperatures = {},
                                                           std::span<const double> otherTemperatures = {});

/// Log acceptance ratio for X_A ↔ X_B. Empty when the pair is incompatible
/// (no RNG should be consumed). Exposed so tests can check the ratio.
std::optional<double> ParallelTemperingLogAcceptance(const System &systemA, const System &systemB);

/// Solute tempering: E_target(X) - E_holder(X) for the configuration X held by 'holder' (the energy change
/// of switching the solute Hamiltonian at fixed positions). Zero when the replicas share the same solute
/// scaling; empty when an energy routine reports an overlap. Exposed so tests can check it.
std::optional<double> ParallelTemperingSoluteHamiltonianChange(const System &holder, const System &target);

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
