module;

export module ideal_gas_reference;

import std;

import component;
import system;

export namespace Preparation
{
/// How long the isolated molecule is grown for its ideal-gas reference.
///
/// The grows are made in chunks, each with its own generator (seeded with the number of the chunk) and its own copy
/// of the component: the chunks are independent of each other and of the thread that makes them, so that the
/// reference does not depend on the number of threads. The number of grows is doubled, from the minimum up to the
/// maximum, until the standard error of the mean weight is below the tolerance.
struct IdealGasSettings
{
  std::size_t growsPerChunk{2500};
  std::size_t minimumNumberOfGrows{20000};
  std::size_t maximumNumberOfGrows{1280000};
  /// The standard error of the mean Rosenbluth weight over the mean. It is a systematic shift of beta mu for the
  /// whole run that divides by the weight. A component that was given a tolerance of its own
  /// ('Component::idealGasRosenbluthWeightTolerance') is computed to that one.
  double relativeTolerance{5.0e-3};
};

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
  /// The standard error of the mean weight over the mean weight.
  double relativeError{0.0};
  /// The number of grows the weight is the mean of (the failed ones, of weight zero, among them).
  std::size_t numberOfGrows{0};
  /// Whether 'relativeError' came below the tolerance before the maximum number of grows.
  bool converged{false};
};

/// A molecule with conformations that is exchanged with a reservoir. Its insertions draw (or grow) the conformation
/// from the intramolecular Boltzmann distribution, so that a replica samples the activity beta f / q(beta).
bool isFlexibleOpenComponent(const Component& component);

/// Whether the insertions and deletions of component 'componentId' of 'system' divide beta f by its ideal-gas
/// Rosenbluth weight: those that grow the molecule do, with the weight that was given also for a rigid molecule,
/// those that place a drawn conformation do not. (A move that takes other components along counts for every
/// component: which ones it takes is not looked at.) Where they do not, a weight that was given is in no sampled
/// distribution, and is not in the activity that two replicas are compared with.
bool insertionsDivideByIdealGasRosenbluthWeight(const System& system, std::size_t componentId);

/// Whether a move of 'system' divides by the ideal-gas Rosenbluth weight of component 'componentId': a molecule
/// with conformations that is inserted, deleted, probed or converted with configurational bias. (A rigid molecule
/// and a single bead have the weight one.)
bool needsIdealGasRosenbluthWeight(const System& system, std::size_t componentId);

/// The ideal-gas Rosenbluth weight of component 'componentId' at the temperature of 'system': the mean weight of
/// grows of the isolated molecule, with the growth settings, the force field and the conformation reservoirs of
/// 'system' ('logPartitionFunction' is zero).
IdealGasReference idealGasRosenbluthWeight(const System& system, std::size_t componentId,
                                           const IdealGasSettings& settings = {});

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
/// is wide). The generators are seeded with fixed numbers: the same input gives the same references, on any number
/// of threads.
std::vector<IdealGasReference> idealGasReferences(const System& system, std::size_t componentId,
                                                  std::span<const double> temperatures,
                                                  const IdealGasSettings& settings = {});

/// Gives every component of 'systems' that needs one ('needsIdealGasRosenbluthWeight') and has none its ideal-gas
/// Rosenbluth weight at the temperature of its system, with 'Component::idealGasRosenbluthWeightEstimate' telling
/// that it was computed. A weight that was given is left as it is; one that was computed at another temperature
/// is computed again. Each system is taken on its own: the systems need not be replicas of each other.
void prepareIdealGasRosenbluthWeights(std::span<System> systems, const IdealGasSettings& settings = {});

/// Gives the swappable flexible components of the replicas of a ladder what belongs to the temperature of each:
/// the ideal-gas Rosenbluth weight (the entry of 'Component::idealGasRosenbluthWeightLadder' for the temperature of
/// the replica in 'ladderTemperatures', or the computed one when no weight was given) and, when the replicas are
/// not all at one temperature, 'Component::logInsertionNormalization' (see 'idealGasReferences'; without the weight
/// where the insertions do not divide by it, 'insertionsDivideByIdealGasRosenbluthWeight'). What is then
/// still without a weight is given one by 'prepareIdealGasRosenbluthWeights'.
///
/// The systems are replicas of one system (the same components). 'ladderTemperatures' is the ladder of a driver
/// that replicated one system (empty for systems that were declared one by one, each with its own weight).
/// Returns, per component, the ln of the configuration integral at 'otherTemperatures' on the scale of the
/// normalizations of the replicas (zeros for a component that needs none): the reference of a state point that is
/// reweighted to but not sampled.
std::vector<std::vector<double>> prepareIdealGasReferences(std::span<System> systems,
                                                           std::span<const double> ladderTemperatures = {},
                                                           std::span<const double> otherTemperatures = {},
                                                           const IdealGasSettings& settings = {});
}  // namespace Preparation
