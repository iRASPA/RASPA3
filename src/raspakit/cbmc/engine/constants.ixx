module;

export module cbmc_constants;

import std;

import units;

/**
 * \brief Numerical constants of the CBMC operator engine and its setup-time estimates, in one place.
 *
 * None of them affects the sampled distribution (the samplers are exact or Metropolis-correct for any
 * value); they set fall-backs for under-specified topologies, the resolution of one-off numerical
 * estimates, and the budgets of rejection loops.
 */
export namespace CBMC::Constants
{
/// Bond length (Angstrom) used when a step's bond carries no potential (e.g. a connectivity entry
/// without a matching 'Bonds' term): a fixed C-C single-bond length. Also the value the normalization
/// integrates over in that case.
constexpr double defaultBondLength = 1.54;

/// Junction bend angle used by the rigid tilt when no previous-current-inner bend potential exists:
/// a trigonal 120 degrees.
constexpr double defaultRigidJunctionBendAngle = 120.0 * Units::DegreesToRadians;

/// Length (Angstrom) below which a previous-current vector is treated as degenerate (the two beads
/// coincide) and replaced by the z axis; compared against the vector's length, not its square.
constexpr double degenerateAxisLength = 1e-8;

/// Rejection budget of the exact flexible base sampler; exceeding it signals a pathologically stiff
/// coupling or an unsatisfiable chiral constraint and throws.
constexpr std::size_t baseSamplerMaximumAttempts = 1'000'000;

/// Random re-orientation attempts of the ring seed geometry to recover the declared parity of a
/// chiral centre before giving up on that seed.
constexpr std::size_t ringParityMaximumReorientations = 1000;

/// Grid over [0, pi] for the numerical minimum of a non-harmonic bend potential, and the slack
/// (Kelvin) subtracted because a grid minimum can only overestimate the true minimum.
constexpr std::size_t bendMinimumGridPoints = 8192;
constexpr double bendMinimumSlack = 0.01;

/// Fixed seed of the one-off Monte-Carlo estimate of the base-coupling constants: the estimate must
/// be one frozen number per step signature so grow and retrace share it.
constexpr std::size_t baseCouplingEstimateSeed = 1806;
/// Samples of the minimum search of the coupling energy, and the slack (Kelvin) below that minimum
/// for the rejection reference (see estimateBaseCouplingConstants; exactness does not depend on it).
constexpr std::size_t baseCouplingMinimumSearchSamples = 1uz << 20;
constexpr double baseCouplingReferenceSlack = 1.0;
/// Batching of the clamped-acceptance mean: batches of 'batchSize' samples until the standard error
/// drops below 'relativeTolerance' of the mean, or 'maximumBatches' is reached.
constexpr std::size_t baseCouplingBatchSize = 1uz << 20;
constexpr std::size_t baseCouplingMaximumBatches = 32;
constexpr double baseCouplingRelativeTolerance = 0.005;

/// Roll angles about the junction bond tried when seating a hinged rigid body (a Rosenbluth selection
/// over this grid, randomly offset, seeds the tilt Monte-Carlo).
constexpr std::size_t rigidTiltRollGridPoints = 72;

/// The recoil-growth openness reference ('System::buildRecoilReferenceConformations'): per growth step
/// the maximum intramolecular strain over this many equilibrated ideal-gas conformations of the
/// component, grown once at setup from a generator with this fixed seed. The reference is a constant
/// of the run (grow and retrace divide by the same openness probabilities), so the build must be
/// deterministic and independent of the simulation's random stream.
constexpr std::size_t recoilReferenceConformations = 50;
constexpr std::size_t recoilReferenceSeed = 1867;

/// The closure-guide tables of fixed-endpoint regrowth ('cbmc_closure_guide'): the grid spacing
/// (Angstrom) of the bead-target distance, the fixed seed of the one-off tabulation (one frozen table
/// per path signature, shared by grow and retrace), the number of bond-length pairs the two-bond
/// table integrates over, the number of ideal sub-chains sampled for a longer path, the floor of the
/// guide relative to its maximum (as a log; e^-14 ~ 1e-6), the rejection budget of the torsion
/// sampler of the ideal sub-chain, and the grid of the numerical torsion minimum that sampler needs.
/// None affects the sampled distribution: the guide is divided out of the Rosenbluth weight again.
constexpr double closureGuideSpacing = 0.01;
constexpr std::size_t closureGuideSeed = 1901;
constexpr std::size_t closureGuideBondPairs = 4096;
constexpr std::size_t closureGuideChainSamples = 400'000;
constexpr double closureGuideLogFloor = 14.0;
constexpr std::size_t closureGuideTorsionAttempts = 10'000;
constexpr std::size_t torsionMinimumGridPoints = 3600;
/// Margin (Angstrom) added to the sampled maximum reach of a path when sizing its table, and the
/// number of draws per bond used to estimate that reach.
constexpr double closureGuideReachMargin = 0.2;
constexpr std::size_t closureGuideReachSamples = 4096;
}  // namespace CBMC::Constants
