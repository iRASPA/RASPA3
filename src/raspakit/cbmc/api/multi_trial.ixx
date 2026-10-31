module;

export module cbmc_multi_trial;

import std;

import atom;
import molecule;
import randomnumbers;
import component;
import cbmc;

// Multiple-try configurational-bias growth: the whole-molecule generalization of the multiple-first-bead
// scheme (Esselink, Loyens, Smit, Phys. Rev. E 51, 1560 (1995); the Multiple-Try Metropolis of Liu, Liang,
// Wong, JASA 95, 121 (2000)). A move grows M complete trial molecules, selects one with probability
// W~_i / sum_j W~_j, and accepts with the MEAN weight
//
//   W-bar_new = (1/M) sum_{i=1}^{M} W~_i
//
// in place of the single Rosenbluth weight; the reverse move pairs the retraced weight of the existing
// molecule with M-1 freshly grown reference molecules,
//
//   W-bar_old = (1/M) ( W~_old + sum_{j=1}^{M-1} W~'_j ).
//
// W~_i is the Rosenbluth weight of trial i times the Boltzmann factor of its long-range contributions
// (Ewald Fourier, tail corrections): the selection probability must see every term the acceptance sees,
// otherwise the chain samples the wrong distribution. The move owns the long-range evaluation (it knows
// the system and its structure factors) and hands it to the grower as 'LongRangeCorrection'; the grower
// is agnostic of how it is obtained. With M = 1 the scheme reduces exactly to the ordinary CBMC move.
//
// The grower is an interface so that the moves are device-agnostic: 'cbmc_multi_trial_cpu' is the
// reference implementation (M calls to the single-molecule entry points) that the device implementation
// is held to.
export namespace CBMC
{
/// The long-range energy of inserting 'atoms' into the move's environment (Ewald Fourier + tail), in the
/// internal energy unit; the grower multiplies each trial's weight by exp(-beta * value). Invoked once per
/// constructed trial, and once for the existing molecule on the retrace side.
using LongRangeCorrection = std::function<double(std::span<const Atom>)>;

/// The outcome of growing M trial molecules and selecting one.
struct MultiTrialGrowResult
{
  /// The selected trial: its atoms, external energies (at the full cut-offs), and log Rosenbluth weight
  /// (without the long-range factor).
  GrowResult selected;
  /// The long-range energy of the selected trial as returned by the correction callback.
  double selectedLongRangeEnergy{0.0};
  /// log W-bar = log( (1/M) sum_i W~_i ), the quantity the acceptance rule uses.
  double logMeanWeight{0.0};
  /// Number of trials that could be constructed (non-overlapping, above the minimum factor).
  std::size_t numberOfConstructed{0};
  /// Proposal correction of a grower that proposes with a different (cheaper) potential than the one the
  /// system accounts for: -beta (U_target(selected) - U_proposal(selected)), to be ADDED to the log acceptance.
  /// The selection and 'logMeanWeight' are then those of the proposal potential, 'selected.energies' those of
  /// the target potential, and the correction restores detailed balance for the target (the proposal kernel is
  /// reversible for its own potential; the factor is the ratio of the two Boltzmann factors at the end points,
  /// as in the dual cut-off scheme). Zero for the CPU grower, which proposes with the target potential.
  double logProposalCorrection{0.0};
};

/// The outcome of retracing an existing molecule together with its M-1 reference trials.
struct MultiTrialRetraceResult
{
  /// The retrace of the existing molecule: energies and log Rosenbluth weight (without the long-range factor).
  RetraceResult old;
  /// The long-range energy of the existing molecule as returned by the correction callback.
  double oldLongRangeEnergy{0.0};
  /// log W-bar = log( (1/M) ( W~_old + sum_j W~'_j ) ).
  double logMeanWeight{0.0};
  /// Number of reference trials that could be constructed.
  std::size_t numberOfConstructed{0};
  /// The retrace-side counterpart of 'MultiTrialGrowResult::logProposalCorrection':
  /// +beta (U_target(old) - U_proposal(old)), to be ADDED to the log acceptance. Zero for the CPU grower. When the
  /// proposal potential gives the existing configuration zero weight (the grow side could never have produced it),
  /// the correction is 'rejectedLogProposalCorrection': the move's acceptance probability underflows to zero.
  double logProposalCorrection{0.0};
};

/// The log proposal correction of a configuration the proposal potential can not produce: a finite stand-in for
/// minus infinity (the code is compiled with fast-math, so an infinity could not be relied on), chosen so that it
/// dominates every other finite term of a log acceptance and its exponential underflows to exactly zero.
inline constexpr double rejectedLogProposalCorrection = -1.0e300;

/**
 * \brief Grows M complete trial molecules and selects one (or retraces one against M-1 references).
 *
 * Two pairs of entry points: a NEW molecule (insertion / deletion; the trials use the multiple-first-bead
 * scheme and carry 'identity') and the REGROWTH of an existing molecule (reinsertion; the trials use the
 * Esselink reinsertion first-bead scheme, and the retrace of the old configuration takes the retained
 * partial first-bead weight of the selected trial). Every implementation must draw the same random
 * experiment per trial as the single-molecule entry points of 'cbmc', so that M = 1 reproduces them.
 */
class MultiTrialGrower
{
 public:
  virtual ~MultiTrialGrower() = default;

  /// Grows M new molecules of 'component' and selects one with probability proportional to W~_i.
  /// Returns std::nullopt when no trial could be constructed (an ordinary rejection).
  [[nodiscard]] virtual std::optional<MultiTrialGrowResult> growNew(RandomNumber &random, const GrowContext &context,
                                                                    const Component &component,
                                                                    const NewMoleculeIdentity &identity,
                                                                    std::size_t numberOfTrials,
                                                                    const LongRangeCorrection &longRange) const = 0;

  /// Retraces the existing molecule 'moleculeAtoms' (multiple-first-bead weighting) and grows M-1 reference
  /// molecules with 'referenceIdentity' against the context's background. The context must exclude the
  /// existing molecule from the reference grows (the same-id rule does when 'referenceIdentity.moleculeId'
  /// is the molecule's own id).
  [[nodiscard]] virtual MultiTrialRetraceResult retraceNew(RandomNumber &random, const GrowContext &context,
                                                           const Component &component,
                                                           std::span<const Atom> moleculeAtoms,
                                                           const NewMoleculeIdentity &referenceIdentity,
                                                           std::size_t numberOfTrials,
                                                           const LongRangeCorrection &longRange) const = 0;

  /// Regrows the existing molecule M times (reinsertion first-bead scheme) and selects one trial.
  [[nodiscard]] virtual std::optional<MultiTrialGrowResult> regrow(RandomNumber &random, const GrowContext &context,
                                                                   const Component &component,
                                                                   const Molecule &molecule,
                                                                   std::span<const Atom> moleculeAtoms,
                                                                   std::size_t numberOfTrials,
                                                                   const LongRangeCorrection &longRange) const = 0;

  /// Retraces the existing molecule with the reinsertion scheme ('storedR' is the retained partial weight
  /// of the selected trial of the matching 'regrow') and regrows M-1 reference trials.
  [[nodiscard]] virtual MultiTrialRetraceResult retraceRegrow(RandomNumber &random, const GrowContext &context,
                                                              const Component &component,
                                                              const Molecule &molecule,
                                                              std::span<const Atom> moleculeAtoms, double storedR,
                                                              std::size_t numberOfTrials,
                                                              const LongRangeCorrection &longRange) const = 0;
};

/// log( sum_i exp(x_i) ) over the finite entries; -infinity when none is finite.
/// Every entry must be a finite log-weight (the callers only pass the weights of the trials that were actually
/// constructed; the code is compiled with fast-math, so infinities can not serve as "no trial" markers).
/// Precondition: 'values' is not empty.
[[nodiscard]] inline double logSumExp(std::span<const double> values) noexcept
{
  const double maximum = *std::ranges::max_element(values);
  double sum = 0.0;
  for (double value : values) sum += std::exp(value - maximum);
  return maximum + std::log(sum);
}

/// Draws an index with probability exp(logWeights[i]) / sum_j exp(logWeights[j]). Precondition: at least one
/// entry, all finite.
[[nodiscard]] inline std::size_t selectByLogWeight(RandomNumber &random, std::span<const double> logWeights)
{
  const double logTotal = logSumExp(logWeights);
  const double target = random.uniform();
  double cumulative = 0.0;
  for (std::size_t i = 0; i != logWeights.size(); ++i)
  {
    cumulative += std::exp(logWeights[i] - logTotal);
    if (target < cumulative) return i;
  }
  return logWeights.size() - 1;  // round-off: the cumulative sum fell short of one
}
}  // namespace CBMC
