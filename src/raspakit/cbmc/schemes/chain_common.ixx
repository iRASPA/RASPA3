module;

export module cbmc_chain_common;

import std;

import atom;
import component;
import running_energy;
import cbmc_results;
import cbmc_grow_context;
import cbmc_grow_step;

// What the two chain schemes (configurational-bias and recoil growth) share: the energy evaluation of
// one trial direction, the not-sampled internal-term factor, the per-step overlap guard, the
// accumulation of the chain's log weight and external energies, the assembly of the results (internal
// energies, molecule record), and the error for an existing configuration that overlaps. Each scheme
// keeps only its own per-step Rosenbluth factor.
export namespace CBMC
{
/// The energy of one trial direction of a step, split into the external (external-field, framework,
/// inter-molecular) part that the chain accumulates and the intramolecular van der Waals and Coulomb
/// part that only enters the per-step Boltzmann factor.
struct StepTrialEnergy
{
  RunningEnergy external{};
  RunningEnergy intra{};
  [[nodiscard]] double potentialEnergy() const noexcept
  {
    return external.potentialEnergy() + intra.potentialEnergy();
  }
};

/**
 * \brief Evaluates one trial direction of a step: 'positions' (in 'step.nextBeads' order) against the
 * context's background and against the placed part of the chain.
 *
 * std::nullopt when the trial lies in a blocked pocket or overlaps (a closed direction). The step's
 * next-beads of 'chainAtoms' are used as scratch for the intramolecular evaluation and are restored
 * before returning, so the chain is unchanged for the caller.
 */
[[nodiscard]] std::optional<StepTrialEnergy> evaluateStepTrial(const GrowContext &context, const Component &component,
                                                               const GrowStep &step, std::vector<Atom> &chainAtoms,
                                                               std::span<const Atom> positions,
                                                               std::optional<std::size_t> skipBackgroundMolecule);

/**
 * \brief The log factor -beta u_unsampled of a step's classically not-sampled internal terms
 * (Urey-Bradley, inversion/out-of-plane bends, improper torsions, cross terms), evaluated with the
 * step's beads in place in 'chainAtoms'.
 *
 * Zero for flexible attach steps: their base sampler and torsion-spin selection handle these terms
 * inside the operator engine ('stepHandlesUnsampledInternalTerms'), so applying the factor again
 * would double count. Rigid-body, ring-closure, and seed steps fold them in here.
 */
double unsampledInternalLogFactor(double beta, const GrowStep &step, const std::vector<Atom> &chainAtoms);

/**
 * \brief Running totals of a chain grow or retrace: the log Rosenbluth weight and the external
 * energies summed over the steps.
 */
struct ChainAccumulator
{
  double logRosenbluthWeight{0.0};
  RunningEnergy externalEnergies{};

  /**
   * \brief Adds a grown step. Returns false when the step's Rosenbluth factor is below
   * 'minimumRosenbluthFactor' (the caller abandons the grow).
   *
   * The overlap guard is on the per-step factor, not the running product: the cumulative weight of
   * a long chain decays roughly as f^N (f ~ 0.1 per bead for a chain with intra 1-4 charges), so a
   * chain of a few hundred beads is below any fixed absolute threshold in every attempt and growth
   * would always "dead-end" even though no step overlaps. A genuine overlap still trips the per-step
   * test, since a surviving trial near the 'energyOverlapCriteria' contributes exp(-beta E) << the
   * threshold. The retrace carries no guard, so grow and retrace stay symmetric.
   */
  [[nodiscard]] bool addGrownStep(double stepLogWeight, const RunningEnergy &stepExternalEnergy,
                                  double minimumRosenbluthFactor);

  /// Adds a retraced step (no guard: the old configuration always has a weight).
  void addRetracedStep(double stepLogWeight, const RunningEnergy &stepExternalEnergy);
};

/**
 * \brief The result of a completed grow: the accumulated external energies plus all internal
 * interactions of the grown chain (including the cross terms, recomputed on the final positions),
 * and a valid molecule record (centre of mass, and for a fully rigid component the orientation
 * quaternion recovered from the grown positions).
 */
GrowResult finishGrownChain(const Component &component, std::vector<Atom> chainAtoms,
                            const ChainAccumulator &accumulated);

/// The result of a completed retrace: the accumulated external energies plus the internal energies
/// of the existing chain.
RetraceResult finishRetracedChain(const Component &component, std::span<const Atom> chainAtoms,
                                  const ChainAccumulator &accumulated);

/**
 * \brief Throws the error for an existing configuration that overlaps during a retrace.
 *
 * The old configuration is an accepted state of the simulation: it can not overlap. An overlap means
 * the system is inconsistent (a restart or initial configuration with overlapping molecules, a
 * molecule inside a blocked pocket, or a force field / scaling that changed since the molecule was
 * placed), and no weight is defined for it; weighing it anyway would drive the acceptance with a bogus
 * W_old, so the scheme fails loudly instead. 'scheme' names the scheme in the message.
 */
[[noreturn]] void throwExistingConfigurationOverlaps(std::string_view scheme, const Component &component,
                                                     std::size_t stepIndex, const GrowStep &step);
}  // namespace CBMC
