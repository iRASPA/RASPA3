module;

export module cbmc_multi_trial_device;

import std;

import atom;
import molecule;
import randomnumbers;
import component;
import spatial_decomposition_settings;
import cbmc;
import cbmc_multi_trial;
import cbmc_device_environment;

/**
 * \file multi_trial_device.ixx
 * \brief The device implementation of the multiple-try grower for rigid molecules.
 *
 * Proposal scheme (the whole thing is one Multiple-Try Metropolis proposal, reversible for the DEVICE potential):
 *
 *  - M groups of k complete rigid trial molecules, each at a uniform random position with a uniform random
 *    orientation, are generated on the host and evaluated on the device in one batch (external energies in
 *    single precision: Lennard-Jones + real-space Ewald against the framework and the other molecules).
 *  - Within a group one trial is Rosenbluth-selected with the device energies; the group's weight is
 *    W_i = (1/k) sum_j exp(-beta U_ij) over its non-overlapping trials (so a group is a k-trial rigid insertion,
 *    the one-stage analogue of the multiple-first-bead + orientation stages of the CPU grower).
 *  - The move's long-range callback (Ewald Fourier, tail; double, on the host) multiplies every group weight,
 *    W~_i = W_i exp(-beta dU_lr,i), and one group is selected with probability W~_i / sum_j W~_j; the mean weight
 *    (1/M) sum_i W~_i is the quantity of the acceptance rule (cbmc_multi_trial).
 *  - The retrace side evaluates the existing molecule with the SAME kernel, pairs it with k-1 fresh trials for
 *    its group weight and M-1 reference groups, exactly as the grow side.
 *
 * Precision ("option 2"): the selection and the mean weights use the device energies; the ONE selected trial
 * (and the retraced molecule) is re-evaluated on the host in double with the ordinary CBMC external-energy
 * routine, which is what the running energies receive. The acceptance is corrected by the ratio of the double
 * and the device Boltzmann factors at the end points ('logProposalCorrection'), which makes the composite move
 * sample the double-precision distribution exactly (the proposal is reversible for the device potential, and the
 * correction is the Metropolis factor between the two potentials, as in the dual cut-off scheme).
 *
 * Scope: fully rigid components of at most 'DeviceGrowEnvironment::maximumAtomsPerSet' atoms, unit scaling (no
 * fractional molecules), no external field, no cross-link tethers, force fields of plain 12-6 Lennard-Jones and
 * Ewald (or no) electrostatics. Everything else is forwarded to the CPU grower, so a system can mix supported and
 * unsupported components (every component is always grown by one and the same scheme, which is all detailed
 * balance needs). k is 'GrowthSettings::numberOfTrialDirections'.
 */
export namespace CBMC
{
class DeviceMultiTrialGrower final : public MultiTrialGrower
{
 public:
  /// Creates the device environment; throws when no device of the kind is available.
  explicit DeviceMultiTrialGrower(PairDevice device);

  PairDevice device() const;
  std::string deviceName() const;

  /// Whether the device scheme grows 'component' in 'context' (else 'reason' says why the CPU grower is used).
  static bool supports(const GrowContext &context, const Component &component, std::string &reason);

  [[nodiscard]] std::optional<MultiTrialGrowResult> growNew(RandomNumber &random, const GrowContext &context,
                                                            const Component &component,
                                                            const NewMoleculeIdentity &identity,
                                                            std::size_t numberOfTrials,
                                                            const LongRangeCorrection &longRange) const override;

  [[nodiscard]] MultiTrialRetraceResult retraceNew(RandomNumber &random, const GrowContext &context,
                                                   const Component &component, std::span<const Atom> moleculeAtoms,
                                                   const NewMoleculeIdentity &referenceIdentity,
                                                   std::size_t numberOfTrials,
                                                   const LongRangeCorrection &longRange) const override;

  [[nodiscard]] std::optional<MultiTrialGrowResult> regrow(RandomNumber &random, const GrowContext &context,
                                                           const Component &component, const Molecule &molecule,
                                                           std::span<const Atom> moleculeAtoms,
                                                           std::size_t numberOfTrials,
                                                           const LongRangeCorrection &longRange) const override;

  [[nodiscard]] MultiTrialRetraceResult retraceRegrow(RandomNumber &random, const GrowContext &context,
                                                      const Component &component, const Molecule &molecule,
                                                      std::span<const Atom> moleculeAtoms, double storedR,
                                                      std::size_t numberOfTrials,
                                                      const LongRangeCorrection &longRange) const override;

 private:
  // the environment holds the device state of the current evaluation; the grower's interface is const
  mutable std::unique_ptr<DeviceGrowEnvironment> environment;
  // One grower per device kind is shared by every system of the process. The multi-system drivers (parallel
  // tempering, parallel TMMC, ...) run their replicas on concurrent threads, so the device section of a batch
  // (background upload + evaluation, which rewrites the environment's buffers and staging) is serialized here;
  // the host side of a batch (trial placement, selection, double re-evaluation) only touches locals and runs
  // concurrently.
  mutable std::mutex environmentMutex;

  std::optional<MultiTrialGrowResult> growBatch(RandomNumber &random, const GrowContext &context,
                                                const Component &component, std::span<const Atom> templateAtoms,
                                                std::size_t numberOfTrials, const LongRangeCorrection &longRange) const;
  MultiTrialRetraceResult retraceBatch(RandomNumber &random, const GrowContext &context, const Component &component,
                                       std::span<const Atom> moleculeAtoms, std::span<const Atom> templateAtoms,
                                       std::size_t numberOfTrials, const LongRangeCorrection &longRange) const;
};

/// The process-wide device grower of a device kind (created on first use; throws when unavailable).
const MultiTrialGrower &deviceMultiTrialGrower(PairDevice device);
}  // namespace CBMC
