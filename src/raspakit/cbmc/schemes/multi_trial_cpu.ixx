module;

export module cbmc_multi_trial_cpu;

import std;

import atom;
import molecule;
import randomnumbers;
import component;
import cbmc;
export import cbmc_multi_trial;

// The reference implementation of the multiple-try grower: M independent calls to the single-molecule
// entry points of 'cbmc', each trial corrected by the move's long-range callback before the selection.
// It is the definition of the sampled distribution that any faster implementation (the device grower)
// is checked against, and with M = 1 it reproduces the ordinary CBMC moves draw for draw.
export namespace CBMC
{
class CPUMultiTrialGrower final : public MultiTrialGrower
{
 public:
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
};

/// The process-wide CPU grower (stateless).
[[nodiscard]] const MultiTrialGrower &cpuMultiTrialGrower();
}  // namespace CBMC
