module;

module cbmc_multi_trial_cpu;

import std;

import atom;
import molecule;
import randomnumbers;
import component;
import cbmc;
import cbmc_multi_trial;

namespace
{
// Grows the M trials with 'growOne', folds the long-range factor into each weight, and selects one.
// Shared by the new-molecule and the regrowth entry points; they differ only in how a trial is grown.
std::optional<CBMC::MultiTrialGrowResult> growAndSelect(
    RandomNumber &random, double beta, std::size_t numberOfTrials, const CBMC::LongRangeCorrection &longRange,
    const std::function<std::optional<CBMC::GrowResult>(RandomNumber &)> &growOne)
{
  const std::size_t M = std::max<std::size_t>(1, numberOfTrials);

  // Only the trials that were actually constructed are kept (an overlapping grow returns no result and carries
  // zero weight); the mean weight still divides by M.
  std::vector<CBMC::GrowResult> trials;
  std::vector<double> longRangeEnergies;
  std::vector<double> logWeights;
  trials.reserve(M);
  longRangeEnergies.reserve(M);
  logWeights.reserve(M);

  for (std::size_t i = 0; i != M; ++i)
  {
    std::optional<CBMC::GrowResult> trial = growOne(random);
    if (!trial.has_value()) continue;
    const double longRangeEnergy = longRange(trial->atoms);
    longRangeEnergies.push_back(longRangeEnergy);
    logWeights.push_back(trial->logRosenbluthWeight - beta * longRangeEnergy);
    trials.push_back(std::move(trial.value()));
  }
  if (trials.empty()) return std::nullopt;

  const std::size_t selected = CBMC::selectByLogWeight(random, logWeights);
  return CBMC::MultiTrialGrowResult{.selected = std::move(trials[selected]),
                                    .selectedLongRangeEnergy = longRangeEnergies[selected],
                                    .logMeanWeight = CBMC::logSumExp(logWeights) - std::log(static_cast<double>(M)),
                                    .numberOfConstructed = trials.size()};
}

// Retraces the existing molecule with 'retraceOld' and grows M-1 references with 'growOne'.
CBMC::MultiTrialRetraceResult retraceWithReference(
    RandomNumber &random, double beta, std::span<const Atom> moleculeAtoms, std::size_t numberOfTrials,
    const CBMC::LongRangeCorrection &longRange, const std::function<CBMC::RetraceResult(RandomNumber &)> &retraceOld,
    const std::function<std::optional<CBMC::GrowResult>(RandomNumber &)> &growOne)
{
  const std::size_t M = std::max<std::size_t>(1, numberOfTrials);
  std::vector<double> logWeights;
  logWeights.reserve(M);

  CBMC::RetraceResult old = retraceOld(random);
  const double oldLongRange = longRange(moleculeAtoms);
  logWeights.push_back(old.logRosenbluthWeight - beta * oldLongRange);

  std::size_t constructed = 0;
  for (std::size_t j = 1; j != M; ++j)
  {
    std::optional<CBMC::GrowResult> reference = growOne(random);
    if (!reference.has_value()) continue;
    ++constructed;
    logWeights.push_back(reference->logRosenbluthWeight - beta * longRange(reference->atoms));
  }

  return CBMC::MultiTrialRetraceResult{.old = old,
                                       .oldLongRangeEnergy = oldLongRange,
                                       .logMeanWeight = CBMC::logSumExp(logWeights) - std::log(static_cast<double>(M)),
                                       .numberOfConstructed = constructed};
}
}  // namespace

std::optional<CBMC::MultiTrialGrowResult> CBMC::CPUMultiTrialGrower::growNew(
    RandomNumber &random, const GrowContext &context, const Component &component, const NewMoleculeIdentity &identity,
    std::size_t numberOfTrials, const LongRangeCorrection &longRange) const
{
  return growAndSelect(random, context.beta, numberOfTrials, longRange,
                       [&](RandomNumber &rng) { return CBMC::growNewMolecule(rng, context, component, identity); });
}

CBMC::MultiTrialRetraceResult CBMC::CPUMultiTrialGrower::retraceNew(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<const Atom> moleculeAtoms,
    const NewMoleculeIdentity &referenceIdentity, std::size_t numberOfTrials,
    const LongRangeCorrection &longRange) const
{
  return retraceWithReference(
      random, context.beta, moleculeAtoms, numberOfTrials, longRange,
      [&](RandomNumber &rng) { return CBMC::retraceMolecule(rng, context, component, moleculeAtoms); },
      [&](RandomNumber &rng) { return CBMC::growNewMolecule(rng, context, component, referenceIdentity); });
}

std::optional<CBMC::MultiTrialGrowResult> CBMC::CPUMultiTrialGrower::regrow(
    RandomNumber &random, const GrowContext &context, const Component &component, const Molecule &molecule,
    std::span<const Atom> moleculeAtoms, std::size_t numberOfTrials, const LongRangeCorrection &longRange) const
{
  return growAndSelect(random, context.beta, numberOfTrials, longRange,
                       [&](RandomNumber &rng)
                       {
                         return CBMC::regrowMolecule(rng, context, component, molecule, moleculeAtoms,
                                                     {.firstBead = FirstBeadScheme::Reinsertion});
                       });
}

CBMC::MultiTrialRetraceResult CBMC::CPUMultiTrialGrower::retraceRegrow(
    RandomNumber &random, const GrowContext &context, const Component &component, const Molecule &molecule,
    std::span<const Atom> moleculeAtoms, double storedR, std::size_t numberOfTrials,
    const LongRangeCorrection &longRange) const
{
  return retraceWithReference(
      random, context.beta, moleculeAtoms, numberOfTrials, longRange,
      [&](RandomNumber &rng)
      {
        return CBMC::retraceMolecule(rng, context, component, moleculeAtoms,
                                     {.firstBead = FirstBeadScheme::Reinsertion, .storedR = storedR});
      },
      [&](RandomNumber &rng)
      {
        return CBMC::regrowMolecule(rng, context, component, molecule, moleculeAtoms,
                                    {.firstBead = FirstBeadScheme::Reinsertion});
      });
}

const CBMC::MultiTrialGrower &CBMC::cpuMultiTrialGrower()
{
  static const CPUMultiTrialGrower grower{};
  return grower;
}
