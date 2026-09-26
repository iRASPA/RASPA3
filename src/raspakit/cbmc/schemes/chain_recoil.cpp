module;

module cbmc_chain_recoil;

import std;

import randomnumbers;
import component;
import atom;
import running_energy;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_results;
import cbmc_grow_context;
import cbmc_grow_step;
import cbmc_operators;
import cbmc_chain_common;
import cbmc_statistics;

namespace
{
using namespace CBMC;

// What every level of the recoil search reads: the growth environment, the component, the plan, and
// the per-step openness reference. The trial count 'k' and recoil length 'l' are the context's
// settings ('recoilGrowthNumberOfTrialDirections', 'recoilGrowthMaximumRecoilLength').
//
// Trial energies are 'evaluateStepTrial' of the chain schemes (nullopt on a hard overlap: a 'closed'
// direction). The chain is edited in place throughout the recursion (a few beads per level, guarded
// by CBMC::ScratchBeads) rather than copied per trial: for a polymer of N beads with k trial directions
// and recoil length l, per-trial chain copies made every grow O(N^2).
struct RecoilContext
{
  const GrowContext &env;
  const Component &component;
  const std::vector<GrowStep> &steps;
  // Per-step openness reference energy (see 'openProbability' below); cached per plan by the
  // component, alongside the plan itself.
  const std::vector<double> &referenceStepEnergies;

  [[nodiscard]] std::size_t numberOfTrialDirections() const noexcept
  {
    return env.settings.recoilGrowthNumberOfTrialDirections;
  }
  [[nodiscard]] std::size_t recoilLength() const noexcept { return env.settings.recoilGrowthMaximumRecoilLength; }
};

// The open/closed test of recoil growth is an absolute Boltzmann filter, and the weight divides by the
// openness probability of the selected direction, so any fixed per-step energy reference is formally
// valid (grow and retrace use the same function). The reference matters in practice: a molecule whose
// correctly grown chains carry systematic positive non-bonded strain per placed bead — crowded
// united-atom beads at branch points, 1-4/1-5 intramolecular Coulomb between partial charges; a few
// hundred to a few thousand kelvin — has most correctly-placed beads test 'closed' against a zero
// reference, and the recoil search backtracks essentially forever. Measuring each step against the same
// step evaluated in the component's recoil reference conformations (equilibrated ideal-gas conformations
// built once at setup) removes exactly the molecule's own intrinsic strain, while a genuine overlap
// (1e6 K) still tests closed. Per step, the reference is the MAXIMUM energy over the conformations:
// every reference conformation is a valid equilibrated chain, so any strain it exhibits at a step is by
// definition acceptable there, and a placement at that level must test open with probability one. (The
// per-step minimum would be the floor of the strain distribution: typical good placements sit above it,
// each step then tests open with probability well below one, and the compounded attrition over the
// hundreds of steps of a polymer aborts essentially every grow.) For an unstrained chain the reference
// is ~0 and the standard test is recovered; no per-molecule tuning is needed. The reference is a pure
// function of the plan and the reference conformations, so the component computes and caches it per
// plan ('Component::recoilReferenceStepEnergies').
double openProbability(const RecoilContext &ctx, std::size_t seg, double potentialEnergy)
{
  return std::min(1.0, std::exp(-ctx.env.beta * (potentialEnergy - ctx.referenceStepEnergies[seg])));
}

// Test whether an open pathway ('feeler') of 'depth' more steps can be grown starting at step 'seg'.
//
// A trial direction is 'available' when it is open and a feeler of 'recoilLength - 1' steps can be
// grown from it. The growth ('growRecursive') decides the availability of the directions it tries by
// attempting to grow the chain itself: a direction that is open but whose sub-tree dead-ends before
// reaching 'recoilLength - 1' steps ahead is unavailable, exactly the negation of this feeler. Both
// tests are exhaustive depth-first searches over the same 'numberOfTrialDirections' per level with the
// same openness test, and they MUST draw their trial beads from the same generator ('generateRecoilTrial'
// with its torsion-selected spin): the number of available directions m_i enters the acceptance ratio,
// and the retrace decides every alternative with this feeler while the grow decides the tried-and-failed
// directions with the growth attempt. If the two probes sampled the spin differently, the probability of
// finding a direction available would differ between the two, and m_i would be biased between grow and
// retrace.
//
// 'atoms' is scratch: the feeler places its trial beads in place and restores them, so the chain is
// unchanged on return.
bool feelerExists(RandomNumber &random, const RecoilContext &ctx, std::size_t seg, std::size_t depth,
                  std::vector<Atom> &atoms)
{
  if (seg >= ctx.steps.size()) return true;
  if (depth == 0) return true;

  const GrowStep &step = ctx.steps[seg];
  const ScratchBeads scratch(atoms, step);

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections(); ++j)
  {
    StepTrial trial = generateRecoilTrial(random, ctx.env.settings, ctx.env.beta, ctx.component, atoms, step);

    std::optional<StepTrialEnergy> energy = evaluateStepTrial(ctx.env, ctx.component, step, atoms, trial.positions);
    if (!energy.has_value()) continue;

    if (random.uniform() < openProbability(ctx, seg, energy->potentialEnergy()))
    {
      placeStepBeads(atoms, step, trial.positions);
      const bool found = feelerExists(random, ctx, seg + 1, depth - 1, atoms);
      scratch.restore();
      if (found) return true;
    }
  }

  return false;
}

// Whether a trial direction that tested open at step 'seg' is 'available': a feeler of the remaining
// recoil length can be grown from it. Places the trial beads in 'atoms' for the feeler and restores
// them afterwards.
bool feelerExistsFrom(RandomNumber &random, const RecoilContext &ctx, std::size_t seg, std::vector<Atom> &atoms,
                      std::span<const Atom> trialPositions)
{
  const GrowStep &step = ctx.steps[seg];
  const ScratchBeads scratch(atoms, step);
  placeStepBeads(atoms, step, trialPositions);
  return feelerExists(random, ctx, seg + 1, ctx.recoilLength() - 1, atoms);
}

// The number of available directions m_i of step 'seg': the selected direction (always counted) plus
// those of the freshly generated alternatives 'firstAlternative' .. k-1 that are open and have a
// feeler. On the grow the alternatives are the directions the growth never tried ('triedCount' ..);
// on the retrace the old configuration is direction 0 and all others (1 ..) are alternatives. Both
// probe with the same generator and feeler, so m_i is the same random experiment in either direction.
// 'atoms' holds the chain with the selected direction in place at 'seg'; the feelers use the beads
// beyond it as scratch and restore them.
std::size_t countAvailableDirections(RandomNumber &random, const RecoilContext &ctx, std::size_t seg,
                                     std::vector<Atom> &atoms, std::size_t firstAlternative)
{
  const GrowStep &step = ctx.steps[seg];
  std::size_t available = 1;
  for (std::size_t j = firstAlternative; j < ctx.numberOfTrialDirections(); ++j)
  {
    StepTrial alternative = generateRecoilTrial(random, ctx.env.settings, ctx.env.beta, ctx.component, atoms, step);

    std::optional<StepTrialEnergy> energy = evaluateStepTrial(ctx.env, ctx.component, step, atoms, alternative.positions);
    if (!energy.has_value()) continue;

    if (random.uniform() >= openProbability(ctx, seg, energy->potentialEnergy())) continue;

    if (feelerExistsFrom(random, ctx, seg, atoms, alternative.positions)) ++available;
  }
  return available;
}

// Log of the per-step recoil factor m_i / k * exp(-beta u_i) / p_open * w_torsion * exp(-beta u_unsampled),
// evaluated with the step's beads in place in 'atoms' for the not-sampled internal terms.
double stepLogWeight(const RecoilContext &ctx, std::size_t seg, std::size_t available, double potentialEnergy,
                     double openProbabilityOfSelected, double torsionWeight, const std::vector<Atom> &atoms)
{
  return std::log(static_cast<double>(available) / static_cast<double>(ctx.numberOfTrialDirections())) -
         ctx.env.beta * potentialEnergy - std::log(openProbabilityOfSelected) + std::log(torsionWeight) +
         unsampledInternalLogFactor(ctx.env.beta, ctx.steps[seg], atoms);
}

enum class GrowOutcome
{
  Complete,
  DeadEnd,
  Discard
};

// What the growth records per step once the chain is complete: the selected direction, its energy and
// openness probability, and how many directions were tried at that step (the last one is the selected).
struct GrowRecord
{
  StepTrial selected{};
  StepTrialEnergy energy{};
  double openProbability{1.0};
  std::size_t triedCount{0};
};

GrowOutcome growRecursive(RandomNumber &random, const RecoilContext &ctx, std::size_t seg, std::vector<Atom> &atoms,
                          std::size_t &maxHead, std::vector<GrowRecord> &records)
{
  if (seg == ctx.steps.size()) return GrowOutcome::Complete;

  const GrowStep &step = ctx.steps[seg];

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections(); ++j)
  {
    StepTrial trial = generateRecoilTrial(random, ctx.env.settings, ctx.env.beta, ctx.component, atoms, step);

    std::optional<StepTrialEnergy> energy = evaluateStepTrial(ctx.env, ctx.component, step, atoms, trial.positions);
    if (!energy.has_value()) continue;

    const double open = openProbability(ctx, seg, energy->potentialEnergy());
    if (random.uniform() >= open) continue;

    // The placed beads stay in the chain when the growth below completes (the chain then holds the
    // complete grown molecule); otherwise the guard restores them for the next direction.
    ScratchBeads scratch(atoms, step);
    placeStepBeads(atoms, step, trial.positions);
    maxHead = std::max(maxHead, seg);

    const GrowOutcome result = growRecursive(random, ctx, seg + 1, atoms, maxHead, records);
    if (result == GrowOutcome::Complete)
    {
      scratch.keep();
      records[seg] = {std::move(trial), energy.value(), open, j + 1};
      return GrowOutcome::Complete;
    }

    scratch.restore();

    if (result == GrowOutcome::Discard) return GrowOutcome::Discard;

    // The sub-tree reached 'recoilLength - 1' steps ahead, so this direction was available and the
    // growth committed to it; recoiling past a committed direction would break the scheme.
    if (maxHead + 1 >= seg + ctx.recoilLength()) return GrowOutcome::Discard;
  }

  return GrowOutcome::DeadEnd;
}

RecoilContext makeRecoilContext(const GrowContext &context, const Component &component,
                                const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  return RecoilContext{context, component, component.growthPlan(beadsAlreadyPlaced),
                       component.recoilReferenceStepEnergies(beadsAlreadyPlaced)};
}
}  // namespace

[[nodiscard]] std::optional<CBMC::GrowResult> CBMC::growChainRecoil(RandomNumber &random, const GrowContext &context,
                                                                    const Component &component,
                                                                    std::span<const Atom> moleculeAtoms,
                                                                    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const RecoilContext ctx = makeRecoilContext(context, component, beadsAlreadyPlaced);

  std::vector<Atom> chainAtoms(moleculeAtoms.begin(), moleculeAtoms.end());
  std::vector<GrowRecord> records(ctx.steps.size());
  std::size_t maxHead = 0;

  RecoilGrowthStatistics &statistics = component.recoilGrowthStatistics;
  statistics.grows += 1.0;

  switch (growRecursive(random, ctx, 0, chainAtoms, maxHead, records))
  {
    case GrowOutcome::DeadEnd:
      statistics.deadEnds += 1.0;
      return std::nullopt;
    case GrowOutcome::Discard:
      statistics.discarded += 1.0;
      return std::nullopt;
    case GrowOutcome::Complete:
      statistics.completed += 1.0;
      break;
  }

  // 'chainAtoms' now holds the complete grown molecule; the alternatives of every step are probed
  // a posteriori against it (the feelers use the beads beyond the step as scratch).
  ChainAccumulator chain{};
  for (std::size_t seg = 0; seg != ctx.steps.size(); ++seg)
  {
    const GrowRecord &record = records[seg];

    // The directions tried and failed by the growth are unavailable by construction; only the
    // untried ones are probed.
    const std::size_t available = countAvailableDirections(random, ctx, seg, chainAtoms, record.triedCount);
    statistics.availableDirectionsSum += static_cast<double>(available);
    statistics.growSteps += 1.0;

    const double logWeight = stepLogWeight(ctx, seg, available, record.energy.potentialEnergy(),
                                           record.openProbability, record.selected.torsionWeight, chainAtoms);
    if (!chain.addGrownStep(logWeight, record.energy.external, context.settings.minimumRosenbluthFactor))
    {
      return std::nullopt;
    }
  }

  return finishGrownChain(component, std::move(chainAtoms), chain);
}

[[nodiscard]] CBMC::RetraceResult CBMC::retraceChainRecoil(RandomNumber &random, const GrowContext &context,
                                                           const Component &component,
                                                           std::span<const Atom> moleculeAtoms,
                                                           const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const RecoilContext ctx = makeRecoilContext(context, component, beadsAlreadyPlaced);

  std::vector<Atom> oldAtoms(moleculeAtoms.begin(), moleculeAtoms.end());
  ChainAccumulator chain{};

  for (std::size_t seg = 0; seg != ctx.steps.size(); ++seg)
  {
    const GrowStep &step = ctx.steps[seg];

    // The old configuration must not overlap: the recoil weight divides by its openness probability,
    // which would be zero, so no weight is defined for it.
    const std::optional<StepTrialEnergy> oldEnergy =
        evaluateStepTrial(ctx.env, ctx.component, step, oldAtoms, stepBeadPositions(oldAtoms, step));
    if (!oldEnergy.has_value()) throwExistingConfigurationOverlaps("Recoil growth", component, seg, step);
    const double oldPotential = oldEnergy->potentialEnergy();

    const double torsionWeight =
        oldConfigurationTorsionWeight(random, ctx.env.settings, ctx.env.beta, component, oldAtoms, step);

    // The old configuration is trial direction 0 (always counted available); the remaining k - 1
    // directions are generated and probed exactly as on the grow. No per-step guard on the retrace.
    const std::size_t available = countAvailableDirections(random, ctx, seg, oldAtoms, 1);

    chain.addRetracedStep(stepLogWeight(ctx, seg, available, oldPotential, openProbability(ctx, seg, oldPotential),
                                        torsionWeight, oldAtoms),
                          oldEnergy->external);
  }

  return finishRetracedChain(component, moleculeAtoms, chain);
}
