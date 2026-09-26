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

// A generated trial direction: candidate positions of the step's next-beads plus the torsion
// Rosenbluth weight accumulated while selecting the torsion rotation.
using Trial = CBMC::StepTrial;
using Step = CBMC::GrowStep;
using TrialEnergy = CBMC::StepTrialEnergy;

// Bundles the (mostly constant) data needed to evaluate energies and grow feelers.
struct RecoilContext
{
  const CBMC::GrowContext &env;
  const Component &component;
  std::size_t numberOfTrialDirections;  // k
  std::size_t recoilLength;             // l
  const std::vector<Step> &steps;
  // Per-step openness reference energy (see 'openProbability' below); cached per plan by the
  // component, alongside the plan itself.
  const std::vector<double> &referenceStepEnergies;
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
static double openProbability(const RecoilContext &ctx, std::size_t seg, double potentialEnergy)
{
  return std::min(1.0, std::exp(-ctx.env.beta * (potentialEnergy - ctx.referenceStepEnergies[seg])));
}

// Energy of placing 'trialPositions' for the step's next-beads: the shared trial evaluation of the
// chain schemes. Returns nullopt on a hard overlap (a 'closed' direction). 'atoms' is used as scratch
// and left unchanged. The chain is edited in place throughout the recursion (a few beads per level,
// guarded by CBMC::ScratchBeads) rather than copied per trial: for a polymer of N beads with k trial
// directions and recoil length l, per-trial chain copies made every grow O(N^2).
static std::optional<TrialEnergy> computeTrialEnergy(const RecoilContext &ctx, const Step &step,
                                                     std::vector<Atom> &atoms, std::span<const Atom> trialPositions)
{
  return CBMC::evaluateStepTrial(ctx.env, ctx.component, step, atoms, trialPositions);
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
static bool feelerExists(RandomNumber &random, const RecoilContext &ctx, std::size_t seg, std::size_t depth,
                         std::vector<Atom> &atoms)
{
  if (seg >= ctx.steps.size()) return true;
  if (depth == 0) return true;

  const Step &step = ctx.steps[seg];
  const CBMC::ScratchBeads scratch(atoms, step);

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections; ++j)
  {
    Trial trial = CBMC::generateRecoilTrial(random, ctx.env.settings, ctx.env.beta, ctx.component, atoms, step);

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, atoms, trial.positions);
    if (!energy.has_value()) continue;

    double open_probability = openProbability(ctx, seg, energy->potentialEnergy());
    if (random.uniform() < open_probability)
    {
      CBMC::placeStepBeads(atoms, step, trial.positions);
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
static bool feelerExistsFrom(RandomNumber &random, const RecoilContext &ctx, std::size_t seg,
                             std::vector<Atom> &atoms, std::span<const Atom> trialPositions)
{
  const Step &step = ctx.steps[seg];
  const CBMC::ScratchBeads scratch(atoms, step);
  CBMC::placeStepBeads(atoms, step, trialPositions);
  return feelerExists(random, ctx, seg + 1, ctx.recoilLength - 1, atoms);
}

// The number of available directions m_i of step 'seg': the selected direction (always counted) plus
// those of the freshly generated alternatives 'firstAlternative' .. k-1 that are open and have a
// feeler. On the grow the alternatives are the directions the growth never tried ('triedCount' ..);
// on the retrace the old configuration is direction 0 and all others (1 ..) are alternatives. Both
// probe with the same generator and feeler, so m_i is the same random experiment in either direction.
// 'atoms' holds the chain with the selected direction in place at 'seg'; the feelers use the beads
// beyond it as scratch and restore them.
static std::size_t countAvailableDirections(RandomNumber &random, const RecoilContext &ctx, std::size_t seg,
                                            std::vector<Atom> &atoms, std::size_t firstAlternative)
{
  const Step &step = ctx.steps[seg];
  std::size_t available = 1;
  for (std::size_t j = firstAlternative; j < ctx.numberOfTrialDirections; ++j)
  {
    Trial alternative = CBMC::generateRecoilTrial(random, ctx.env.settings, ctx.env.beta, ctx.component, atoms, step);

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, atoms, alternative.positions);
    if (!energy.has_value()) continue;

    if (random.uniform() >= openProbability(ctx, seg, energy->potentialEnergy())) continue;

    if (feelerExistsFrom(random, ctx, seg, atoms, alternative.positions)) ++available;
  }
  return available;
}

// Log of the per-step recoil factor m_i / k * exp(-beta u_i) / p_open * w_torsion * exp(-beta u_unsampled),
// evaluated with the step's beads in place in 'atoms' for the not-sampled internal terms.
static double stepLogWeight(const RecoilContext &ctx, std::size_t seg, std::size_t available, double potentialEnergy,
                            double open_probability, double torsionWeight, const std::vector<Atom> &atoms)
{
  return std::log(static_cast<double>(available) / static_cast<double>(ctx.numberOfTrialDirections)) -
         ctx.env.beta * potentialEnergy - std::log(open_probability) + std::log(torsionWeight) +
         CBMC::unsampledInternalLogFactor(ctx.env.beta, ctx.steps[seg], atoms);
}

enum class GrowOutcome { Complete, DeadEnd, Discard };

struct GrowRecord
{
  Trial selected{};
  TrialEnergy energy{};
  double openProbability{1.0};
  std::size_t triedCount{0};
};

static GrowOutcome growRecursive(RandomNumber &random, const RecoilContext &ctx, std::size_t seg,
                                std::vector<Atom> &atoms, std::size_t &maxHead, std::vector<GrowRecord> &records)
{
  if (seg == ctx.steps.size()) return GrowOutcome::Complete;

  const Step &step = ctx.steps[seg];

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections; ++j)
  {
    Trial trial = CBMC::generateRecoilTrial(random, ctx.env.settings, ctx.env.beta, ctx.component, atoms, step);

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, atoms, trial.positions);
    if (!energy.has_value()) continue;

    double open_probability = openProbability(ctx, seg, energy->potentialEnergy());
    if (random.uniform() >= open_probability) continue;

    // The placed beads stay in the chain when the growth below completes (the chain then holds the
    // complete grown molecule); otherwise the guard restores them for the next direction.
    CBMC::ScratchBeads scratch(atoms, step);
    CBMC::placeStepBeads(atoms, step, trial.positions);
    maxHead = std::max(maxHead, seg);

    GrowOutcome result = growRecursive(random, ctx, seg + 1, atoms, maxHead, records);
    if (result == GrowOutcome::Complete)
    {
      scratch.keep();
      records[seg] = {std::move(trial), energy.value(), open_probability, j + 1};
      return GrowOutcome::Complete;
    }

    scratch.restore();

    if (result == GrowOutcome::Discard) return GrowOutcome::Discard;

    if (maxHead + 1 >= seg + ctx.recoilLength) return GrowOutcome::Discard;
  }

  return GrowOutcome::DeadEnd;
}

static RecoilContext makeRecoilContext(const CBMC::GrowContext &context, const Component &component,
                                       const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const CBMC::GrowthSettings &settings = context.settings;
  return RecoilContext{context,
                       component,
                       settings.recoilGrowthNumberOfTrialDirections,
                       settings.recoilGrowthMaximumRecoilLength,
                       component.growthPlan(beadsAlreadyPlaced),
                       component.recoilReferenceStepEnergies(beadsAlreadyPlaced)};
}

[[nodiscard]] std::optional<CBMC::GrowResult> CBMC::growChainRecoil(RandomNumber &random, const GrowContext &context,
                                                                    const Component &component,
                                                                    std::span<const Atom> molecule_atoms,
                                                                    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const RecoilContext ctx = makeRecoilContext(context, component, beadsAlreadyPlaced);

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<GrowRecord> records(ctx.steps.size());
  std::size_t maxHead = 0;

  if (growRecursive(random, ctx, 0, chain_atoms, maxHead, records) != GrowOutcome::Complete) return std::nullopt;

  // 'chain_atoms' now holds the complete grown molecule; the alternatives of every step are probed
  // a posteriori against it (the feelers use the beads beyond the step as scratch).
  ChainAccumulator chain{};
  for (std::size_t seg = 0; seg != ctx.steps.size(); ++seg)
  {
    const GrowRecord &record = records[seg];

    // The directions tried and failed by the growth are unavailable by construction; only the
    // untried ones are probed.
    const std::size_t available = countAvailableDirections(random, ctx, seg, chain_atoms, record.triedCount);

    const double step_log_weight = stepLogWeight(ctx, seg, available, record.energy.potentialEnergy(),
                                                 record.openProbability, record.selected.torsionWeight, chain_atoms);
    if (!chain.addGrownStep(step_log_weight, record.energy.external, context.settings.minimumRosenbluthFactor))
    {
      return std::nullopt;
    }
  }

  return finishGrownChain(component, std::move(chain_atoms), chain);
}

[[nodiscard]] CBMC::RetraceResult CBMC::retraceChainRecoil(RandomNumber &random, const GrowContext &context,
                                                           const Component &component,
                                                           std::span<const Atom> molecule_atoms,
                                                           const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const RecoilContext ctx = makeRecoilContext(context, component, beadsAlreadyPlaced);

  std::vector<Atom> old_atoms(molecule_atoms.begin(), molecule_atoms.end());
  ChainAccumulator chain{};

  for (std::size_t seg = 0; seg != ctx.steps.size(); ++seg)
  {
    const Step &step = ctx.steps[seg];

    // The old configuration must not overlap: the recoil weight divides by its openness probability,
    // which would be zero, so no weight is defined for it.
    const std::optional<TrialEnergy> old_energy =
        computeTrialEnergy(ctx, step, old_atoms, CBMC::stepBeadPositions(old_atoms, step));
    if (!old_energy.has_value()) throwExistingConfigurationOverlaps("Recoil growth", component, seg, step);
    const double old_potential = old_energy->potentialEnergy();

    const double torsion_weight =
        CBMC::oldConfigurationTorsionWeight(random, ctx.env.settings, ctx.env.beta, component, old_atoms, step);

    // The old configuration is trial direction 0 (always counted available); the remaining k - 1
    // directions are generated and probed exactly as on the grow. No per-step guard on the retrace.
    const std::size_t available = countAvailableDirections(random, ctx, seg, old_atoms, 1);

    chain.addRetracedStep(stepLogWeight(ctx, seg, available, old_potential, openProbability(ctx, seg, old_potential),
                                        torsion_weight, old_atoms),
                          old_energy->external);
  }

  return finishRetracedChain(component, molecule_atoms, chain);
}
