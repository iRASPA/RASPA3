module;

module cbmc_recoil_growth;

import std;

import randomnumbers;
import units;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
import double3x3;
import simulationbox;
import energy_status;
import forcefield;
import running_energy;
import framework;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_growth_context;
import cbmc_growth_plan;
import cbmc_operators;

// A generated trial direction: candidate positions of the step's next-beads plus the torsion
// Rosenbluth weight accumulated while selecting the torsion rotation.
using Trial = CBMC::StepTrial;
using Step = CBMC::GrowStep;

// Bundles the (mostly constant) data needed to evaluate energies and grow feelers.
struct RecoilContext
{
  const CBMC::GrowContext &env;
  const Component &component;
  std::make_signed_t<std::size_t> skipBackgroundMolecule;
  std::size_t numberOfTrialDirections;  // k
  std::size_t recoilLength;             // l
  const std::vector<Step> &steps;
  // Per-step openness reference energy (see 'openProbability' below).
  std::vector<double> referenceStepEnergies{};
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
// is ~0 and the standard test is recovered; no per-molecule tuning is needed.
static std::vector<double> computeReferenceStepEnergies(const Component &component, const std::vector<Step> &steps)
{
  std::vector<double> reference(steps.size());
  for (std::size_t seg = 0; seg != steps.size(); ++seg)
  {
    double referenceEnergy;
    if (component.recoilReferenceConformations.empty())
    {
      // No reference conformations built (e.g. a unit test constructing the context directly): fall
      // back to the component's declared geometry.
      referenceEnergy =
          steps[seg].intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(component.atoms).potentialEnergy();
    }
    else
    {
      referenceEnergy = 0.0;
      for (const std::vector<Atom> &conformation : component.recoilReferenceConformations)
      {
        referenceEnergy = std::max(
            referenceEnergy,
            steps[seg].intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(conformation).potentialEnergy());
      }
    }
    reference[seg] = std::max(0.0, referenceEnergy);
  }
  return reference;
}

static double openProbability(const RecoilContext &ctx, std::size_t seg, double potentialEnergy)
{
  return std::min(1.0, std::exp(-ctx.env.beta * (potentialEnergy - ctx.referenceStepEnergies[seg])));
}

// Energy of a trial placement, split into external (non-bonded) and intramolecular vdW/Coulomb.
struct TrialEnergy
{
  RunningEnergy external{};
  RunningEnergy intra{};
  double potentialEnergy() const { return external.potentialEnergy() + intra.potentialEnergy(); }
};

// Compute the energy of placing 'trialPositions' for the step's next-beads. Returns nullopt on a
// hard overlap (a 'closed' direction).
static std::optional<TrialEnergy> computeTrialEnergy(const RecoilContext &ctx, const Step &step,
                                                     const std::vector<Atom> &contextAtoms,
                                                     const std::vector<Atom> &trialPositions)
{
  std::vector<std::vector<Atom>> trialPositionSets{trialPositions};

  std::vector<CBMC::ChainTrial> external = CBMC::computeExternalNonOverlappingEnergies(
      ctx.env, ctx.component, trialPositionSets, -1, ctx.skipBackgroundMolecule);

  if (external.empty()) return std::nullopt;

  std::vector<Atom> candidate(contextAtoms.begin(), contextAtoms.end());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
  {
    candidate[step.nextBeads[k]] = trialPositions[k];
  }
  RunningEnergy intra = step.intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(candidate);

  return TrialEnergy{external.front().energy, intra};
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
static bool feelerExists(RandomNumber &random, const RecoilContext &ctx, std::size_t seg, std::size_t depth,
                         const std::vector<Atom> &contextAtoms)
{
  if (seg >= ctx.steps.size()) return true;
  if (depth == 0) return true;

  const Step &step = ctx.steps[seg];

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections; ++j)
  {
    Trial trial = CBMC::generateRecoilTrial(random, ctx.env.forceField, ctx.env.beta, ctx.component, contextAtoms, step);

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, contextAtoms, trial.positions);
    if (!energy.has_value()) continue;

    double open_probability = openProbability(ctx, seg, energy->potentialEnergy());
    if (random.uniform() < open_probability)
    {
      std::vector<Atom> next_atoms(contextAtoms.begin(), contextAtoms.end());
      for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
      {
        next_atoms[step.nextBeads[k]] = trial.positions[k];
      }
      if (feelerExists(random, ctx, seg + 1, depth - 1, next_atoms)) return true;
    }
  }

  return false;
}

enum class GrowResult { Complete, DeadEnd, Discard };

struct GrowRecord
{
  Trial selected{};
  TrialEnergy energy{};
  double openProbability{1.0};
  std::size_t triedCount{0};
};

static GrowResult growRecursive(RandomNumber &random, const RecoilContext &ctx, std::size_t seg,
                                std::vector<Atom> &atoms, std::size_t &maxHead, std::vector<GrowRecord> &records)
{
  if (seg == ctx.steps.size()) return GrowResult::Complete;

  const Step &step = ctx.steps[seg];

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections; ++j)
  {
    Trial trial = CBMC::generateRecoilTrial(random, ctx.env.forceField, ctx.env.beta, ctx.component, atoms, step);

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, atoms, trial.positions);
    if (!energy.has_value()) continue;

    double open_probability = openProbability(ctx, seg, energy->potentialEnergy());
    if (random.uniform() >= open_probability) continue;

    std::vector<Atom> saved(step.nextBeads.size());
    for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
    {
      saved[k] = atoms[step.nextBeads[k]];
      atoms[step.nextBeads[k]] = trial.positions[k];
    }
    maxHead = std::max(maxHead, seg);

    GrowResult result = growRecursive(random, ctx, seg + 1, atoms, maxHead, records);
    if (result == GrowResult::Complete)
    {
      records[seg] = {trial, energy.value(), open_probability, j + 1};
      return GrowResult::Complete;
    }

    for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
    {
      atoms[step.nextBeads[k]] = saved[k];
    }

    if (result == GrowResult::Discard) return GrowResult::Discard;

    if (maxHead + 1 >= seg + ctx.recoilLength) return GrowResult::Discard;
  }

  return GrowResult::DeadEnd;
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growRecoilGrowthMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced, std::make_signed_t<std::size_t> skipBackgroundMolecule)
{
  const ForceField &forceField = context.forceField;

  const std::vector<Step> &steps = component.growthPlan(beadsAlreadyPlaced);
  RecoilContext ctx{context,
                    component,
                    skipBackgroundMolecule,
                    std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
                    std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength),
                    steps,
                    computeReferenceStepEnergies(component, steps)};

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<GrowRecord> records(ctx.steps.size());
  std::size_t maxHead = 0;

  if (growRecursive(random, ctx, 0, chain_atoms, maxHead, records) != GrowResult::Complete) return std::nullopt;

  double chain_rosenbluth_weight = 1.0;
  double chain_log_rosenbluth_weight = 0.0;
  RunningEnergy chain_external_energies{};

  for (std::size_t seg = 0; seg != ctx.steps.size(); ++seg)
  {
    const Step &step = ctx.steps[seg];
    const GrowRecord &record = records[seg];

    std::size_t numberOfFeelers = 1;
    for (std::size_t j = record.triedCount; j < ctx.numberOfTrialDirections; ++j)
    {
      Trial alternative =
          CBMC::generateRecoilTrial(random, ctx.env.forceField, ctx.env.beta, ctx.component, chain_atoms, step);

      std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, chain_atoms, alternative.positions);
      if (!energy.has_value()) continue;

      double open_probability = openProbability(ctx, seg, energy->potentialEnergy());
      if (random.uniform() >= open_probability) continue;

      std::vector<Atom> feeler_atoms(chain_atoms.begin(), chain_atoms.end());
      for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
      {
        feeler_atoms[step.nextBeads[k]] = alternative.positions[k];
      }

      if (feelerExists(random, ctx, seg + 1, ctx.recoilLength - 1, feeler_atoms)) ++numberOfFeelers;
    }

    double step_weight = static_cast<double>(numberOfFeelers) / static_cast<double>(ctx.numberOfTrialDirections) *
                         std::exp(-ctx.env.beta * record.energy.potentialEnergy()) / record.openProbability *
                         record.selected.torsionWeight;

    chain_external_energies += record.energy.external;

    // Fold this step's not-sampled internal terms into the weight (mirrors the CBMC path). Flexible
    // attach steps handle these terms inside the operator engine, so skip them here (no double count).
    for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
    {
      chain_atoms[step.nextBeads[k]] = record.selected.positions[k];
    }
    if (!CBMC::stepHandlesUnsampledInternalTerms(step))
    {
      RunningEnergy stepUnsampled = step.intra.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms);
      step_weight *= std::exp(-ctx.env.beta * stepUnsampled.potentialEnergy());
    }

    // Per-step overlap guard: the cumulative weight of a long chain is below any fixed absolute
    // threshold (it decays exponentially with chain length), so guarding the running product would
    // reject every grow of a long polymer. Mirrors the CBMC insertion path.
    if (step_weight < forceField.minimumRosenbluthFactor) return std::nullopt;
    chain_rosenbluth_weight *= step_weight;
    // The per-step factor is bounded below by the guard, so its log is finite; the log sum stays exact
    // where the raw product of a long chain underflows to zero.
    chain_log_rosenbluth_weight += std::log(step_weight);
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

  // Center of mass, and for a fully rigid component the orientation quaternion recovered from the
  // grown positions (used downstream to regenerate the atoms of rigid molecules).
  Molecule molecule = component.createMoleculeRecord(chain_atoms);

  return ChainGrowData(molecule, chain_atoms, chain_external_energies + internal_energies, chain_rosenbluth_weight,
                       0.0, chain_log_rosenbluth_weight);
}

[[nodiscard]] ChainRetraceData CBMC::retraceRecoilGrowthMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced) noexcept
{
  const ForceField &forceField = context.forceField;

  const std::vector<Step> &steps = component.growthPlan(beadsAlreadyPlaced);
  RecoilContext ctx{context,
                    component,
                    -1,
                    std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
                    std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength),
                    steps,
                    computeReferenceStepEnergies(component, steps)};

  std::vector<Atom> old_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double chain_rosenbluth_weight = 1.0;
  double chain_log_rosenbluth_weight = 0.0;
  RunningEnergy chain_external_energies{};

  for (std::size_t seg = 0; seg != ctx.steps.size(); ++seg)
  {
    const Step &step = ctx.steps[seg];

    std::vector<Atom> old_positions(step.nextBeads.size());
    for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
    {
      old_positions[k] = old_atoms[step.nextBeads[k]];
    }

    std::optional<TrialEnergy> old_energy = computeTrialEnergy(ctx, step, old_atoms, old_positions);
    TrialEnergy selected_energy = old_energy.value_or(TrialEnergy{});
    double selected_potential = selected_energy.potentialEnergy();
    double open_probability = openProbability(ctx, seg, selected_potential);

    double torsion_weight =
        CBMC::oldConfigurationTorsionWeight(random, ctx.env.forceField, ctx.env.beta, component, old_atoms, step);

    std::size_t numberOfFeelers = 1;
    if (ctx.numberOfTrialDirections > 1)
    {
      for (std::size_t j = 1; j < ctx.numberOfTrialDirections; ++j)
      {
        Trial trial =
            CBMC::generateRecoilTrial(random, ctx.env.forceField, ctx.env.beta, ctx.component, old_atoms, step);

        std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, old_atoms, trial.positions);
        if (!energy.has_value()) continue;

        double alternative_open_probability = openProbability(ctx, seg, energy->potentialEnergy());
        if (random.uniform() >= alternative_open_probability) continue;

        std::vector<Atom> next_atoms(old_atoms.begin(), old_atoms.end());
        for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
        {
          next_atoms[step.nextBeads[k]] = trial.positions[k];
        }

        if (feelerExists(random, ctx, seg + 1, ctx.recoilLength - 1, next_atoms)) ++numberOfFeelers;
      }
    }

    chain_rosenbluth_weight *= static_cast<double>(numberOfFeelers) /
                                static_cast<double>(ctx.numberOfTrialDirections) *
                                std::exp(-ctx.env.beta * selected_potential) / open_probability * torsion_weight;

    chain_external_energies += selected_energy.external;

    double stepUnsampledEnergy = 0.0;
    if (!CBMC::stepHandlesUnsampledInternalTerms(step))
    {
      stepUnsampledEnergy = step.intra.computeInternalEnergiesNotSampledDuringGrowth(old_atoms).potentialEnergy();
    }
    chain_rosenbluth_weight *= std::exp(-ctx.env.beta * stepUnsampledEnergy);

    // Log of the same per-step factor (retrace has no per-step guard, so the raw product of a long
    // chain underflows; the log sum stays exact).
    chain_log_rosenbluth_weight +=
        std::log(static_cast<double>(numberOfFeelers) / static_cast<double>(ctx.numberOfTrialDirections)) -
        ctx.env.beta * selected_potential - std::log(open_probability) + std::log(torsion_weight) -
        ctx.env.beta * stepUnsampledEnergy;
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(old_atoms);

  return ChainRetraceData(chain_external_energies + internal_energies, chain_rosenbluth_weight, 0.0,
                          chain_log_rosenbluth_weight);
}
