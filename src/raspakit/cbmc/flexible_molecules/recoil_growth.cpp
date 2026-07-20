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
};

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

    double open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
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

    double open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
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

  RecoilContext ctx{context,
                    component,
                    skipBackgroundMolecule,
                    std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
                    std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength),
                    component.growthPlan(beadsAlreadyPlaced)};

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<GrowRecord> records(ctx.steps.size());
  std::size_t maxHead = 0;

  if (growRecursive(random, ctx, 0, chain_atoms, maxHead, records) != GrowResult::Complete) return std::nullopt;

  double chain_rosenbluth_weight = 1.0;
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

      double open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
      if (random.uniform() >= open_probability) continue;

      std::vector<Atom> feeler_atoms(chain_atoms.begin(), chain_atoms.end());
      for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
      {
        feeler_atoms[step.nextBeads[k]] = alternative.positions[k];
      }

      if (feelerExists(random, ctx, seg + 1, ctx.recoilLength - 1, feeler_atoms)) ++numberOfFeelers;
    }

    chain_rosenbluth_weight *= static_cast<double>(numberOfFeelers) /
                                static_cast<double>(ctx.numberOfTrialDirections) *
                                std::exp(-ctx.env.beta * record.energy.potentialEnergy()) / record.openProbability *
                                record.selected.torsionWeight;

    chain_external_energies += record.energy.external;

    // Fold this step's not-sampled internal terms into the weight (mirrors the CBMC path).
    for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
    {
      chain_atoms[step.nextBeads[k]] = record.selected.positions[k];
    }
    RunningEnergy stepUnsampled = step.intra.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms);
    chain_rosenbluth_weight *= std::exp(-ctx.env.beta * stepUnsampled.potentialEnergy());

    if (chain_rosenbluth_weight < forceField.minimumRosenbluthFactor) return std::nullopt;
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

  // Keep this thermalized, non-overlapping conformation as the warm start for the next grow.
  component.warmStartConformation = chain_atoms;

  // Center of mass, and for a fully rigid component the orientation quaternion recovered from the
  // grown positions (used downstream to regenerate the atoms of rigid molecules).
  Molecule molecule = component.createMoleculeRecord(chain_atoms);

  return ChainGrowData(molecule, chain_atoms, chain_external_energies + internal_energies, chain_rosenbluth_weight,
                       0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceRecoilGrowthMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced) noexcept
{
  const ForceField &forceField = context.forceField;

  RecoilContext ctx{context,
                    component,
                    -1,
                    std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
                    std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength),
                    component.growthPlan(beadsAlreadyPlaced)};

  std::vector<Atom> old_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double chain_rosenbluth_weight = 1.0;
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
    double open_probability = std::min(1.0, std::exp(-ctx.env.beta * selected_potential));

    double torsion_weight =
        CBMC::oldConfigurationTorsionWeight(random, ctx.env.forceField, ctx.env.beta, old_atoms, step);

    std::size_t numberOfFeelers = 1;
    if (ctx.numberOfTrialDirections > 1)
    {
      for (std::size_t j = 1; j < ctx.numberOfTrialDirections; ++j)
      {
        Trial trial =
            CBMC::generateRecoilTrial(random, ctx.env.forceField, ctx.env.beta, ctx.component, old_atoms, step);

        std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, step, old_atoms, trial.positions);
        if (!energy.has_value()) continue;

        double alternative_open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
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

    RunningEnergy stepUnsampled = step.intra.computeInternalEnergiesNotSampledDuringGrowth(old_atoms);
    chain_rosenbluth_weight *= std::exp(-ctx.env.beta * stepUnsampled.potentialEnergy());
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(old_atoms);

  return ChainRetraceData(chain_external_energies + internal_energies, chain_rosenbluth_weight, 0.0);
}
