module;

module cbmc_chain_cbmc;

import std;

import randomnumbers;
import component;
import molecule;
import atom;
import forcefield;
import running_energy;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_results;
import cbmc_growth_context;
import cbmc_growth_plan;
import cbmc_operators;
import cbmc_external_energy;

namespace
{
/// A non-overlapping trial direction of a step: its positions, its external energy, and the torsion
/// Rosenbluth weight accumulated while generating it.
struct ChainTrialTorsion
{
  std::vector<Atom> positions;
  RunningEnergy energy;
  double torsionWeight;
};

// The trial directions of a step against the background: the operator engine's trials, filtered to the
// non-overlapping ones and paired with their external energies. 'firstSurvived' records explicitly
// whether trial direction 0 passed the filter (the retrace needs to know that its old configuration,
// which it pins as trial 0, is still among the survivors).
struct EvaluatedTrials
{
  std::vector<ChainTrialTorsion> trials{};
  bool firstSurvived{false};
};

EvaluatedTrials externalEnergiesOfTrials(const CBMC::GrowContext &context, const Component &component,
                                         std::vector<CBMC::StepTrial> stepTrials,
                                         std::optional<std::size_t> skipBackgroundMolecule)
{
  EvaluatedTrials evaluated{};
  evaluated.trials.reserve(stepTrials.size());
  for (std::size_t i = 0; i != stepTrials.size(); ++i)
  {
    std::optional<RunningEnergy> energy = CBMC::computeExternalNonOverlappingEnergy(
        context, component, stepTrials[i].positions, skipBackgroundMolecule);
    if (!energy.has_value()) continue;
    if (i == 0) evaluated.firstSurvived = true;
    evaluated.trials.push_back({std::move(stepTrials[i].positions), energy.value(), stepTrials[i].torsionWeight});
  }
  return evaluated;
}

struct StepWeight
{
  std::size_t selected;  ///< Index into the step's non-overlapping trials.
  double logWeight;      ///< Log of the step's Rosenbluth factor.
};

// The per-step Rosenbluth factor shared by grow and retrace,
//
//   w_torsion(selected) * sum_j exp(-beta u_j) / k * exp(-beta u_unsampled),
//
// with u_j the external plus the intramolecular van der Waals and Coulomb energy of trial j (each
// trial placed in 'chainAtoms' for that evaluation) and k the number of trial directions. On the grow
// the selected trial is drawn with probability exp(-beta u_j) / sum; on the retrace it is trial 0, the
// old configuration. The selected trial's positions are left in 'chainAtoms'.
//
// The unsampled internal terms (Urey-Bradley, inversion/out-of-plane bend, improper torsion, cross
// terms) of flexible attach steps are handled inside the operator engine (base-coupling rejection plus
// the torsion-spin selection); rigid-body, ring-closure, and seed steps fold them into this factor.
//
// The Boltzmann sum is evaluated as log-sum-exp so the log stays exact even where the raw factor or
// the running product of a long chain would underflow.
StepWeight stepWeight(RandomNumber &random, double beta, std::size_t numberOfTrialDirections,
                      const CBMC::GrowStep &step, std::vector<Atom> &chainAtoms,
                      const std::vector<ChainTrialTorsion> &trials, bool retrace)
{
  const std::vector<std::size_t> &nextBeads = step.nextBeads;

  std::vector<double> logBoltzmannFactors{};
  logBoltzmannFactors.reserve(trials.size());
  for (const ChainTrialTorsion &trial : trials)
  {
    for (std::size_t k = 0; k != nextBeads.size(); ++k) chainAtoms[nextBeads[k]] = trial.positions[k];
    RunningEnergy intraEnergy = step.intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(chainAtoms);
    logBoltzmannFactors.push_back(-beta * (trial.energy.potentialEnergy() + intraEnergy.potentialEnergy()));
  }

  const std::size_t selected = retrace ? 0 : CBMC::selectTrialPosition(random, logBoltzmannFactors);
  const ChainTrialTorsion &selectedTrial = trials[selected];
  for (std::size_t k = 0; k != nextBeads.size(); ++k) chainAtoms[nextBeads[k]] = selectedTrial.positions[k];

  double unsampledEnergy = 0.0;
  if (!CBMC::stepHandlesUnsampledInternalTerms(step))
  {
    unsampledEnergy = step.intra.computeInternalEnergiesNotSampledDuringGrowth(chainAtoms).potentialEnergy();
  }

  const double maxLogBoltzmannFactor = *std::max_element(logBoltzmannFactors.begin(), logBoltzmannFactors.end());
  const double logRosenbluthSum =
      maxLogBoltzmannFactor +
      std::log(std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                               [&](const double &acc, const double &logBoltzmannFactor)
                               { return acc + std::exp(logBoltzmannFactor - maxLogBoltzmannFactor); }));

  return {selected, std::log(selectedTrial.torsionWeight) + logRosenbluthSum -
                        std::log(static_cast<double>(numberOfTrialDirections)) - beta * unsampledEnergy};
}
}  // namespace

[[nodiscard]] std::optional<ChainGrowData> CBMC::growFlexibleMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced, std::optional<std::size_t> skipBackgroundMolecule)
{
  const ForceField &forceField = context.forceField;
  const double beta = context.beta;

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double chain_log_rosenbluth_weight = 0.0;
  RunningEnergy chain_external_energies{};

  // Deterministic growth plan over the fragment graph (flexible beads, hinged rigid bodies, and
  // ring-closure of cyclic clusters), shared with the retrace so grow and retrace are reversible.
  const std::vector<GrowStep> &plan = component.growthPlan(beadsAlreadyPlaced);

  for (const GrowStep &step : plan)
  {
    // All trial directions of this step (the operator engine handles the seed / attach / ring-closure
    // cases, the rigid-body tilt, and the coupled-decoupled torsion selection).
    std::vector<ChainTrialTorsion> trials =
        externalEnergiesOfTrials(context, component,
                                 generateGrowTrials(random, forceField, beta, component, chain_atoms, step,
                                                    forceField.numberOfTrialDirections),
                                 skipBackgroundMolecule)
            .trials;
    if (trials.empty()) return std::nullopt;

    const StepWeight weight =
        stepWeight(random, beta, forceField.numberOfTrialDirections, step, chain_atoms, trials, false);

    chain_external_energies += trials[weight.selected].energy;

    // Overlap guard on the per-step factor, not the running product. The cumulative weight of a long
    // chain decays roughly as f^N (f ~ 0.1 per bead for a chain with intra 1-4 charges), so a chain of
    // a few hundred beads is below any fixed absolute threshold in every attempt: growth would always
    // "dead-end" even though no step overlaps. A genuine overlap still trips the per-step test, since a
    // surviving trial near the 'energyOverlapCriteria' contributes exp(-beta*E) << the threshold. The
    // retrace path carries no guard, so grow and retrace stay symmetric.
    if (std::exp(weight.logWeight) < forceField.minimumRosenbluthFactor) return std::nullopt;
    chain_log_rosenbluth_weight += weight.logWeight;
  }

  // Recompute all the internal interactions (including the cross-terms) for the returned energy.
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

  // Build a valid molecule record: center of mass, and for a fully rigid component the orientation
  // quaternion recovered from the grown positions (used by translation/rotation moves and rigid-body
  // molecular dynamics to regenerate the atoms).
  Molecule molecule = component.createMoleculeRecord(chain_atoms);

  return ChainGrowData(molecule, std::move(chain_atoms), chain_external_energies + internal_energies,
                       chain_log_rosenbluth_weight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceFlexibleMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const ForceField &forceField = context.forceField;
  const double beta = context.beta;

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double chain_log_rosenbluth_weight = 0.0;
  RunningEnergy chain_external_energies{};

  // Same deterministic growth plan as the insertion so grow and retrace are exactly reversible.
  const std::vector<GrowStep> &plan = component.growthPlan(beadsAlreadyPlaced);

  for (std::size_t seg = 0; seg != plan.size(); ++seg)
  {
    const GrowStep &step = plan[seg];

    // The old positions of this step's beads are trial direction 0 of the retrace (read from
    // 'chain_atoms', which holds the old configuration of every bead at this point).
    EvaluatedTrials evaluated =
        externalEnergiesOfTrials(context, component,
                                 generateRetraceTrials(random, forceField, beta, component, chain_atoms, step,
                                                       forceField.numberOfTrialDirections),
                                 std::nullopt);

    // The old configuration is an accepted state of the simulation: it can not overlap, so it survives
    // the overlap filter as trial direction 0. An overlap means the system is inconsistent and no weight
    // is defined for it; fail loudly rather than weigh a bogus W_old.
    if (!evaluated.firstSurvived)
    {
      std::string beads{};
      for (std::size_t bead : step.nextBeads) beads += std::format(" {}", bead);
      throw std::runtime_error(std::format(
          "CBMC: the existing configuration of component '{}' overlaps at growth step {} (bead(s){}); the "
          "retrace of an overlapping molecule has no defined weight. The simulation state is inconsistent "
          "(overlapping molecules in the initial/restart configuration, a molecule inside a blocked pocket, or a "
          "force field or scaling changed after placement).\n",
          component.name, seg, beads));
    }

    const std::vector<ChainTrialTorsion> &trials = evaluated.trials;
    const StepWeight weight =
        stepWeight(random, beta, forceField.numberOfTrialDirections, step, chain_atoms, trials, true);

    chain_external_energies += trials.front().energy;
    chain_log_rosenbluth_weight += weight.logWeight;
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  return ChainRetraceData(chain_external_energies + internal_energies, chain_log_rosenbluth_weight, 0.0);
}
