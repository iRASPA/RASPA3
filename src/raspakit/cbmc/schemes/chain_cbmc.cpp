module;

module cbmc_chain_cbmc;

import std;

import randomnumbers;
import component;
import molecule;
import atom;
import running_energy;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_results;
import cbmc_grow_context;
import cbmc_grow_step;
import cbmc_operators;
import cbmc_external_energy;
import cbmc_chain_common;

namespace
{
/// A non-overlapping trial direction of a step: its positions, its energies, and the torsion Rosenbluth
/// weight accumulated while generating it.
struct EvaluatedTrial
{
  std::vector<Atom> positions;
  CBMC::StepTrialEnergy energy;
  double torsionWeight;
};

// The trial directions of a step: the operator engine's trials, filtered to the non-overlapping ones
// and paired with their energies. 'firstSurvived' records explicitly whether trial direction 0 passed
// the filter (the retrace needs to know that its old configuration, which it pins as trial 0, is still
// among the survivors).
struct EvaluatedTrials
{
  std::vector<EvaluatedTrial> trials{};
  bool firstSurvived{false};
};

EvaluatedTrials evaluateTrials(const CBMC::GrowContext &context, const Component &component,
                               const CBMC::GrowStep &step, std::vector<Atom> &chainAtoms,
                               std::vector<CBMC::StepTrial> stepTrials)
{
  EvaluatedTrials evaluated{};
  evaluated.trials.reserve(stepTrials.size());
  for (std::size_t i = 0; i != stepTrials.size(); ++i)
  {
    std::optional<CBMC::StepTrialEnergy> energy =
        CBMC::evaluateStepTrial(context, component, step, chainAtoms, stepTrials[i].positions);
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
// with u_j the external plus the intramolecular van der Waals and Coulomb energy of trial j and k the
// number of trial directions. On the grow the selected trial is drawn with probability
// exp(-beta u_j) / sum; on the retrace it is trial 0, the old configuration. The selected trial's
// positions are left in 'chainAtoms'.
//
// The Boltzmann sum is evaluated as log-sum-exp so the log stays exact even where the raw factor or
// the running product of a long chain would underflow.
StepWeight stepWeight(RandomNumber &random, double beta, std::size_t numberOfTrialDirections,
                      const CBMC::GrowStep &step, std::vector<Atom> &chainAtoms,
                      const std::vector<EvaluatedTrial> &trials, bool retrace)
{
  std::vector<double> logBoltzmannFactors{};
  logBoltzmannFactors.reserve(trials.size());
  for (const EvaluatedTrial &trial : trials)
  {
    logBoltzmannFactors.push_back(-beta * trial.energy.potentialEnergy());
  }

  const std::size_t selected = retrace ? 0 : CBMC::selectTrialPosition(random, logBoltzmannFactors);
  const EvaluatedTrial &selectedTrial = trials[selected];
  CBMC::placeStepBeads(chainAtoms, step, selectedTrial.positions);

  const double maxLogBoltzmannFactor = *std::max_element(logBoltzmannFactors.begin(), logBoltzmannFactors.end());
  const double logRosenbluthSum =
      maxLogBoltzmannFactor +
      std::log(std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                               [&](const double &acc, const double &logBoltzmannFactor)
                               { return acc + std::exp(logBoltzmannFactor - maxLogBoltzmannFactor); }));

  return {selected, std::log(selectedTrial.torsionWeight) + logRosenbluthSum -
                        std::log(static_cast<double>(numberOfTrialDirections)) +
                        CBMC::unsampledInternalLogFactor(beta, step, chainAtoms)};
}
}  // namespace

[[nodiscard]] std::optional<CBMC::GrowResult> CBMC::growChainCBMC(RandomNumber &random, const GrowContext &context,
                                                                  const Component &component,
                                                                  std::span<const Atom> molecule_atoms,
                                                                  const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const GrowthSettings &settings = context.settings;
  const double beta = context.beta;

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
  ChainAccumulator chain{};

  // Deterministic growth plan over the fragment graph (flexible beads, hinged rigid bodies, and
  // ring-closure of cyclic clusters), shared with the retrace so grow and retrace are reversible.
  const std::vector<GrowStep> &plan = component.growthPlan(beadsAlreadyPlaced);

  for (const GrowStep &step : plan)
  {
    // All trial directions of this step (the operator engine handles the seed / attach / ring-closure
    // cases, the rigid-body tilt, and the coupled-decoupled torsion selection).
    std::vector<EvaluatedTrial> trials =
        evaluateTrials(context, component, step, chain_atoms,
                       generateGrowTrials(random, settings, beta, component, chain_atoms, step,
                                          settings.numberOfTrialDirections))
            .trials;
    if (trials.empty()) return std::nullopt;

    const StepWeight weight =
        stepWeight(random, beta, settings.numberOfTrialDirections, step, chain_atoms, trials, false);
    if (!chain.addGrownStep(weight.logWeight, trials[weight.selected].energy.external,
                            settings.minimumRosenbluthFactor))
    {
      return std::nullopt;
    }
  }

  return finishGrownChain(component, std::move(chain_atoms), chain);
}

[[nodiscard]] CBMC::RetraceResult CBMC::retraceChainCBMC(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<const Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const GrowthSettings &settings = context.settings;
  const double beta = context.beta;

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
  ChainAccumulator chain{};

  // Same deterministic growth plan as the insertion so grow and retrace are exactly reversible.
  const std::vector<GrowStep> &plan = component.growthPlan(beadsAlreadyPlaced);

  for (std::size_t seg = 0; seg != plan.size(); ++seg)
  {
    const GrowStep &step = plan[seg];

    // The old positions of this step's beads are trial direction 0 of the retrace (read from
    // 'chain_atoms', which holds the old configuration of every bead at this point).
    EvaluatedTrials evaluated =
        evaluateTrials(context, component, step, chain_atoms,
                       generateRetraceTrials(random, settings, beta, component, chain_atoms, step,
                                             settings.numberOfTrialDirections));

    // The old configuration must survive the overlap filter as trial direction 0.
    if (!evaluated.firstSurvived) throwExistingConfigurationOverlaps("CBMC", component, seg, step);

    const std::vector<EvaluatedTrial> &trials = evaluated.trials;
    const StepWeight weight =
        stepWeight(random, beta, settings.numberOfTrialDirections, step, chain_atoms, trials, true);
    chain.addRetracedStep(weight.logWeight, trials.front().energy.external);
  }

  return finishRetracedChain(component, molecule_atoms, chain);
}
