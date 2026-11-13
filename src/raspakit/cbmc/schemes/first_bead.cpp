module;

module cbmc_first_bead;

import std;

import cbmc_util;
import atom;
import randomnumbers;
import cbmc_results;
import cbmc_external_energy;
import cbmc_grow_context;
import running_energy;
import component;

namespace
{
using namespace CBMC;

/// External energy of the existing first bead; throws when it overlaps (an accepted configuration can
/// not overlap, so a weight is undefined: see the error contract in the 'cbmc' module).
RunningEnergy existingFirstBeadEnergy(const GrowContext &context, const Component &component, const Atom &atom)
{
  const std::optional<RunningEnergy> energy =
      computeExternalNonOverlappingEnergy(context, component, std::span<const Atom>(&atom, 1));
  if (!energy.has_value())
  {
    throw std::runtime_error(
        std::format("CBMC retrace: the first bead of the existing molecule (component '{}', molecule {}) overlaps "
                    "with its environment; the configuration is inconsistent and has no Rosenbluth weight",
                    component.name, atom.moleculeId));
  }
  return energy.value();
}

/// log(sum_i exp(x_i)), shifted by the largest term. The Rosenbluth sums here are formed in log space: the
/// Boltzmann factor of a single bead at a strongly repulsive site (-beta U below about -745, well inside the
/// overlap criterion) is zero in double precision, and a sum of such terms, or of an existing bead alone, gave
/// log(0) = -infinity, an acceptance probability of 0/0 or x/0 in the move that used it. (The same helper
/// exists in the multi-trial module, which sits above this one in the module graph.) Precondition: not empty.
double logSumExp(std::span<const double> values)
{
  const double maximum = *std::ranges::max_element(values);
  double sum = 0.0;
  for (double value : values) sum += std::exp(value - maximum);
  return maximum + std::log(sum);
}

/// The log-Boltzmann factors -beta U of a set of trial positions (overlapping trials are not among them).
std::vector<double> logBoltzmannFactorsOf(const GrowContext& context, const std::vector<FirstBeadTrial>& trials)
{
  std::vector<double> logFactors(trials.size());
  std::transform(trials.begin(), trials.end(), logFactors.begin(),
                 [&](const FirstBeadTrial& trial) { return -context.beta * trial.energy.potentialEnergy(); });
  return logFactors;
}

/// Rosenbluth selection among random trial positions: the selected trial, the logarithm of the sum of the
/// Boltzmann factors, and the log-Boltzmann factor of the selected trial. std::nullopt when every trial
/// overlaps or the weight falls below 'minimumRosenbluthFactor'.
struct MultipleFirstBeadSelection
{
  FirstBeadTrial selected;
  double logRosenbluthSum;
  /// Log-sum of the Boltzmann factors of the trials other than the selected one; empty when there is none.
  std::optional<double> logRemainingSum;
};

std::optional<MultipleFirstBeadSelection> selectAmongRandomPositions(RandomNumber &random, const GrowContext &context,
                                                                     const Component &component, const Atom &atom)
{
  std::vector<Atom> trialPositions(context.settings.numberOfFirstBeadPositions, atom);
  for (Atom &trial : trialPositions) trial.position = context.simulationBox.randomPosition(random);

  const std::vector<FirstBeadTrial> trials = computeExternalNonOverlappingEnergies(context, component, trialPositions);
  if (trials.empty()) return std::nullopt;

  const std::vector<double> logBoltzmannFactors = logBoltzmannFactorsOf(context, trials);

  const auto [selected, logRosenbluthSum] = selectTrialPosition(random, logBoltzmannFactors);

  // compared as logarithms: the exponential of the log-sum of a strongly attractive site overflows
  if (context.settings.minimumRosenbluthFactor > 0.0 &&
      logRosenbluthSum < std::log(context.settings.minimumRosenbluthFactor))
  {
    return std::nullopt;
  }

  // summed in log space over the other trials themselves: the sum minus the selected factor cancelled to zero when
  // the selected trial dominated, and the sum itself overflowed for a strongly attractive site
  std::vector<double> remaining{};
  remaining.reserve(logBoltzmannFactors.size());
  for (std::size_t i = 0; i != logBoltzmannFactors.size(); ++i)
  {
    if (i != selected) remaining.push_back(logBoltzmannFactors[i]);
  }
  std::optional<double> logRemainingSum{};
  if (!remaining.empty()) logRemainingSum = logSumExp(remaining);

  return MultipleFirstBeadSelection{trials[selected], logRosenbluthSum, logRemainingSum};
}

// --- Multiple first bead: 'numberOfFirstBeadPositions' random positions, weight = sum / number ------

std::optional<FirstBeadData> growMultipleFirstBead(RandomNumber &random, const GrowContext &context,
                                                   const Component &component, const Atom &atom)
{
  const std::optional<MultipleFirstBeadSelection> selection = selectAmongRandomPositions(random, context, component, atom);
  if (!selection) return std::nullopt;

  return FirstBeadData(selection->selected.position, selection->selected.energy,
                       selection->logRosenbluthSum - std::log(double(context.settings.numberOfFirstBeadPositions)),
                       std::nullopt);
}

FirstBeadData retraceMultipleFirstBead(RandomNumber &random, const GrowContext &context, const Component &component,
                                       const Atom &atom)
{
  // The existing bead is trial 0; the remaining positions are drawn at random.
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);

  std::vector<Atom> trialPositions(context.settings.numberOfFirstBeadPositions - 1, atom);
  for (Atom &trial : trialPositions) trial.position = context.simulationBox.randomPosition(random);

  const std::vector<FirstBeadTrial> trials = computeExternalNonOverlappingEnergies(context, component, trialPositions);

  // the existing bead is always among the terms: a sum of Boltzmann factors that underflows (a strongly
  // repulsive existing bead, every random trial overlapping) is kept in log space
  std::vector<double> logBoltzmannFactors = logBoltzmannFactorsOf(context, trials);
  logBoltzmannFactors.push_back(-context.beta * oldEnergy.potentialEnergy());

  return FirstBeadData(atom, oldEnergy,
                       logSumExp(logBoltzmannFactors) - std::log(double(context.settings.numberOfFirstBeadPositions)),
                       std::nullopt);
}

// --- Multiple-first-bead reinsertion (Esselink et al.): retains the partial weight r ----------------

std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(RandomNumber &random, const GrowContext &context,
                                                              const Component &component, const Atom &atom)
{
  const std::optional<MultipleFirstBeadSelection> selection = selectAmongRandomPositions(random, context, component, atom);
  if (!selection) return std::nullopt;

  // r = w(n) - exp(-beta U[h_n]), Eq. 16 of Esselink et al.: the sum over the trials other than the selected one,
  // kept as its logarithm (see FirstBeadData::logStoredR), empty when there are none
  return FirstBeadData(selection->selected.position, selection->selected.energy,
                       selection->logRosenbluthSum - std::log(double(context.settings.numberOfFirstBeadPositions)),
                       selection->logRemainingSum);
}

FirstBeadData retraceMultipleFirstBeadReinsertion(const GrowContext& context, const Component& component,
                                                  const Atom& atom, std::optional<double> logStoredR)
{
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);

  // w(o) = exp(-beta u(o)) + r, Eq. 18 of Esselink et al., in log space: the Boltzmann factor of a strongly
  // repulsive existing bead is zero in double precision, and with r = 0 (no other trial) the sum was log(0)
  std::vector<double> logTerms{-context.beta * oldEnergy.potentialEnergy()};
  if (logStoredR.has_value()) logTerms.push_back(logStoredR.value());

  return FirstBeadData(atom, oldEnergy,
                       logSumExp(logTerms) - std::log(double(context.settings.numberOfFirstBeadPositions)));
}

// --- Pinned first bead: a single trial at the given position, its Boltzmann factor as weight --------

std::optional<FirstBeadData> growPinnedFirstBead(const GrowContext &context, const Component &component,
                                                 const Atom &atom)
{
  const std::optional<RunningEnergy> energy =
      computeExternalNonOverlappingEnergy(context, component, std::span<const Atom>(&atom, 1));
  if (!energy) return std::nullopt;

  // A single trial: the weight is the Boltzmann factor itself, already available as its logarithm.
  const double logBoltzmannFactor = -context.beta * energy->potentialEnergy();
  if (context.settings.minimumRosenbluthFactor > 0.0 &&
      logBoltzmannFactor < std::log(context.settings.minimumRosenbluthFactor))
  {
    return std::nullopt;
  }

  return FirstBeadData(atom, energy.value(), logBoltzmannFactor);
}

FirstBeadData retracePinnedFirstBead(const GrowContext &context, const Component &component, const Atom &atom)
{
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);
  return FirstBeadData(atom, oldEnergy, -context.beta * oldEnergy.potentialEnergy());
}

// --- Fixed first bead: a single trial at the given position with weight one -------------------------

std::optional<FirstBeadData> growFixedFirstBead(const GrowContext &context, const Component &component,
                                                const Atom &atom)
{
  const std::optional<RunningEnergy> energy =
      computeExternalNonOverlappingEnergy(context, component, std::span<const Atom>(&atom, 1));
  if (!energy) return std::nullopt;

  return FirstBeadData(atom, energy.value(), 0.0);
}

FirstBeadData retraceFixedFirstBead(const GrowContext &context, const Component &component, const Atom &atom)
{
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);
  return FirstBeadData(atom, oldEnergy, 0.0);
}
}  // namespace

std::optional<CBMC::FirstBeadData> CBMC::growFirstBead(RandomNumber &random, const GrowContext &context,
                                                       const Component &component, const Atom &firstBead,
                                                       FirstBeadScheme scheme)
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return growMultipleFirstBead(random, context, component, firstBead);
    case FirstBeadScheme::Reinsertion:
      return growMultipleFirstBeadReinsertion(random, context, component, firstBead);
    case FirstBeadScheme::Pinned:
      return growPinnedFirstBead(context, component, firstBead);
    case FirstBeadScheme::Fixed:
      return growFixedFirstBead(context, component, firstBead);
    case FirstBeadScheme::AlreadyPlaced:
      break;
  }
  throw std::invalid_argument("CBMC: no first-bead stage for FirstBeadScheme::AlreadyPlaced");
}

CBMC::FirstBeadData CBMC::retraceFirstBead(RandomNumber& random, const GrowContext& context, const Component& component,
                                           const Atom& firstBead, FirstBeadScheme scheme,
                                           std::optional<double> logStoredR)
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return retraceMultipleFirstBead(random, context, component, firstBead);
    case FirstBeadScheme::Reinsertion:
      return retraceMultipleFirstBeadReinsertion(context, component, firstBead, logStoredR);
    case FirstBeadScheme::Pinned:
      return retracePinnedFirstBead(context, component, firstBead);
    case FirstBeadScheme::Fixed:
      return retraceFixedFirstBead(context, component, firstBead);
    case FirstBeadScheme::AlreadyPlaced:
      break;
  }
  throw std::invalid_argument("CBMC: no first-bead stage for FirstBeadScheme::AlreadyPlaced");
}
