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

/// Sum of the Boltzmann factors of a set of trial positions (overlapping trials contribute zero).
double sumOfBoltzmannFactors(const GrowContext &context, const std::vector<FirstBeadTrial> &trials)
{
  return std::accumulate(trials.begin(), trials.end(), 0.0,
                         [&](double acc, const FirstBeadTrial &trial)
                         { return acc + std::exp(-context.beta * trial.energy.potentialEnergy()); });
}

/// Rosenbluth selection among random trial positions: the selected trial, the sum of the Boltzmann
/// factors, and the Boltzmann factor of the selected trial. std::nullopt when every trial overlaps
/// or the weight falls below 'minimumRosenbluthFactor'.
struct MultipleFirstBeadSelection
{
  FirstBeadTrial selected;
  double rosenbluthSum;
  double selectedBoltzmannFactor;
};

std::optional<MultipleFirstBeadSelection> selectAmongRandomPositions(RandomNumber &random, const GrowContext &context,
                                                                     const Component &component, const Atom &atom)
{
  std::vector<Atom> trialPositions(context.settings.numberOfFirstBeadPositions, atom);
  for (Atom &trial : trialPositions) trial.position = context.simulationBox.randomPosition(random);

  const std::vector<FirstBeadTrial> trials = computeExternalNonOverlappingEnergies(context, component, trialPositions);
  if (trials.empty()) return std::nullopt;

  std::vector<double> logBoltzmannFactors(trials.size());
  std::transform(trials.begin(), trials.end(), logBoltzmannFactors.begin(),
                 [&](const FirstBeadTrial &trial) { return -context.beta * trial.energy.potentialEnergy(); });

  const std::size_t selected = selectTrialPosition(random, logBoltzmannFactors);

  const double rosenbluthSum = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                               [](double acc, double logFactor) { return acc + std::exp(logFactor); });
  if (rosenbluthSum < context.settings.minimumRosenbluthFactor) return std::nullopt;

  return MultipleFirstBeadSelection{trials[selected], rosenbluthSum, std::exp(logBoltzmannFactors[selected])};
}

// --- Multiple first bead: 'numberOfFirstBeadPositions' random positions, weight = sum / number ------

std::optional<FirstBeadData> growMultipleFirstBead(RandomNumber &random, const GrowContext &context,
                                                   const Component &component, const Atom &atom)
{
  const std::optional<MultipleFirstBeadSelection> selection = selectAmongRandomPositions(random, context, component, atom);
  if (!selection) return std::nullopt;

  return FirstBeadData(selection->selected.position, selection->selected.energy,
                       std::log(selection->rosenbluthSum / double(context.settings.numberOfFirstBeadPositions)), 0.0);
}

FirstBeadData retraceMultipleFirstBead(RandomNumber &random, const GrowContext &context, const Component &component,
                                       const Atom &atom)
{
  // The existing bead is trial 0; the remaining positions are drawn at random.
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);

  std::vector<Atom> trialPositions(context.settings.numberOfFirstBeadPositions - 1, atom);
  for (Atom &trial : trialPositions) trial.position = context.simulationBox.randomPosition(random);

  const std::vector<FirstBeadTrial> trials = computeExternalNonOverlappingEnergies(context, component, trialPositions);

  const double rosenbluthSum =
      std::exp(-context.beta * oldEnergy.potentialEnergy()) + sumOfBoltzmannFactors(context, trials);

  return FirstBeadData(atom, oldEnergy, std::log(rosenbluthSum / double(context.settings.numberOfFirstBeadPositions)),
                       0.0);
}

// --- Multiple-first-bead reinsertion (Esselink et al.): retains the partial weight r ----------------

std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(RandomNumber &random, const GrowContext &context,
                                                              const Component &component, const Atom &atom)
{
  const std::optional<MultipleFirstBeadSelection> selection = selectAmongRandomPositions(random, context, component, atom);
  if (!selection) return std::nullopt;

  // r = w(n) - exp(-beta U[h_n]), Eq. 16 of Esselink et al. (kept linear, see FirstBeadData::storedR)
  const double storedR = selection->rosenbluthSum - selection->selectedBoltzmannFactor;

  return FirstBeadData(selection->selected.position, selection->selected.energy,
                       std::log(selection->rosenbluthSum / double(context.settings.numberOfFirstBeadPositions)),
                       storedR);
}

FirstBeadData retraceMultipleFirstBeadReinsertion(const GrowContext &context, const Component &component,
                                                  const Atom &atom, double storedR)
{
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);

  // w(o) = exp(-beta u(o)) + r, Eq. 18 of Esselink et al.
  const double rosenbluthSum = std::exp(-context.beta * oldEnergy.potentialEnergy()) + storedR;

  return FirstBeadData(atom, oldEnergy, std::log(rosenbluthSum / double(context.settings.numberOfFirstBeadPositions)),
                       0.0);
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
  if (std::exp(logBoltzmannFactor) < context.settings.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(atom, energy.value(), logBoltzmannFactor, 0.0);
}

FirstBeadData retracePinnedFirstBead(const GrowContext &context, const Component &component, const Atom &atom)
{
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);
  return FirstBeadData(atom, oldEnergy, -context.beta * oldEnergy.potentialEnergy(), 0.0);
}

// --- Fixed first bead: a single trial at the given position with weight one -------------------------

std::optional<FirstBeadData> growFixedFirstBead(const GrowContext &context, const Component &component,
                                                const Atom &atom)
{
  const std::optional<RunningEnergy> energy =
      computeExternalNonOverlappingEnergy(context, component, std::span<const Atom>(&atom, 1));
  if (!energy) return std::nullopt;

  return FirstBeadData(atom, energy.value(), 0.0, 0.0);
}

FirstBeadData retraceFixedFirstBead(const GrowContext &context, const Component &component, const Atom &atom)
{
  const RunningEnergy oldEnergy = existingFirstBeadEnergy(context, component, atom);
  return FirstBeadData(atom, oldEnergy, 0.0, 0.0);
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

CBMC::FirstBeadData CBMC::retraceFirstBead(RandomNumber &random, const GrowContext &context,
                                           const Component &component, const Atom &firstBead,
                                           FirstBeadScheme scheme, double storedR)
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return retraceMultipleFirstBead(random, context, component, firstBead);
    case FirstBeadScheme::Reinsertion:
      return retraceMultipleFirstBeadReinsertion(context, component, firstBead, storedR);
    case FirstBeadScheme::Pinned:
      return retracePinnedFirstBead(context, component, firstBead);
    case FirstBeadScheme::Fixed:
      return retraceFixedFirstBead(context, component, firstBead);
    case FirstBeadScheme::AlreadyPlaced:
      break;
  }
  throw std::invalid_argument("CBMC: no first-bead stage for FirstBeadScheme::AlreadyPlaced");
}
