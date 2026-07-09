module;

module cbmc_multiple_first_bead;

import std;

import cbmc_util;
import atom;
import randomnumbers;
import cbmc_first_bead_data;
import cbmc_interactions;
import cbmc_growth_context;
import running_energy;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber& random, const GrowContext& context, const Component& component, const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions(context.forceField.numberOfFirstBeadPositions, atom);

  // create trial positions randomly in the simulation box
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = context.simulationBox.randomPosition(random); });

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions);

  // if all positions over lap return failure
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltzmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                 [&](const FirstBeadTrial& v) { return -context.beta * v.energy.potentialEnergy(); });

  std::size_t selected = selectTrialPosition(random, logBoltzmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [&](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  if (RosenbluthWeight < context.forceField.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[selected].position, externalEnergies[selected].energy,
                       RosenbluthWeight / double(context.forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceMultipleFirstBeadSwapDeletion(RandomNumber& random,
                                                                       const GrowContext& context,
                                                                       const Component& component,
                                                                       const Atom atom) noexcept
{
  std::vector<Atom> trialPositions(context.forceField.numberOfFirstBeadPositions, atom);

  // set the trial positions of the first bead randomly in the simulation box for the 1..N_trial atomns, but leave the
  // first as the old
  std::for_each(trialPositions.begin() + 1, trialPositions.end(),
                [&](Atom& a) { a.position = context.simulationBox.randomPosition(random); });

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions);

  std::vector<double> logBoltzmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltzmannFactors),
                 [&](const FirstBeadTrial& v) { return -context.beta * v.energy.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  return FirstBeadData(atom, externalEnergies[0].energy,
                       RosenbluthWeight / double(context.forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMultipleFirstBeadReinsertion(
    RandomNumber& random, const GrowContext& context, const Component& component, const Atom& atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<Atom> trialPositions(context.forceField.numberOfFirstBeadPositions, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = context.simulationBox.randomPosition(random); });

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions, skipBackgroundMolecule);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltzmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                 [&](const FirstBeadTrial& v) { return -context.beta * v.energy.potentialEnergy(); });

  std::size_t selected = CBMC::selectTrialPosition(random, logBoltzmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [&](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  if (RosenbluthWeight < context.forceField.minimumRosenbluthFactor) return std::nullopt;

  // r=w(n)-exp(-beta U[h_n]) Eq.16 from Esselink et al.
  double storedR = RosenbluthWeight - std::exp(logBoltzmannFactors[selected]);

  return FirstBeadData(externalEnergies[selected].position, externalEnergies[selected].energy,
                       RosenbluthWeight / double(context.forceField.numberOfFirstBeadPositions), storedR);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::retraceMultipleFirstBeadReinsertion(
    [[maybe_unused]] RandomNumber& random, const GrowContext& context, const Component& component, const Atom& atom,
    double storedR, std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions, skipBackgroundMolecule);
  if (externalEnergies.empty())
  {
    return std::nullopt;
  }

  std::vector<double> logBoltzmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltzmannFactors),
                 [&](const FirstBeadTrial& v) { return -context.beta * v.energy.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  // w(o)=exp(-beta u(o))+r  Eq. 18 from Esselink et al.
  return FirstBeadData(atom, externalEnergies[0].energy,
                       (RosenbluthWeight + storedR) / double(context.forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMultipleFirstBeadPartialInsertion(
    const GrowContext& context, const Component& component, const Atom& atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions, skipBackgroundMolecule);

  if (externalEnergies.empty()) return std::nullopt;

  double logBoltzmannFactor = -context.beta * externalEnergies[0].energy.potentialEnergy();
  double RosenbluthWeight = std::exp(logBoltzmannFactor);

  if (RosenbluthWeight < context.forceField.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[0].position, externalEnergies[0].energy, RosenbluthWeight, 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceMultipleFirstBeadPartialDeletion(const GrowContext& context,
                                                                          const Component& component,
                                                                          const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions);

  double logBoltzmannFactor = -context.beta * externalEnergies[0].energy.potentialEnergy();
  double RosenbluthWeight = std::exp(logBoltzmannFactor);

  return FirstBeadData(atom, externalEnergies[0].energy, RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growFirstBeadAtFixedPosition(const GrowContext& context,
                                                                              const Component& component,
                                                                              const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  return FirstBeadData(externalEnergies[0].position, externalEnergies[0].energy, 1.0, 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceFirstBeadAtFixedPosition(const GrowContext& context,
                                                                  const Component& component, const Atom atom) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<FirstBeadTrial> externalEnergies =
      computeExternalNonOverlappingEnergies(context, component, trialPositions);

  return FirstBeadData(atom, externalEnergies[0].energy, 1.0, 0.0);
}
