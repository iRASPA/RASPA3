module;

module cbmc_multiple_first_bead;

import std;

import cbmc_util;
import atom;
import randomnumbers;
import cbmc_first_bead_data;
import cbmc_interactions;
import running_energy;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber& random, const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions(forceField.numberOfFirstBeadPositions, atom);

  // create trial positions randomly in the simulation box
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, 
      externalFieldInterpolationGrid, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  // if all positions over lap return failure
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltzmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  std::size_t selected = selectTrialPosition(random, logBoltzmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [&](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second,
                       RosenbluthWeight / double(forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceMultipleFirstBeadSwapDeletion(
    RandomNumber& random, const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom) noexcept
{
  std::vector<Atom> trialPositions(forceField.numberOfFirstBeadPositions, atom);

  // set the trial positions of the first bead randomly in the simulation box for the 1..N_trial atomns, but leave the
  // first as the old
  std::for_each(trialPositions.begin() + 1, trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, externalFieldInterpolationGrid,
      framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  std::vector<double> logBoltzmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltzmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  return FirstBeadData(atom, externalEnergies[0].second,
                       RosenbluthWeight / double(forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMultipleFirstBeadReinsertion(
    RandomNumber& random, const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<Atom> trialPositions(forceField.numberOfFirstBeadPositions, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, 
      externalFieldInterpolationGrid, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions, skipBackgroundMolecule);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltzmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  std::size_t selected = CBMC::selectTrialPosition(random, logBoltzmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [&](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  // r=w(n)-exp(-beta U[h_n]) Eq.16 from Esselink et al.
  double storedR = RosenbluthWeight - std::exp(logBoltzmannFactors[selected]);

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second,
                       RosenbluthWeight / double(forceField.numberOfFirstBeadPositions), storedR);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::retraceMultipleFirstBeadReinsertion(
    [[maybe_unused]] RandomNumber& random, const Component& component, bool hasExternalField,
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom, double storedR,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, externalFieldInterpolationGrid, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions, skipBackgroundMolecule);
  if (externalEnergies.empty())
  {
    return std::nullopt;
  }

  std::vector<double> logBoltzmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltzmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [](const double& acc, const double& logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  // w(o)=exp(-beta u(o))+r  Eq. 18 from Esselink et al.
  return FirstBeadData(atom, externalEnergies[0].second,
                       (RosenbluthWeight + storedR) / double(forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMultipleFirstBeadPartialInsertion(
    const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom, std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, externalFieldInterpolationGrid,
      framework, frameworkAtoms, moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions,
      skipBackgroundMolecule);

  if (externalEnergies.empty()) return std::nullopt;

  double logBoltzmannFactor = -beta * externalEnergies[0].second.potentialEnergy();
  double RosenbluthWeight = std::exp(logBoltzmannFactor);

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[0].first, externalEnergies[0].second, RosenbluthWeight, 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceMultipleFirstBeadPartialDeletion(
    const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, externalFieldInterpolationGrid,
      framework, frameworkAtoms, moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  double logBoltzmannFactor = -beta * externalEnergies[0].second.potentialEnergy();
  double RosenbluthWeight = std::exp(logBoltzmannFactor);

  return FirstBeadData(atom, externalEnergies[0].second, RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growFirstBeadAtFixedPosition(
    const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid>& externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, externalFieldInterpolationGrid,
      framework, frameworkAtoms, moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  return FirstBeadData(externalEnergies[0].first, externalEnergies[0].second, 1.0, 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceFirstBeadAtFixedPosition(
    const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<InterpolationEnergyGrid>& externalFieldInterpolationGrid,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, [[maybe_unused]] double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom) noexcept
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, externalFieldInterpolationGrid,
      framework, frameworkAtoms, moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  return FirstBeadData(atom, externalEnergies[0].second, 1.0, 0.0);
}
