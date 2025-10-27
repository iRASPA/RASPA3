module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <optional>
#include <print>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>
#endif

module cbmc_multiple_first_bead;

#ifdef USE_STD_IMPORT
import std;
#endif

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
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions(forceField.numberOfFirstBeadPositions, atom);

  // create trial positions randomly in the simulation box
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  // if all positions over lap return failure
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  std::size_t selected = selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [&](const double& acc, const double& logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second,
                       RosenbluthWeight / double(forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceMultipleFirstBeadSwapDeletion(
    RandomNumber& random, const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom) noexcept
{
  std::vector<Atom> trialPositions(forceField.numberOfFirstBeadPositions, atom);

  // set the trial positions of the first bead randomly in the simulation box for the 1..N_trial atomns, but leave the
  // first as the old
  std::for_each(trialPositions.begin() + 1, trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [](const double& acc, const double& logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  return FirstBeadData(atom, externalEnergies[0].second,
                       RosenbluthWeight / double(forceField.numberOfFirstBeadPositions), 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMultipleFirstBeadReinsertion(
    RandomNumber& random, const Component& component, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom) noexcept
{
  std::vector<Atom> trialPositions(forceField.numberOfFirstBeadPositions, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  std::size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [&](const double& acc, const double& logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  // r=w(n)-exp(-beta U[h_n]) Eq.16 from Esselink et al.
  double storedR = RosenbluthWeight - std::exp(logBoltmannFactors[selected]);

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second,
                       RosenbluthWeight / double(forceField.numberOfFirstBeadPositions), storedR);
}

[[nodiscard]] FirstBeadData CBMC::retraceMultipleFirstBeadReinsertion(
    [[maybe_unused]] RandomNumber& random, const Component& component, bool hasExternalField,
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom& atom, double storedR)
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(
      component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions);
  if (externalEnergies.empty())
  {
    throw std::runtime_error(
        "[retraceMultipleFirstBeadReinsertion]: all overlap, \
                              including existing configuration\n");
  }

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [](const double& acc, const double& logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  // w(o)=exp(-beta u(o))+r  Eq. 18 from Esselink et al.
  return FirstBeadData(atom, externalEnergies[0].second,
                       (RosenbluthWeight + storedR) / double(forceField.numberOfFirstBeadPositions), 0.0);
}
