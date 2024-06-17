module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module cbmc_multiple_first_bead;

#ifndef USE_LEGACY_HEADERS
import <optional>;
import <cmath>;
import <vector>;
import <algorithm>;
import <numeric>;
import <iterator>;
import <span>;
import <tuple>;
#endif

import cbmc_util;
import atom;
import randomnumbers;
import cbmc_first_bead_data;
import cbmc_interactions;
import running_energy;
import forcefield;
import simulationbox;

[[nodiscard]] std::optional<FirstBeadData> CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber& random, bool hasExternalField, const ForceField& forceField, const SimulationBox& simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOff,
    double cutOffCoulomb, const Atom& atom, size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPositions(numberOfTrialDirections, atom);

  // create trial positions randomly in the simulation box
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies =
      computeExternalNonOverlappingEnergies(hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms,
                                            cutOff, cutOffCoulomb, trialPositions);

  // if all positions over lap return failure
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  size_t selected = selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                        [&](const double& acc, const double& logBoltmannFactor)
                                        { return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second,
                       RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

[[nodiscard]] FirstBeadData CBMC::retraceRigidMultipleFirstBeadSwapDeletion(
    RandomNumber& random, bool hasExternalField, const ForceField& forcefield, const SimulationBox& simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOff,
    double cutOffCoulomb, const Atom atom, double scaling, size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPositions(numberOfTrialDirections, atom);
  for (Atom& trialPosition : trialPositions)
  {
    trialPosition.setScaling(scaling);
  }

  // set the trial positions of the first bead randomly in the simulation box for the 1..N_trial atomns, but leave the
  // first as the old
  std::for_each(trialPositions.begin() + 1, trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies =
      computeExternalNonOverlappingEnergies(hasExternalField, forcefield, simulationBox, frameworkAtoms, moleculeAtoms,
                                            cutOff, cutOffCoulomb, trialPositions);

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight =
      std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                  [](const double& acc, const double& logBoltmannFactor) { return acc + std::exp(logBoltmannFactor); });

  return FirstBeadData(atom, externalEnergies[0].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> CBMC::growRigidMultipleFirstBeadReinsertion(
    RandomNumber& random, bool hasExternalField, const ForceField& forceField, const SimulationBox& simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOff,
    double cutOffCoulomb, const Atom& atom, size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPositions(numberOfTrialDirections, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
                [&](Atom& a) { a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies =
      computeExternalNonOverlappingEnergies(hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms,
                                            cutOff, cutOffCoulomb, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                        [&](const double& acc, const double& logBoltmannFactor)
                                        { return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  // r=w(n)-exp(-beta U[h_n]) Eq.16 from Esselink et al.
  double storedR = RosenbluthWeight - std::exp(logBoltmannFactors[selected]);

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second,
                       RosenbluthWeight / double(numberOfTrialDirections), storedR);
}

[[nodiscard]] FirstBeadData CBMC::retraceRigidMultipleFirstBeadReinsertion(
    [[maybe_unused]] RandomNumber& random, bool hasExternalField, const ForceField& forceField,
    const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
    double beta, double cutOff, double cutOffCoulomb, const Atom& atom, double storedR, size_t numberOfTrialDirections)
{
  std::vector<Atom> trialPositions({atom});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies =
      computeExternalNonOverlappingEnergies(hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms,
                                            cutOff, cutOffCoulomb, trialPositions);
  if (externalEnergies.empty())
  {
    throw std::runtime_error(
        "[retraceMultipleFirstBeadReinsertion]: all overlap, \
                              including existing configuration\n");
  }

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<Atom, RunningEnergy>& v) { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight =
      std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                  [](const double& acc, const double& logBoltmannFactor) { return acc + std::exp(logBoltmannFactor); });

  // w(o)=exp(-beta u(o))+r  Eq. 18 from Esselink et al.
  return FirstBeadData(atom, externalEnergies[0].second, (RosenbluthWeight + storedR) / double(numberOfTrialDirections),
                       0.0);
}
