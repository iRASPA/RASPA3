module;

module system;

import randomnumbers;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import energy_status;
import cbmc;
import cbmc_growing_status;
import forcefield;
import energy_factor;
import running_energy;

import <vector>;
import <tuple>;
import <optional>;
import <span>;

import <iostream>;
import <algorithm>;
import <numeric>;

[[nodiscard]] std::optional<ChainData> System::growMoleculeReinsertion(double cutOff, double cutOffCoulomb, size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule, std::span<Atom> molecule) const noexcept
{
  std::vector<Atom> atoms = components[selectedComponent].copiedAtoms(molecule);
  size_t startingBead = components[selectedComponent].startingBead;

  std::optional<FirstBeadData> const firstBeadData = growMultipleFirstBeadReinsertion(cutOff, cutOffCoulomb, atoms[startingBead]);

  if (!firstBeadData) return std::nullopt;

  std::for_each(atoms.begin(), atoms.end(), [&](Atom& atom) {atom.position += firstBeadData->atom.position; });

  if(molecule.size() == 1)
  {
    return ChainData({firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, firstBeadData->storedR);
  }

  std::optional<ChainData> rigidRotationData = growChain(cutOff, cutOffCoulomb, startingBead, atoms);
  if (!rigidRotationData) return std::nullopt;
  return ChainData(rigidRotationData->atom, firstBeadData->energies + rigidRotationData->energies, firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, firstBeadData->storedR);
}

[[nodiscard]] ChainData System::retraceMoleculeReinsertion(double cutOff, double cutOffCoulomb, size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule, std::span<Atom> molecule, double storedR) const noexcept
{
  size_t startingBead = components[selectedComponent].startingBead;

  const FirstBeadData firstBeadData = retraceMultipleFirstBeadReinsertion(cutOff, cutOffCoulomb, molecule[startingBead], storedR);

  if(molecule.size() == 1)
  {
    return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainData rigidRotationData = retraceChainReinsertion(cutOff, cutOffCoulomb, startingBead, molecule);
  return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), firstBeadData.energies + rigidRotationData.energies, firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}



[[nodiscard]] std::optional<FirstBeadData> System::growMultipleFirstBeadReinsertion(double cutOff, double cutOffCoulomb, const Atom& atom) const noexcept
{
  std::vector<Atom> trialPositions(numberOfTrialDirections, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
      [&](Atom& a) {a.position = simulationBox.randomPosition();});

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(cutOff, cutOffCoulomb, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(),
      std::back_inserter(logBoltmannFactors), [this](const std::pair<Atom, RunningEnergy>& v) {return -Beta * v.second.total(); });

  size_t selected = selectTrialPosition(logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [&](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < minimumRosenbluthFactor) return std::nullopt;

  // r=w(n)-exp(-beta U[h_n]) Eq.16 from Esselink et al.
  double storedR = RosenbluthWeight - std::exp(logBoltmannFactors[selected]);

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second, RosenbluthWeight / double(numberOfTrialDirections), storedR);
}

[[nodiscard]] FirstBeadData System::retraceMultipleFirstBeadReinsertion(double cutOff, double cutOffCoulomb, const Atom& atom, double storedR) const noexcept
{
  std::vector<Atom> trialPositions({ atom });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(cutOff, cutOffCoulomb, trialPositions);

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
      [this](const std::pair<Atom, RunningEnergy>& v) {return -Beta * v.second.total(); });

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  // w(o)=exp(-beta u(o))+r  Eq. 18 from Esselink et al.
  return FirstBeadData(atom, externalEnergies[0].second, (RosenbluthWeight + storedR) / double(numberOfTrialDirections), 0.0);
}

[[nodiscard]] ChainData System::retraceChainReinsertion(double cutOff, double cutOffCoulomb, size_t startingBead, std::span<Atom> molecule) const noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule.begin(), molecule.end());
  std::vector<std::vector<Atom>> trialPositions = { trialPosition };

  for (size_t i = 1; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(rotateRandomlyAround(trialPosition, startingBead));
  };

  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(cutOff, cutOffCoulomb, trialPositions, std::make_signed_t<std::size_t>(startingBead));

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies),
      std::back_inserter(logBoltmannFactors), [this](const std::pair<std::vector<Atom>, RunningEnergy>& v) {return -Beta * v.second.total(); });

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  return ChainData(trialPositions[0], externalEnergies[0].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

