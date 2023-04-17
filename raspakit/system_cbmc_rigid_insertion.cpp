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


[[nodiscard]] std::optional<ChainData> System::growMoleculeSwapInsertion(double cutOff, double cutOffCoulomb, size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule, double scaling, [[maybe_unused]] std::vector<Atom> atoms) const noexcept
{
  //std::vector<Atom> atoms = components[selectedComponent].newAtoms(scaling, selectedMolecule);
  for (Atom& atom : atoms)
  {
    atom.setScaling(scaling);
  }
  size_t startingBead = components[selectedComponent].startingBead;

  std::optional<FirstBeadData> const firstBeadData = growMultipleFirstBeadSwapInsertion(cutOff, cutOffCoulomb, atoms[startingBead]);

  if (!firstBeadData) return std::nullopt;

  std::for_each(atoms.begin(), atoms.end(), [&](Atom& atom) {atom.position += firstBeadData->atom.position; });

  if(atoms.size() == 1)
  {
    return ChainData({firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  std::optional<ChainData> const rigidRotationData = growChain(cutOff, cutOffCoulomb, startingBead, atoms);
  
  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->atom, firstBeadData->energies + rigidRotationData->energies, firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<FirstBeadData> System::growMultipleFirstBeadSwapInsertion(double cutOff, double cutOffCoulomb, const Atom& atom) const noexcept
{
  std::vector<Atom> trialPositions(numberOfTrialDirections, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
      [&](Atom& a) {a.position = simulationBox.randomPosition(); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(cutOff, cutOffCoulomb, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(),
      std::back_inserter(logBoltmannFactors), [this](const std::pair<Atom,RunningEnergy>& v) {return -beta * v.second.total(); });


  size_t selected = selectTrialPosition(logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [&](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

[[nodiscard]] std::optional<ChainData> System::growChain(double cutOff, double cutOffCoulomb, size_t startingBead, std::vector<Atom> molecule) const noexcept
{
  std::vector<std::vector<Atom>> trialPositions{};

  for(size_t i = 0; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(rotateRandomlyAround(molecule, startingBead));
  };
  
  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(cutOff, cutOffCoulomb, trialPositions, std::make_signed_t<std::size_t>(startingBead));
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(),
      std::back_inserter(logBoltmannFactors), [&](const std::pair<std::vector<Atom>, RunningEnergy>& v) {return -beta * v.second.total(); });

  size_t selected = selectTrialPosition(logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < minimumRosenbluthFactor) return std::nullopt;

  return ChainData(externalEnergies[selected].first, externalEnergies[selected].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}
