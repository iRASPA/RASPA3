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

import <vector>;
import <tuple>;
import <optional>;
import <span>;

import <iostream>;
import <algorithm>;
import <numeric>;


// system_cbmc_rigid.cpp

// LogBoltzmannFactors are (-Beta U)
size_t System::selectTrialPosition(RandomNumber &random, std::vector <double> LogBoltzmannFactors) const noexcept
{
  std::vector<double> ShiftedBoltzmannFactors(LogBoltzmannFactors.size());

  // Energies are always bounded from below [-U_max, infinity>
  // Find the lowest energy value, i.e. the largest value of (-Beta U)
  double largest_value = *std::max_element(LogBoltzmannFactors.begin(), LogBoltzmannFactors.end());

  // Standard trick: shift the Boltzmann factors down to avoid numerical problems
  // The largest value of 'ShiftedBoltzmannFactors' will be 1 (which corresponds to the lowest energy).
  double SumShiftedBoltzmannFactors = 0.0;
  for (size_t i = 0; i < LogBoltzmannFactors.size(); ++i)
  {
    ShiftedBoltzmannFactors[i] = exp(LogBoltzmannFactors[i] - largest_value);
    SumShiftedBoltzmannFactors += ShiftedBoltzmannFactors[i];
  }

  // select the Boltzmann factor
  size_t selected = 0;
  double cumw = ShiftedBoltzmannFactors[0];
  double ws = random.uniform() * SumShiftedBoltzmannFactors;
  while (cumw < ws)
    cumw += ShiftedBoltzmannFactors[++selected];

  return selected;
}

[[nodiscard]] std::optional<ChainData> 
System::growMoleculeSwapInsertion(RandomNumber &random, Component::GrowType growType, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                  size_t selectedMolecule, double scaling, std::vector<Atom> atoms) const noexcept
{
  switch(growType)
  {
    default:
    return growRigidMoleculeSwapInsertion(random, cutOff, cutOffCoulomb, selectedComponent, selectedMolecule, scaling, atoms);
  }
}

[[nodiscard]] std::optional<ChainData> 
System::growMoleculeReinsertion(RandomNumber &random, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                size_t selectedMolecule, std::span<Atom> molecule) const noexcept
{
  return growRigidMoleculeReinsertion(random, cutOff, cutOffCoulomb, selectedComponent, selectedMolecule, molecule);
}
[[nodiscard]] ChainData 
System::retraceMoleculeReinsertion(RandomNumber &random, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                   size_t selectedMolecule, std::span<Atom> molecule, double storedR) const noexcept
{
  return retraceRigidMoleculeReinsertion(random, cutOff, cutOffCoulomb, selectedComponent, selectedMolecule, molecule, storedR);
}

[[nodiscard]] ChainData 
System::retraceMoleculeSwapDeletion(RandomNumber &random, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                     size_t selectedMolecule, std::span<Atom> molecule, double scaling, double storedR) const noexcept
{
  return retraceRigidMoleculeSwapDeletion(random, cutOff,cutOffCoulomb, selectedComponent, selectedMolecule, molecule, scaling, storedR);
}


[[nodiscard]] std::optional<FirstBeadData> System::growMoleculeMultipleFirstBeadSwapInsertion(RandomNumber &random, double cutOff, double cutOffCoulomb, const Atom& atom) const noexcept
{
  std::vector<Atom> trialPositions(numberOfTrialDirections, atom);
  std::for_each(trialPositions.begin(), trialPositions.end(),
      [&](Atom& a) {a.position = simulationBox.randomPosition(random); });

  const std::vector<std::pair<Atom, RunningEnergy>> externalEnergies = computeExternalNonOverlappingEnergies(cutOff, cutOffCoulomb, trialPositions);

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(),
      std::back_inserter(logBoltmannFactors), [this](const std::pair<Atom,RunningEnergy>& v) {return -beta * v.second.total(); });


  size_t selected = selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [&](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < minimumRosenbluthFactor) return std::nullopt;

  return FirstBeadData(externalEnergies[selected].first, externalEnergies[selected].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

