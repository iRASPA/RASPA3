module;

module cbmc_rigid_insertion;

import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <iostream>;
import <algorithm>;
import <numeric>;
import <type_traits>;

import randomnumbers;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import energy_status;
import cbmc_growing_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import forcefield;
import energy_factor;
import running_energy;
import component;


[[nodiscard]] std::optional<ChainData>                                                                                
growRigidMoleculeChain(RandomNumber &random, bool hasExternalField,  const ForceField &forceField, const SimulationBox &simulationBox,        
                       std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,                     
                       double beta, double cutOff, double cutOffCoulomb, size_t startingBead,                         
                       std::vector<Atom> molecule, size_t numberOfTrialDirections) noexcept;   

[[nodiscard]] std::optional<ChainData> 
CBMC::growRigidMoleculeSwapInsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                     const ForceField &forceField, const SimulationBox &simulationBox, 
                                     std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                     double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                     [[maybe_unused]] size_t selectedMolecule, double scaling, 
                                     std::vector<Atom> atoms, size_t numberOfTrialDirections) noexcept
{
  for (Atom& atom : atoms)
  {
    atom.setScaling(scaling);
  }
  size_t startingBead = components[selectedComponent].startingBead;

  std::optional<FirstBeadData> 
    const firstBeadData = CBMC::growMoleculeMultipleFirstBeadSwapInsertion(random, hasExternalField, forceField, simulationBox,
                                                frameworkAtoms, moleculeAtoms, beta, cutOff, cutOffCoulomb, 
                                                atoms[startingBead], numberOfTrialDirections);

  if (!firstBeadData) return std::nullopt;

  std::for_each(atoms.begin(), atoms.end(), [&](Atom& atom) {atom.position += firstBeadData->atom.position; });

  if(atoms.size() == 1)
  {
    return ChainData({firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  std::optional<ChainData> const rigidRotationData = 
    CBMC::growRigidMoleculeChain(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, 
                           cutOffCoulomb, startingBead, atoms, numberOfTrialDirections);
  
  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->atom, firstBeadData->energies + rigidRotationData->energies, 
                   firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, 0.0);
}


// helper function
[[nodiscard]] std::optional<ChainData> 
CBMC::growRigidMoleculeChain(RandomNumber &random, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                             std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                             double cutOff, double cutOffCoulomb, size_t startingBead, 
                             std::vector<Atom> molecule, size_t numberOfTrialDirections) noexcept
{
  std::vector<std::vector<Atom>> trialPositions{};

  for(size_t i = 0; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, molecule, startingBead));
  };
  
  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>>  externalEnergies = 
    CBMC::computeExternalNonOverlappingEnergies(hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, 
                                cutOff, cutOffCoulomb, trialPositions, std::make_signed_t<std::size_t>(startingBead));
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(),
      std::back_inserter(logBoltmannFactors), [&](const std::pair<std::vector<Atom>, RunningEnergy>& v) 
                                                  {return -beta * v.second.total(); });

  size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return ChainData(externalEnergies[selected].first, externalEnergies[selected].second, 
                   RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}
