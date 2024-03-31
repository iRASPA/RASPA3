module;

module cbmc_rigid_reinsertion;

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
import forcefield;
import energy_factor;
import running_energy;
import component;
import cbmc_growing_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_rigid_insertion;
import cbmc_multiple_first_bead;


[[nodiscard]] ChainData 
retraceRigidChainReinsertion(RandomNumber &random, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
                             std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                             double cutOff, double cutOffCoulomb, size_t startingBead, std::span<Atom> molecule, 
                             size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule.begin(), molecule.end());
  std::vector<std::vector<Atom>> trialPositions = { trialPosition };

  for (size_t i = 1; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, trialPosition, startingBead));
  };

  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies = 
    CBMC::computeExternalNonOverlappingEnergies(hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, cutOff, 
                                        cutOffCoulomb, trialPositions, std::make_signed_t<std::size_t>(startingBead));

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies),
      std::back_inserter(logBoltmannFactors), [&](const std::pair<std::vector<Atom>, RunningEnergy>& v) 
                                                  {return -beta * v.second.total(); });

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  return ChainData(trialPositions[0], externalEnergies[0].second, 
                   RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

[[nodiscard]] std::optional<ChainData> 
CBMC::growRigidMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                   const ForceField &forceField, const SimulationBox &simulationBox, 
                                   std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                   double beta,  double cutOff, double cutOffCoulomb, size_t selectedComponent,
                                   [[maybe_unused]] size_t selectedMolecule, std::span<Atom> molecule, 
                                   size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> atoms = components[selectedComponent].copiedAtoms(molecule);
  size_t startingBead = components[selectedComponent].startingBead; 

  std::optional<FirstBeadData> const firstBeadData = 
    CBMC::growRigidMultipleFirstBeadReinsertion(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, 
                                      beta, cutOff, cutOffCoulomb, atoms[startingBead], numberOfTrialDirections);

  if (!firstBeadData) return std::nullopt;

  std::for_each(atoms.begin(), atoms.end(), [&](Atom& atom) {atom.position += firstBeadData->atom.position; });

  if(molecule.size() == 1)
  {
    return ChainData({firstBeadData->atom}, firstBeadData->energies, 
                     firstBeadData->RosenbluthWeight, firstBeadData->storedR);
  }

  std::optional<ChainData> rigidRotationData = 
    growRigidMoleculeChain(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, 
                                                          cutOffCoulomb, startingBead, atoms, numberOfTrialDirections);
  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->atom, firstBeadData->energies + rigidRotationData->energies, 
                   firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, firstBeadData->storedR);
}

[[nodiscard]] ChainData 
CBMC::retraceRigidMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                const ForceField &forceField, const SimulationBox &simulationBox, 
                                std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                double beta, double cutOff, double cutOffCoulomb, 
                                [[maybe_unused]] size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule, 
                                std::span<Atom> molecule, double storedR, size_t numberOfTrialDirections)
{
  size_t startingBead = components[selectedComponent].startingBead;

  const FirstBeadData firstBeadData = 
    CBMC::retraceRigidMultipleFirstBeadReinsertion(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, 
                             beta, cutOff, cutOffCoulomb, molecule[startingBead], storedR, numberOfTrialDirections);

  if(molecule.size() == 1)
  {
    return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), 
                     firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainData rigidRotationData = 
    retraceRigidChainReinsertion(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, 
                                 cutOffCoulomb, startingBead, molecule, numberOfTrialDirections);

  return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), 
                   firstBeadData.energies + rigidRotationData.energies, 
                   firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}
