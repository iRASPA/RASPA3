module;

module cbmc_rigid_deletion;

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
import energy_status;
import running_energy;
import component;
import cbmc_growing_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_multiple_first_bead;


[[nodiscard]] ChainData                                                                                                 
retraceRigidChain(RandomNumber &random, bool hasExternalField, const ForceField &forcefield, const SimulationBox &simulationBox,               
                  std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta,               
                  double cutOff, double cutOffCoulomb, size_t startingBead, [[maybe_unused]] double scaling,            
                  std::span<Atom> molecule, size_t numberOfTrialDirections) noexcept;


[[nodiscard]] ChainData 
CBMC::retraceRigidMoleculeSwapDeletion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                 const ForceField &forcefield, const SimulationBox &simulationBox, 
                                 std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                 double beta, double cutOff, double cutOffCoulomb, 
                                 [[maybe_unused]] size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule, 
                                 std::span<Atom> molecule, double scaling, double storedR, size_t 
                                 numberOfTrialDirections) noexcept
{
  size_t startingBead = components[selectedComponent].startingBead;

  const FirstBeadData firstBeadData = 
    CBMC::retraceRigidMultipleFirstBeadSwapDeletion(random, hasExternalField, forcefield, simulationBox, frameworkAtoms, moleculeAtoms, 
                      beta, cutOff, cutOffCoulomb, molecule[startingBead], scaling, storedR, numberOfTrialDirections);

  if(molecule.size() == 1)
  {
    return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), 
                     firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  const ChainData rigidRotationData = 
    retraceRigidChain(random, hasExternalField, forcefield, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, cutOffCoulomb, 
                      startingBead, scaling, molecule, numberOfTrialDirections);

  return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), 
                                     firstBeadData.energies + rigidRotationData.energies,
                                     firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}


[[nodiscard]] ChainData 
retraceRigidChain(RandomNumber &random, bool hasExternalField, const ForceField &forcefield, const SimulationBox &simulationBox, 
                  std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                  double cutOff, double cutOffCoulomb, size_t startingBead, [[maybe_unused]] double scaling, 
                  std::span<Atom> molecule, size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule.begin(), molecule.end());
  std::for_each(trialPosition.begin(), trialPosition.end(), [&](Atom& a) {a.setScaling(scaling); });
  std::vector<std::vector<Atom>> trialPositions = { trialPosition };

  for (size_t i = 1; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, trialPosition, startingBead));
  };

  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies = 
    CBMC::computeExternalNonOverlappingEnergies(hasExternalField, forcefield, simulationBox, frameworkAtoms, moleculeAtoms, cutOff, 
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
