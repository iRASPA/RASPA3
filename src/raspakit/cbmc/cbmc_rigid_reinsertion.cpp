module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <tuple>
#include <optional>
#include <span>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <cmath>
#endif

module cbmc_rigid_reinsertion;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <iostream>;
import <algorithm>;
import <numeric>;
import <type_traits>;
import <cmath>;
#endif

import randomnumbers;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
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
    return ChainData(Molecule(double3(), simd_quatd()), {firstBeadData->atom}, firstBeadData->energies, 
                     firstBeadData->RosenbluthWeight, firstBeadData->storedR);
  }

  double scaling = 1.0;
  std::optional<ChainData> rigidRotationData = 
    growRigidMoleculeChainReinsertion(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, 
                                                          cutOffCoulomb, startingBead, atoms, numberOfTrialDirections,
                                                          selectedMolecule, scaling,
                                                          components, selectedComponent);
  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->molecule, rigidRotationData->atom, firstBeadData->energies + rigidRotationData->energies, 
                   firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, firstBeadData->storedR);
}


[[nodiscard]] std::optional<ChainData>
CBMC::growRigidMoleculeChainReinsertion(RandomNumber &random, bool hasExternalField,
                             const ForceField &forceField, const SimulationBox &simulationBox,
                             std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta,
                             double cutOff, double cutOffCoulomb, size_t startingBead,
                             std::vector<Atom> molecule, size_t numberOfTrialDirections,
                             size_t selectedMolecule, double scaling,
                             const std::vector<Component> &components, size_t selectedComponent) noexcept
{
  std::vector<std::pair<Molecule, std::vector<Atom>>> trialPositions{};

  // randomly rotated configurations around the starting bead
  for(size_t i = 0; i < numberOfTrialDirections; ++i)
  {
    simd_quatd orientation = random.randomSimdQuatd();
    std::vector<Atom> randomlyRotatedAtoms = components[selectedComponent].rotatePositions(orientation);
    double3 shift = molecule[startingBead].position - randomlyRotatedAtoms[startingBead].position;
    std::for_each(std::begin(randomlyRotatedAtoms), std::end(randomlyRotatedAtoms), [shift, selectedMolecule, scaling](Atom& atom) {
        atom.position += shift;
        atom.moleculeId = static_cast<uint32_t>(selectedMolecule);
        atom.groupId = static_cast<uint8_t>(0);
        atom.setScaling(scaling);
    });

    trialPositions.push_back({Molecule(shift, orientation), randomlyRotatedAtoms});
  };

  const std::vector<std::tuple<Molecule, std::vector<Atom>, RunningEnergy>>  externalEnergies =
    CBMC::computeExternalNonOverlappingEnergies(hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms,
                                cutOff, cutOffCoulomb, trialPositions, std::make_signed_t<std::size_t>(startingBead));
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(),
      std::back_inserter(logBoltmannFactors), [&](const std::tuple<Molecule, std::vector<Atom>, RunningEnergy>& v)
                                                  {return -beta * std::get<2>(v).total(); });

  size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
      [](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;


  return ChainData(std::get<0>(externalEnergies[selected]),
                   std::get<1>(externalEnergies[selected]),
                   std::get<2>(externalEnergies[selected]),
                   RosenbluthWeight / double(numberOfTrialDirections), 0.0);
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
    return ChainData(Molecule(double3(), simd_quatd()), std::vector<Atom>(molecule.begin(), molecule.end()), 
                     firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainData rigidRotationData = 
    retraceRigidChainReinsertion(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, 
                                 cutOffCoulomb, startingBead, molecule, numberOfTrialDirections);

  return ChainData(Molecule(double3(), simd_quatd()),
                   std::vector<Atom>(molecule.begin(), molecule.end()), 
                   firstBeadData.energies + rigidRotationData.energies, 
                   firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainData 
CBMC::retraceRigidChainReinsertion(RandomNumber &random, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
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

  return ChainData(Molecule(double3(), simd_quatd()), trialPositions[0], externalEnergies[0].second, 
                   RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

