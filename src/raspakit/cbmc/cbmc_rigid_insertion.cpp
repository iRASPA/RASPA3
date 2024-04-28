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

module cbmc_rigid_insertion;

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


// atoms is a recentered copy of the molecule (recentered around the starting bead)
[[nodiscard]] std::optional<ChainData> 
CBMC::growRigidMoleculeSwapInsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                     const ForceField &forceField, const SimulationBox &simulationBox, 
                                     std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                     double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                     size_t selectedMolecule, double scaling, 
                                     size_t numberOfTrialDirections) noexcept
{
  size_t startingBead = components[selectedComponent].startingBead;
  Atom firstBead = components[selectedComponent].atoms[startingBead];
  firstBead.moleculeId = static_cast<uint32_t>(selectedMolecule);
  firstBead.groupId = static_cast<uint8_t>(0);
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData = 
    CBMC::growMoleculeMultipleFirstBeadSwapInsertion(random, hasExternalField, forceField, simulationBox,
                                                frameworkAtoms, moleculeAtoms, beta, cutOff, cutOffCoulomb, 
                                                firstBead, numberOfTrialDirections);

  if (!firstBeadData) return std::nullopt;

  if(components[selectedComponent].atoms.size() == 1)
  {
    return ChainData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0)), 
                     {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> atoms = components[selectedComponent].atoms;
  std::for_each(atoms.begin(), atoms.end(), [&](Atom& atom) {
      atom.position += firstBeadData->atom.position - components[selectedComponent].atoms[startingBead].position;
      atom.moleculeId = static_cast<uint32_t>(selectedMolecule);
      atom.groupId = static_cast<uint8_t>(0);
      atom.setScaling(scaling);
    });


  std::optional<ChainData> const rigidRotationData = 
    CBMC::growRigidMoleculeChainInsertion(random, hasExternalField, forceField, simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOff, 
                           cutOffCoulomb, startingBead, atoms, numberOfTrialDirections,
                           selectedMolecule, scaling,
                           components, selectedComponent);
  
  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->molecule, rigidRotationData->atom, firstBeadData->energies + rigidRotationData->energies, 
                   firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, 0.0);
}


[[nodiscard]] std::optional<ChainData> 
CBMC::growRigidMoleculeChainInsertion(RandomNumber &random, bool hasExternalField,
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
