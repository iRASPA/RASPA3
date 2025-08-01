module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

module cbmc_rigid_reinsertion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import component;
import atom;
import molecule;
import double3;
import simd_quatd;
import double3x3;
import simulationbox;
import energy_status;
import forcefield;
import energy_factor;
import running_energy;
import framework;
import component;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_rigid_insertion;
import cbmc_multiple_first_bead;
import interpolation_energy_grid;

[[nodiscard]] std::optional<ChainData> CBMC::growRigidMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, [[maybe_unused]] std::size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, std::size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> atoms = components[selectedComponent].copiedAtoms(molecule_atoms);
  std::size_t startingBead = components[selectedComponent].startingBead;

  std::optional<FirstBeadData> const firstBeadData = CBMC::growRigidMultipleFirstBeadReinsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, atoms[startingBead],
      numberOfTrialDirections);

  if (!firstBeadData) return std::nullopt;

  std::for_each(atoms.begin(), atoms.end(), [&](Atom &atom) { atom.position += firstBeadData->atom.position; });

  if (molecule_atoms.size() == 1)
  {
    return ChainData(Molecule(firstBeadData->atom.position, simd_quatd(), component.totalMass, component.componentId,
                              component.definedAtoms.size()),
                     {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight,
                     firstBeadData->storedR);
  }

  std::optional<ChainData> rigidRotationData = growRigidMoleculeChainReinsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, startingBead, molecule, atoms,
      components, selectedComponent, numberOfTrialDirections);
  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->molecule, rigidRotationData->atom,
                   firstBeadData->energies + rigidRotationData->energies,
                   firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, firstBeadData->storedR);
}

// helper function
[[nodiscard]] std::optional<ChainData> CBMC::growRigidMoleculeChainReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t startingBead, [[maybe_unused]] Molecule &molecule,
    std::vector<Atom> molecule_atoms, const std::vector<Component> &components, std::size_t selectedComponent,
    std::size_t numberOfTrialDirections) noexcept
{
  std::vector<std::pair<Molecule, std::vector<Atom>>> trialPositions{};

  // randomly rotated configurations around the starting bead
  for (std::size_t i = 0; i < numberOfTrialDirections; ++i)
  {
    simd_quatd orientation = random.randomSimdQuatd();
    std::vector<Atom> randomlyRotatedAtoms = components[selectedComponent].rotatePositions(orientation);
    double3 shift = molecule_atoms[startingBead].position - randomlyRotatedAtoms[startingBead].position;

    for (std::size_t j = 0; j < randomlyRotatedAtoms.size(); ++j)
    {
      randomlyRotatedAtoms[j].position += shift;
      randomlyRotatedAtoms[j].charge = molecule_atoms[j].charge;
      randomlyRotatedAtoms[j].scalingVDW = molecule_atoms[j].scalingVDW;
      randomlyRotatedAtoms[j].scalingCoulomb = molecule_atoms[j].scalingCoulomb;
      randomlyRotatedAtoms[j].moleculeId = molecule_atoms[j].moleculeId;
      randomlyRotatedAtoms[j].componentId = molecule_atoms[j].componentId;
      randomlyRotatedAtoms[j].groupId = molecule_atoms[j].groupId;
      randomlyRotatedAtoms[j].isFractional = molecule_atoms[j].isFractional;
    }

    trialPositions.push_back(
        {Molecule(shift, orientation, component.totalMass, component.componentId, component.definedAtoms.size()),
         randomlyRotatedAtoms});
  };

  const std::vector<std::tuple<Molecule, std::vector<Atom>, RunningEnergy>> externalEnergies =
      CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forceField, simulationBox,
                                                  interpolationGrids, framework, frameworkAtoms, moleculeAtoms,
                                                  cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions,
                                                  std::make_signed_t<std::size_t>(startingBead));

  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                 [&](const std::tuple<Molecule, std::vector<Atom>, RunningEnergy> &v)
                 { return -beta * std::get<2>(v).potentialEnergy(); });

  std::size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [](const double &acc, const double &logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return ChainData(std::get<0>(externalEnergies[selected]), std::get<1>(externalEnergies[selected]),
                   std::get<2>(externalEnergies[selected]), RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}

[[nodiscard]] ChainData CBMC::retraceRigidMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, [[maybe_unused]] std::size_t selectedComponent, [[maybe_unused]] std::size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, double storedR, std::size_t numberOfTrialDirections)
{
  std::size_t startingBead = components[selectedComponent].startingBead;

  const FirstBeadData firstBeadData = CBMC::retraceRigidMultipleFirstBeadReinsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms[startingBead], storedR,
      numberOfTrialDirections);

  if (molecule_atoms.size() == 1)
  {
    return ChainData(molecule, std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end()), firstBeadData.energies,
                     firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainData rigidRotationData = retraceRigidChainReinsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, startingBead, molecule, molecule_atoms,
      numberOfTrialDirections);

  return ChainData(molecule, std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end()),
                   firstBeadData.energies + rigidRotationData.energies,
                   firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainData CBMC::retraceRigidChainReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t startingBead, Molecule &molecule, std::span<Atom> molecule_atoms,
    std::size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<std::vector<Atom>> trialPositions = {trialPosition};

  for (std::size_t i = 1; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, trialPosition, startingBead));
  };

  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies =
      CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forceField, simulationBox,
                                                  interpolationGrids, framework, frameworkAtoms, moleculeAtoms,
                                                  cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions,
                                                  std::make_signed_t<std::size_t>(startingBead));

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<std::vector<Atom>, RunningEnergy> &v)
                 { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [](const double &acc, const double &logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  return ChainData(molecule, trialPositions[0], externalEnergies[0].second,
                   RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}
