module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

module cbmc_rigid_insertion;

#ifndef USE_LEGACY_HEADERS
import std;
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
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import forcefield;
import energy_factor;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;

// atoms is a recentered copy of the molecule (recentered around the starting bead)
[[nodiscard]] std::optional<ChainData> CBMC::growRigidMoleculeSwapInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional, std::size_t numberOfTrialDirections) noexcept
{
  std::size_t startingBead = components[selectedComponent].startingBead;
  Atom firstBead = components[selectedComponent].atoms[startingBead];
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, firstBead, numberOfTrialDirections);

  if (!firstBeadData) return std::nullopt;

  if (components[selectedComponent].atoms.size() == 1)
  {
    return ChainData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                              component.totalMass, component.componentId, component.definedAtoms.size()),
                     {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> atoms = components[selectedComponent].atoms;
  std::for_each(atoms.begin(), atoms.end(),
                [&](Atom &atom)
                {
                  atom.position +=
                      firstBeadData->atom.position - components[selectedComponent].atoms[startingBead].position;
                  atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
                  atom.groupId = groupId;
                  atom.isFractional = isFractional;
                  atom.setScaling(scaling);
                });

  std::optional<ChainData> const rigidRotationData = CBMC::growRigidMoleculeChainInsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, startingBead, atoms,
      numberOfTrialDirections, selectedMolecule, scaling, groupId, isFractional, components, selectedComponent);

  if (!rigidRotationData) return std::nullopt;

  return ChainData(rigidRotationData->molecule, rigidRotationData->atom,
                   firstBeadData->energies + rigidRotationData->energies,
                   firstBeadData->RosenbluthWeight * rigidRotationData->RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainData> CBMC::growRigidMoleculeChainInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t startingBead, std::vector<Atom> molecule, std::size_t numberOfTrialDirections,
    std::size_t selectedMolecule, double scaling, bool groupId, bool isFractional,
    const std::vector<Component> &components, std::size_t selectedComponent) noexcept
{
  std::vector<std::pair<Molecule, std::vector<Atom>>> trialPositions{};

  // randomly rotated configurations around the starting bead
  for (std::size_t i = 0; i < numberOfTrialDirections; ++i)
  {
    simd_quatd orientation = random.randomSimdQuatd();
    std::vector<Atom> randomlyRotatedAtoms = components[selectedComponent].rotatePositions(orientation);
    double3 shift = molecule[startingBead].position - randomlyRotatedAtoms[startingBead].position;
    std::for_each(std::begin(randomlyRotatedAtoms), std::end(randomlyRotatedAtoms),
                  [shift, selectedMolecule, scaling, groupId, isFractional](Atom &atom)
                  {
                    atom.position += shift;
                    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
                    atom.groupId = groupId;
                    atom.isFractional = isFractional;
                    atom.setScaling(scaling);
                  });

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
