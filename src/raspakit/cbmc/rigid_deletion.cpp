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

module cbmc_rigid_deletion;

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
import energy_status;
import running_energy;
import framework;
import component;
import cbmc_growing_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_multiple_first_bead;
import interpolation_energy_grid;

[[nodiscard]] ChainData retraceRigidChain(RandomNumber &random, const Component &component, bool hasExternalField,
                                          const ForceField &forcefield, const SimulationBox &simulationBox,
                                          const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
                                          const std::optional<Framework> &framework,
                                          std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                          double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
                                          double cutOffCoulomb, size_t startingBead, [[maybe_unused]] double scaling,
                                          std::span<Atom> molecule, size_t numberOfTrialDirections) noexcept;

[[nodiscard]] ChainData CBMC::retraceRigidMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forcefield, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, [[maybe_unused]] size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule,
    std::span<Atom> molecule, double scaling, size_t numberOfTrialDirections) noexcept
{
  size_t startingBead = components[selectedComponent].startingBead;

  const FirstBeadData firstBeadData = CBMC::retraceRigidMultipleFirstBeadSwapDeletion(
      random, component, hasExternalField, forcefield, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule[startingBead], scaling,
      numberOfTrialDirections);

  if (molecule.size() == 1)
  {
    return ChainData(
        Molecule(double3(), simd_quatd(), component.totalMass, component.componentId, component.definedAtoms.size()),
        std::vector<Atom>(molecule.begin(), molecule.end()), firstBeadData.energies, firstBeadData.RosenbluthWeight,
        0.0);
  }

  const ChainData rigidRotationData =
      retraceRigidChain(random, component, hasExternalField, forcefield, simulationBox, interpolationGrids, framework,
                        frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
                        startingBead, scaling, molecule, numberOfTrialDirections);

  return ChainData(
      Molecule(double3(), simd_quatd(), component.totalMass, component.componentId, component.definedAtoms.size()),
      std::vector<Atom>(molecule.begin(), molecule.end()), firstBeadData.energies + rigidRotationData.energies,
      firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainData retraceRigidChain(RandomNumber &random, const Component &component, bool hasExternalField,
                                          const ForceField &forcefield, const SimulationBox &simulationBox,
                                          const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
                                          const std::optional<Framework> &framework,
                                          std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                          double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
                                          double cutOffCoulomb, size_t startingBead, [[maybe_unused]] double scaling,
                                          std::span<Atom> molecule, size_t numberOfTrialDirections) noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule.begin(), molecule.end());
  std::for_each(trialPosition.begin(), trialPosition.end(), [&](Atom &a) { a.setScaling(scaling); });
  std::vector<std::vector<Atom>> trialPositions = {trialPosition};

  for (size_t i = 1; i < numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, trialPosition, startingBead));
  };

  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies =
      CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forcefield, simulationBox,
                                                  interpolationGrids, framework, frameworkAtoms, moleculeAtoms,
                                                  cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions,
                                                  std::make_signed_t<std::size_t>(startingBead));

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<std::vector<Atom>, RunningEnergy> &v)
                 { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight =
      std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                  [](const double &acc, const double &logBoltmannFactor) { return acc + std::exp(logBoltmannFactor); });

  return ChainData(
      Molecule(double3(), simd_quatd(), component.totalMass, component.componentId, component.definedAtoms.size()),
      trialPositions[0], externalEnergies[0].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}
