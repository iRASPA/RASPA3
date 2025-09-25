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

[[nodiscard]] std::optional<ChainGrowData> CBMC::growRigidMoleculeChainInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept
{
  std::vector<std::pair<Molecule, std::vector<Atom>>> trialPositions(forceField.numberOfTrialDirections);

  std::size_t starting_bead = component.startingBead;

  // randomly rotated configurations around the starting bead
  for (std::size_t i = 0; i < forceField.numberOfTrialDirections; ++i)
  {
    simd_quatd orientation = random.randomSimdQuatd();
    std::vector<Atom> randomlyRotatedAtoms = CBMC::rotateRandomlyAround(orientation, molecule_atoms, starting_bead);
    double3 com = component.computeCenterOfMass(randomlyRotatedAtoms);

    trialPositions[i] = {
        Molecule(com, orientation, component.totalMass, component.componentId, component.definedAtoms.size()),
        randomlyRotatedAtoms};
  };

  const std::vector<std::tuple<Molecule, std::vector<Atom>, RunningEnergy>> externalEnergies =
      CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forceField, simulationBox,
                                                  interpolationGrids, framework, frameworkAtomData, moleculeAtomData,
                                                  cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions,
                                                  std::make_signed_t<std::size_t>(component.startingBead));
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                 [&](const std::tuple<Molecule, std::vector<Atom>, RunningEnergy> &v)
                 { return -beta * std::get<2>(v).potentialEnergy(); });

  std::size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

  auto [selected_molecule, selected_atoms, selected_running_energy] = externalEnergies[selected];

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [](const double &acc, const double &logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  if (RosenbluthWeight < forceField.minimumRosenbluthFactor) return std::nullopt;

  return ChainGrowData(selected_molecule, selected_atoms, selected_running_energy,
                       RosenbluthWeight / double(forceField.numberOfTrialDirections), 0.0);
}
