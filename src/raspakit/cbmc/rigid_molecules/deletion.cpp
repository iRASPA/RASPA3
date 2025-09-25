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
import forcefield;
import energy_factor;
import energy_status;
import running_energy;
import framework;
import component;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_multiple_first_bead;
import interpolation_energy_grid;

[[nodiscard]] ChainRetraceData CBMC::retraceRigidMoleculeChainDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<std::vector<Atom>> trialPositions = {trialPosition};

  for (std::size_t i = 1; i < forceField.numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, trialPosition, component.startingBead));
  };

  const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies =
      CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forceField, simulationBox,
                                                  interpolationGrids, framework, frameworkAtomData, moleculeAtomData,
                                                  cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions,
                                                  std::make_signed_t<std::size_t>(component.startingBead));

  std::vector<double> logBoltmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
                 [&](const std::pair<std::vector<Atom>, RunningEnergy> &v)
                 { return -beta * v.second.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                            [](const double &acc, const double &logBoltmannFactor)
                                            { return acc + std::exp(logBoltmannFactor); });

  return ChainRetraceData(externalEnergies[0].second, RosenbluthWeight / double(forceField.numberOfTrialDirections),
                          0.0);
}
