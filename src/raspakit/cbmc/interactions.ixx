module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

export module cbmc_interactions;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import molecule;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import framework;
import component;
import double3;
import double3x3;
import forcefield;
import simulationbox;
import units;
import cbmc_interactions_intermolecular;
import cbmc_interactions_framework_molecule;
import interpolation_energy_grid;

export namespace CBMC
{
bool insideBlockedPockets(const std::optional<Framework> &frameworks, const Component &component,
                          std::span<const Atom> molecule_atoms);

[[nodiscard]] const std::vector<std::pair<Atom, RunningEnergy>> computeExternalNonOverlappingEnergies(
    const Component &component, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    std::vector<Atom> &trialPositions) noexcept;

const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> computeExternalNonOverlappingEnergies(
    const Component &component, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    std::vector<std::vector<Atom>> &trialPositionSets, std::make_signed_t<std::size_t> skip = -1) noexcept;

const std::vector<std::tuple<std::vector<Atom>, RunningEnergy, double>> computeExternalNonOverlappingEnergies(
    const Component &component, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    std::vector<std::vector<Atom>> &trialPositionSets, const std::vector<double> &RosenbluthWeightsTorsion,
    std::make_signed_t<std::size_t> skip = -1) noexcept;

const std::vector<std::tuple<Molecule, std::vector<Atom>, RunningEnergy>> computeExternalNonOverlappingEnergies(
    const Component &component, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    std::vector<std::pair<Molecule, std::vector<Atom>>> &trialPositionSets,
    std::make_signed_t<std::size_t> skip = -1) noexcept;

const std::optional<RunningEnergy> computeExternalNonOverlappingEnergyDualCutOff(
    const Component &component, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    std::vector<Atom> &trialPositionSet) noexcept;
}  // namespace CBMC
