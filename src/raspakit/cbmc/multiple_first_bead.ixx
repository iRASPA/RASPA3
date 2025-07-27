module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_multiple_first_bead;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import randomnumbers;
import cbmc_first_bead_data;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] std::optional<FirstBeadData> growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, std::size_t numberOfTrialDirections) noexcept;

[[nodiscard]] FirstBeadData retraceRigidMultipleFirstBeadSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forcefield,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom, double scaling, std::size_t numberOfTrialDirections) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growRigidMultipleFirstBeadReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, std::size_t numberOfTrialDirections) noexcept;

[[nodiscard]] FirstBeadData retraceRigidMultipleFirstBeadReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, double storedR, std::size_t numberOfTrialDirections);
}  // namespace CBMC
