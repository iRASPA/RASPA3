module;

export module cbmc_multiple_first_bead;

import std;

import atom;
import randomnumbers;
import cbmc_first_bead_data;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;
import cbmc_util;

export namespace CBMC
{
[[nodiscard]] std::optional<FirstBeadData> growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceMultipleFirstBeadSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::optional<FirstBeadData> retraceMultipleFirstBeadReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, double storedR,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadPartialInsertion(
    const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] FirstBeadData retraceMultipleFirstBeadPartialDeletion(
    const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growFirstBeadAtFixedPosition(
    const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceFirstBeadAtFixedPosition(
    const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom) noexcept;
}  // namespace CBMC
