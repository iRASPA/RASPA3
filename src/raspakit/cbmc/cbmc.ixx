module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import molecule;
import double3x3;
import double3;
import randomnumbers;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import cbmc_chain_data;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;

export namespace CBMC
{
// insertion
[[nodiscard]] std::optional<ChainGrowData> growMoleculeSwapInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional) noexcept;

// deletion
[[nodiscard]] ChainRetraceData retraceMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept;

// reinsertion grow
[[nodiscard]] std::optional<ChainGrowData> growMoleculeReinsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms) noexcept;

// reinsertion retrace
[[nodiscard]] ChainRetraceData retraceMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms,
    double storedR) noexcept;

// partial reinsertion grow
[[nodiscard]] std::optional<ChainGrowData> growMoleculePartialReinsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept;

// partial reinsertion retrace
[[nodiscard]] ChainRetraceData retraceMoleculePartialReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept;

// identity change insertion
[[nodiscard]] std::optional<ChainGrowData> growMoleculeIdentityChangeInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional) noexcept;

// identity change deletion
[[nodiscard]] ChainRetraceData retraceMoleculeIdentityChangeDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept;
}  // namespace CBMC
