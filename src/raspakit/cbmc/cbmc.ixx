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
[[nodiscard]] std::optional<ChainData> growMoleculeSwapInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule,
    double scaling, bool groupId, bool isFractional, std::size_t numberOfTrialDirections) noexcept;

// deletion
[[nodiscard]] ChainData retraceMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, std::span<Atom> molecule,
    double scaling, std::size_t numberOfTrialDirections) noexcept;

// reinsertion grow
[[nodiscard]] std::optional<ChainData> growMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, Molecule &molecule,
    std::span<Atom> molecule_atoms, std::size_t numberOfTrialDirections) noexcept;

// reinsertion retrace
[[nodiscard]] ChainData retraceMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, Molecule &molecule,
    std::span<Atom> molecule_atoms, double storedR, std::size_t numberOfTrialDirections) noexcept;
}  // namespace CBMC
