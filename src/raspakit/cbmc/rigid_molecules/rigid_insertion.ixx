module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_rigid_insertion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import framework;
import component;
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] std::optional<ChainData> growRigidMoleculeSwapInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional, std::size_t numberOfTrialDirections) noexcept;
}

namespace CBMC
{
[[nodiscard]] std::optional<ChainData> growRigidMoleculeChainInsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t startingBead, std::vector<Atom> molecule, std::size_t numberOfTrialDirections,
    std::size_t selectedMolecule, double scaling, bool groupId, bool isFractional,
    const std::vector<Component> &components, std::size_t selectedComponent) noexcept;
}
