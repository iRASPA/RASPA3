module;

export module cbmc;

import std;

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
import cbmc_util;

export namespace CBMC
{
// insertion
[[nodiscard]] std::optional<ChainGrowData> growMoleculeSwapInsertion(
    RandomNumber &random, Component &component, std::size_t selectedComponent, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, double scaling, std::uint8_t groupId,
    bool isFractional) noexcept;

// deletion
[[nodiscard]] ChainRetraceData retraceMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept;

// reinsertion grow
[[nodiscard]] std::optional<ChainGrowData> growMoleculeReinsertion(
    RandomNumber &random, Component &component, std::size_t selectedComponent, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms) noexcept;

// reinsertion retrace
[[nodiscard]] std::optional<ChainRetraceData> retraceMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms,
    double storedR) noexcept;

// partial reinsertion grow
[[nodiscard]] std::optional<ChainGrowData> growMoleculePartialReinsertion(
    RandomNumber &random, Component &component, std::size_t selectedComponent, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept;

// partial reinsertion retrace
[[nodiscard]] ChainRetraceData retraceMoleculePartialReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept;

// identity change insertion
[[nodiscard]] std::optional<ChainGrowData> growMoleculeIdentityChangeInsertion(
    RandomNumber &random, Component &component, std::size_t selectedComponent, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, const Atom &oldStartingBead,
    double scaling, std::uint8_t groupId, bool isFractional,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

// identity change deletion
[[nodiscard]] ChainRetraceData retraceMoleculeIdentityChangeDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept;

// distance-biased ion-pair insertion: second molecule with fixed first-bead position
[[nodiscard]] std::optional<ChainGrowData> growMoleculePairSecondSwapInsertion(
    RandomNumber &random, Component &component, std::size_t selectedComponent, bool hasExternalField,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, double3 fixedFirstBeadPosition,
    double scaling, std::uint8_t groupId, bool isFractional) noexcept;

[[nodiscard]] ChainRetraceData retraceMoleculePairSecondSwapDeletion(
    const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept;
}  // namespace CBMC
