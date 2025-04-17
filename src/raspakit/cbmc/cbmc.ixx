module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <optional>;
import <span>;
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

export namespace CBMC
{
// insertion
[[nodiscard]] std::optional<ChainData> growMoleculeSwapInsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta,
    Component::GrowType growType, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    size_t selectedComponent, size_t selectedMolecule, double scaling, size_t groupId,
    size_t numberOfTrialDirections) noexcept;

// deletion
[[nodiscard]] ChainData retraceMoleculeSwapDeletion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
    std::span<Atom> molecule, double scaling, size_t numberOfTrialDirections) noexcept;

// reinsertion grow
[[nodiscard]] std::optional<ChainData> growMoleculeReinsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, size_t numberOfTrialDirections) noexcept;

// reinsertion retrace
[[nodiscard]] ChainData retraceMoleculeReinsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, double storedR, size_t numberOfTrialDirections) noexcept;
}  // namespace CBMC
