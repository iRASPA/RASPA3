module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_rigid_reinsertion;

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
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_interactions;
import framework;
import component;
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] std::optional<ChainData> growRigidMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, Molecule &molecule,
    std::span<Atom> molecule_atoms, size_t numberOfTrialDirections) noexcept;

[[nodiscard]] ChainData retraceRigidMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, [[maybe_unused]] size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, double storedR, size_t numberOfTrialDirections);
}  // namespace CBMC

namespace CBMC
{
[[nodiscard]] std::optional<ChainData> growRigidMoleculeChainReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, size_t startingBead, Molecule &molecule, std::vector<Atom> molecule_atoms,
    const std::vector<Component> &components, size_t selectedComponent, size_t numberOfTrialDirections) noexcept;

[[nodiscard]] ChainData retraceRigidChainReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, size_t startingBead, Molecule &molecule, std::span<Atom> molecule_atoms,
    size_t numberOfTrialDirections) noexcept;

}  // namespace CBMC
