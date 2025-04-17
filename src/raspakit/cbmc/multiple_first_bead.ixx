module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_multiple_first_bead;

#ifndef USE_LEGACY_HEADERS
import <optional>;
import <span>;
import <vector>;
#endif

import atom;
import randomnumbers;
import cbmc_first_bead_data;
import framework;
import component;
import forcefield;
import simulationbox;

export namespace CBMC
{
[[nodiscard]] std::optional<FirstBeadData> growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, size_t numberOfTrialDirections) noexcept;

[[nodiscard]] FirstBeadData retraceRigidMultipleFirstBeadSwapDeletion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom atom, double scaling, size_t numberOfTrialDirections) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growRigidMultipleFirstBeadReinsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, size_t numberOfTrialDirections) noexcept;

[[nodiscard]] FirstBeadData retraceRigidMultipleFirstBeadReinsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, const Atom &atom, double storedR, size_t numberOfTrialDirections);
}  // namespace CBMC
