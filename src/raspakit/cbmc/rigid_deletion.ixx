module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_rigid_deletion;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <optional>;
import <span>;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_interactions;
import framework;
import component;

export namespace CBMC
{
[[nodiscard]] ChainData retraceRigidMoleculeSwapDeletion(
    RandomNumber &random, const Framework &framework, const Component &component,
    bool hasExternalField, const std::vector<Component> &components, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
    double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent,
    size_t selectedMolecule, std::span<Atom> molecule, double scaling, size_t numberOfTrialDirections) noexcept;
}
