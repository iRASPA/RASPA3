module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <numbers>
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_flexible_deletion;

#ifdef USE_STD_IMPORT
import std;
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
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] ChainRetraceData retraceFlexibleMoleculeChainDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forcefield,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::span<Atom> molecule_atoms, const std::vector<std::size_t> beadsAlreadyPlaced) noexcept;
}
