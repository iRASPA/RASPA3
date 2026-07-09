module;

export module cbmc_rigid_insertion;

import std;

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import cbmc_growth_context;
import framework;
import component;
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] std::optional<ChainGrowData> growRigidMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::size_t selectedComponent,
    std::span<Atom> molecule_atoms, std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;
}
