module;

export module cbmc_rigid_deletion;

import std;

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_interactions;
import cbmc_growth_context;
import framework;
import component;
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] ChainRetraceData retraceRigidMoleculeChainDeletion(RandomNumber &random, const GrowContext &context,
                                                                 const Component &component,
                                                                 std::span<Atom> molecule_atoms) noexcept;
}
