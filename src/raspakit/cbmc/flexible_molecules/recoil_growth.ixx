module;

export module cbmc_recoil_growth;

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
// Recoil growth (RG) construction of a flexible molecule chain.
//
// Implements the recoil-growth algorithm for chain molecules with continuous interactions
// (S. Consta, T. J. H. Vlugt, J. Wichers Hoeth, B. Smit and D. Frenkel, Molecular Physics 97,
// 1243-1254 (1999)). The chain is grown segment by segment; at each segment 'k' trial directions
// are generated and a direction is considered 'open' with probability p = min(1, exp(-beta*u)).
// A direction is only used if, in addition to being open, a feeler of length 'l' (the recoil
// length) can be grown ahead of it. The Rosenbluth-like weight
//   W = prod_i [ m_i * max(1, exp(-beta*u_i)) * w_torsion_i ]
// is returned in ChainGrowData::RosenbluthWeight, where m_i is the number of open, feeler-viable
// trial directions at segment i and u_i is the (non-bonded) energy of the selected segment.
//
// The signature matches CBMC::growFlexibleMoleculeChainInsertion so it can be used as a drop-in
// replacement in the CBMC access routines.
[[nodiscard]] std::optional<ChainGrowData> growRecoilGrowthMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced, std::make_signed_t<std::size_t> skipBackgroundMolecule = -1);

// Recoil growth (RG) retrace of the existing flexible molecule chain.
//
// Computes the recoil-growth weight of the current (old) configuration. The existing chain occupies
// one trial direction per segment (always counted as open), and 'k-1' additional trial directions are
// generated to count m_i. The signature matches CBMC::retraceFlexibleMoleculeChainDeletion.
[[nodiscard]] ChainRetraceData retraceRecoilGrowthMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced) noexcept;
}  // namespace CBMC
