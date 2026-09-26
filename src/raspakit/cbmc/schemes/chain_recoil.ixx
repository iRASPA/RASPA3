module;

export module cbmc_chain_recoil;

import std;

import atom;
import randomnumbers;
import component;
import cbmc_results;
import cbmc_grow_context;

export namespace CBMC
{
// Recoil growth (RG) construction of a flexible molecule chain.
//
// Implements the recoil-growth algorithm for chain molecules with continuous interactions
// (S. Consta, T. J. H. Vlugt, J. Wichers Hoeth, B. Smit and D. Frenkel, Molecular Physics 97,
// 1243-1254 (1999)). The chain is grown segment by segment; at each segment 'k' trial directions
// are generated and a direction is considered 'open' with probability
//   p_open = min(1, exp(-beta*(u - u_ref)))
// with u_ref a fixed per-step reference energy (see 'openProbability' in the implementation). A
// direction is only used if, in addition to being open, a feeler of length 'l' (the recoil length)
// can be grown ahead of it. The recoil-growth weight
//   W = prod_i [ (m_i / k) * exp(-beta*u_i) / p_open,i * w_torsion_i * exp(-beta*u_unsampled,i) ]
// is returned as its logarithm in CBMC::GrowResult::logRosenbluthWeight, where m_i is the number of
// available (open and feeler-viable) trial directions at segment i and u_i is the non-bonded energy
// of the selected segment.
//
// Two properties of this weight matter to callers:
//  - It is established as a valid factor of a Metropolis acceptance RATIO W_new / W_old (grow and
//    retrace use identical feelers, so the feeler probabilities cancel in the super-detailed-balance
//    argument of Consta et al.). It is NOT the Rosenbluth weight whose ensemble average is the Widom
//    estimator of the excess chemical potential: Widom sampling always grows with configurational
//    bias ('CBMC::ChainScheme').
//  - The retrace divides by the openness probability of the existing configuration. In a crowded,
//    repulsive environment an accepted configuration can sit well above the reference energy at a
//    step, p_open is then small and W_old large: recoil growth has a higher weight variance than
//    configurational bias there. This is inherent to the open/closed formulation for continuous
//    potentials; the per-step reference (the maximum strain of the ideal-gas reference conformations)
//    removes the molecule's own intrinsic strain from that variance but not the external one.
//
// Same signature as CBMC::growChainCBMC; the 'cbmc' module dispatches on the context's chain scheme.
[[nodiscard]] std::optional<CBMC::GrowResult> growChainRecoil(RandomNumber &random, const GrowContext &context,
                                                              const Component &component,
                                                              std::span<const Atom> moleculeAtoms,
                                                              const std::vector<std::size_t> &beadsAlreadyPlaced);

// Recoil growth (RG) retrace of the existing flexible molecule chain.
//
// Computes the recoil-growth weight of the current (old) configuration. The existing chain occupies
// one trial direction per segment (always counted as open), and 'k-1' additional trial directions are
// generated to count m_i. Same signature as CBMC::retraceChainCBMC.
//
// Throws std::runtime_error when the existing configuration overlaps (hard-core overlap, blocked
// pocket): an accepted state can not overlap, so this signals an inconsistent simulation state rather
// than silently assigning the molecule a weight.
[[nodiscard]] CBMC::RetraceResult retraceChainRecoil(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<const Atom> moleculeAtoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced);
}  // namespace CBMC
