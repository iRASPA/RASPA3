module;

export module cbmc_chain_cbmc;

import std;

import atom;
import randomnumbers;
import component;
import cbmc_results;
import cbmc_growth_context;

// Configurational-bias Monte Carlo growth of the chain beyond the first bead, and the retrace of an
// existing chain. Both walk the same deterministic growth plan ('Component::growthPlan') step by
// step: the operator engine generates the trial directions of a step, the external energies weigh
// them, one is Rosenbluth-selected (the grow) or the old one is kept (the retrace), and the per-step
// factor is accumulated in log space. The recoil-growth alternative is in cbmc_chain_recoil.
export namespace CBMC
{
/**
 * \brief Grows the remaining beads of a molecule with CBMC.
 *
 * 'molecule_atoms' carries the placed beads ('beadsAlreadyPlaced', at least the first bead); the
 * returned molecule has every bead placed. Returns std::nullopt when some step has no non-overlapping
 * trial direction or its factor falls below 'minimumRosenbluthFactor' (an ordinary rejection).
 */
[[nodiscard]] std::optional<ChainGrowData> growFlexibleMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced, std::make_signed_t<std::size_t> skipBackgroundMolecule = -1);

/**
 * \brief The CBMC Rosenbluth weight of an existing molecule beyond its placed beads.
 *
 * Throws std::runtime_error when the existing configuration overlaps: an accepted state can not
 * overlap, so this signals an inconsistent simulation state rather than silently assigning a weight
 * (see the error contract in the 'cbmc' module).
 */
[[nodiscard]] ChainRetraceData retraceFlexibleMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced);
}  // namespace CBMC
