module;

export module cbmc_first_bead;

import std;

import atom;
import randomnumbers;
import cbmc_results;
import component;
import cbmc_grow_context;

// The first-bead schemes. A grow returns std::nullopt when every trial position overlaps (or the
// weight falls below 'minimumRosenbluthFactor'); a retrace throws std::runtime_error when the
// existing first bead overlaps, since an accepted configuration can not overlap (see the error
// contract in the 'cbmc' module).

export namespace CBMC
{
/// Multiple-first-bead scheme: 'numberOfFirstBeadPositions' uniformly random positions in the box;
/// weight = sum of Boltzmann factors / number of positions.
[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBead(RandomNumber &random, const GrowContext &context,
                                                                 const Component &component, const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceMultipleFirstBead(RandomNumber &random, const GrowContext &context,
                                                     const Component &component, const Atom &atom);

/// Multiple-first-bead reinsertion (Esselink et al.): as above against a background without the
/// molecule itself, retaining the partial weight 'storedR' that the retrace needs.
[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(RandomNumber &random,
                                                                            const GrowContext &context,
                                                                            const Component &component,
                                                                            const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceMultipleFirstBeadReinsertion(const GrowContext &context, const Component &component,
                                                                const Atom &atom, double storedR);

/// Pinned first bead: a single trial at the given position, its Boltzmann factor as weight (identity
/// change: the new molecule takes the position of the old one).
[[nodiscard]] std::optional<FirstBeadData> growPinnedFirstBead(const GrowContext &context, const Component &component,
                                                               const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retracePinnedFirstBead(const GrowContext &context, const Component &component,
                                                   const Atom &atom);

/// Fixed first bead: a single trial at the given position with weight one (the caller sampled the
/// position and accounts for its bias, e.g. the distance-biased pair insertion).
[[nodiscard]] std::optional<FirstBeadData> growFixedFirstBead(const GrowContext &context, const Component &component,
                                                              const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceFixedFirstBead(const GrowContext &context, const Component &component,
                                                  const Atom &atom);
}  // namespace CBMC
