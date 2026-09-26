module;

export module cbmc_first_bead;

import std;

import atom;
import randomnumbers;
import cbmc_results;
import component;
import cbmc_grow_context;

// The first-bead stage of a grow or retrace: the two dispatchers over 'CBMC::FirstBeadScheme' (defined
// in 'cbmc_results'). The individual schemes are internal to the implementation. A grow returns
// std::nullopt when every trial position overlaps (or the weight falls below
// 'minimumRosenbluthFactor'); a retrace throws std::runtime_error when the existing first bead
// overlaps, since an accepted configuration can not overlap (see the error contract in the 'cbmc'
// module).

export namespace CBMC
{
/**
 * \brief Places the first bead with 'scheme'. 'firstBead' carries the identity and scaling attributes
 * of the molecule and, for the pinned and fixed schemes, the position.
 *
 * Throws std::invalid_argument for 'AlreadyPlaced', which has no first-bead stage.
 */
[[nodiscard]] std::optional<FirstBeadData> growFirstBead(RandomNumber &random, const GrowContext &context,
                                                         const Component &component, const Atom &firstBead,
                                                         FirstBeadScheme scheme);

/**
 * \brief Weighs the existing first bead with 'scheme'. 'storedR' is the retained partial weight of the
 * matching grow ('Reinsertion' only, ignored otherwise).
 *
 * Throws std::invalid_argument for 'AlreadyPlaced', std::runtime_error when the bead overlaps.
 */
[[nodiscard]] FirstBeadData retraceFirstBead(RandomNumber &random, const GrowContext &context,
                                             const Component &component, const Atom &firstBead,
                                             FirstBeadScheme scheme, double storedR);
}  // namespace CBMC
