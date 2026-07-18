module;

export module cbmc_segments;

import std;

import atom;
import component;
import intra_molecular_potentials;

export namespace CBMC
{
/**
 * \brief A single growth step of the chain.
 *
 * A segment places a set of 'nextBeads' grown from 'currentBead', optionally with a 'previousBead'
 * defining the bend/torsion reference. There are two kinds of segments:
 *  - Flexible: 'nextBeads' are individual beads (a branch point can place several at once) that are
 *    sampled bead-by-bead with the bond/bend/torsion Monte-Carlo scheme.
 *  - Rigid unit: 'nextBeads' are all the (still unplaced) atoms of a rigid group and are placed as
 *    one rigid body. 'rigidGroup' holds the index of that group in 'Component::groups'.
 *
 * The growth sequence is deterministic (it only depends on the connectivity table and the group
 * definition), so grow and retrace generate exactly the same list of segments, as required for
 * detailed balance.
 */
struct GrowSegment
{
  std::optional<std::size_t> previousBead;
  std::size_t currentBead;
  std::vector<std::size_t> nextBeads;
  bool rigidUnit{false};
  std::optional<std::size_t> rigidGroup;  ///< Index into 'Component::groups' when 'rigidUnit'.
  Potentials::IntraMolecularPotentials intra;
};

/**
 * \brief Builds the deterministic sequence of growth segments starting from 'beadsAlreadyPlaced'.
 *
 * When the component has no rigid groups the sequence reproduces exactly the bead-by-bead order of
 * 'ConnectivityTable::nextBeads' (so fully flexible molecules are unaffected). When rigid groups are
 * present each rigid group is emitted as a single rigid-unit segment (RASPA2 'SetGrowingStatus').
 */
std::vector<GrowSegment> buildGrowSegments(const Component &component,
                                           const std::vector<std::size_t> &beadsAlreadyPlaced);
}  // namespace CBMC
