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
 * defining the bend/torsion reference. There are three kinds of segments:
 *  - Flexible: 'nextBeads' are individual beads (a branch point can place several at once) that are
 *    sampled bead-by-bead with the bond/bend/torsion Monte-Carlo scheme.
 *  - Rigid unit: 'nextBeads' are all the (still unplaced) atoms of a rigid group and are placed as
 *    one rigid body. 'rigidGroup' holds the index of that group in 'Component::groups'.
 *  - Ring unit: 'nextBeads' are all the (still unplaced) atoms of a cyclic group and are grown with
 *    ring-closure CBMC (the ring's internal bond/bend/torsion conformation is sampled and the ring
 *    is placed with a coupled-decoupled orientational bias). 'cyclicGroup' holds the group index.
 *    'nextBeads' are ordered along the ring starting from the neighbor of the anchor 'currentBead',
 *    so consecutive beads are bonded and the last bead bonds back to the anchor (the closure bond).
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
  bool ringUnit{false};
  std::optional<std::size_t> rigidGroup;   ///< Index into 'Component::groups' when 'rigidUnit'.
  std::optional<std::size_t> cyclicGroup;  ///< Index into 'Component::groups' when 'ringUnit'.
  Potentials::IntraMolecularPotentials intra;
};

/**
 * \brief Builds the deterministic sequence of growth segments starting from 'beadsAlreadyPlaced'.
 *
 * When the component has no rigid or cyclic groups the sequence reproduces exactly the bead-by-bead
 * order of 'ConnectivityTable::nextBeads' (so fully flexible molecules are unaffected). When rigid
 * groups are present each rigid group is emitted as a single rigid-unit segment (RASPA2
 * 'SetGrowingStatus'); when cyclic groups are present each ring is emitted as a single ring-unit
 * segment grown with ring-closure CBMC.
 */
std::vector<GrowSegment> buildGrowSegments(const Component &component,
                                           const std::vector<std::size_t> &beadsAlreadyPlaced);
}  // namespace CBMC
