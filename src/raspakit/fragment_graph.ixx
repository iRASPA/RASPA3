module;

export module fragment_graph;

import std;

import archive;
import double3;
import atom;
import molecule;
import connectivity_table;
import fragment;

/**
 * \brief The fragment decomposition of a molecule plus the deterministic growth structure over it.
 *
 * A molecule is partitioned into fragments (see 'Fragment'): rigid bodies and single-atom flexible
 * beads. Bonds whose two atoms lie in different fragments form the edges of the fragment graph.
 * The graph is turned into
 *  - a spanning tree rooted at the fragment containing the molecule's starting bead ('parentFragment'
 *    / 'parentBond' / 'growthOrder'), which the CBMC growth follows so that grow and retrace produce
 *    the same sequence (required for detailed balance), and
 *  - a set of closure bonds ('closureBonds'): the non-tree edges. Each closure bond is a ring-closure
 *    constraint. A simple ring produces one closure bond; fused and bridged polycyclic systems
 *    produce several. There is no restriction to simple cycles.
 *
 * The runtime rigid-body dynamics state of the multi-atom fragments lives outside this structure, in
 * 'System::groupData' as one 'GroupState' per rigid fragment (in ascending fragment index order).
 */
export struct FragmentGraph
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::vector<Fragment> fragments{};        ///< Fragments in deterministic order (by lowest atom index).
  std::vector<std::size_t> atomFragmentIds{};  ///< Maps each atom to the index of its fragment.

  std::size_t rootFragment{0};  ///< Fragment containing the molecule's starting bead (spanning-tree root).

  /// Spanning-tree parent of every fragment (nullopt for the root and any disconnected component root).
  std::vector<std::optional<std::size_t>> parentFragment{};
  /// The bond {parentAtom, childAtom} connecting a fragment to its tree parent (nullopt for a root).
  std::vector<std::optional<std::array<std::size_t, 2>>> parentBond{};
  /// Fragments in spanning-tree visitation order (root first); the CBMC growth order.
  std::vector<std::size_t> growthOrder{};

  /// Non-tree inter-fragment bonds: ring-closure constraints, each stored as {atomA, atomB}.
  std::vector<std::array<std::size_t, 2>> closureBonds{};

  /// Cyclic clusters: the 2-edge-connected components of the fragment graph with more than one
  /// fragment. Each cluster is the sorted set of atoms of its fragments. A simple flexible ring is
  /// one cluster with one closure bond; fused and bridged polycyclic systems are one cluster with
  /// several closure bonds. Grown with ring-closure CBMC as a single unit.
  std::vector<std::vector<std::size_t>> cyclicClusters{};
  /// Maps each fragment to its cyclic cluster (-1 when the fragment is not part of any cycle).
  std::vector<std::make_signed_t<std::size_t>> fragmentCyclicClusterIds{};

  bool semiFlexible{false};      ///< Cached: at least one rigid fragment and more than one fragment.
  std::size_t rigidFragmentCount{0};  ///< Cached: number of multi-atom (rigid-body) fragments.

  FragmentGraph() = default;

  /**
   * \brief Deterministically partitions the atoms of a molecule into fragments.
   *
   * Atoms listed together in one entry of 'rigidBodies' become one rigid fragment; every atom not
   * listed becomes its own single-atom (flexible) fragment. Fragments are ordered by their lowest
   * atom index. Throws when a rigid body references an out-of-range atom or when an atom appears in
   * more than one rigid body.
   */
  static std::vector<std::vector<std::size_t>> partitionAtoms(
      std::size_t numberOfBeads, const std::vector<std::vector<std::size_t>> &rigidBodies);

  /**
   * \brief Builds the fragment graph from a connectivity table and an atom partition.
   *
   * Computes each fragment's rigid-body reference data (from 'referencePositions' / 'masses'), the
   * spanning tree rooted at the fragment containing 'startingBead', the growth order, and the closure
   * bonds. 'partition' must cover every atom exactly once (as produced by 'partitionAtoms').
   */
  void build(const ConnectivityTable &connectivity, const std::vector<std::vector<std::size_t>> &partition,
             std::size_t startingBead, std::span<const double3> referencePositions,
             std::span<const double> masses);

  /// Refreshes the cached summary ('semiFlexible', 'rigidFragmentCount'). Call after build/deserialize.
  void finalizeCache();

  /// Index of the fragment that 'atom' belongs to.
  std::size_t fragmentContaining(std::size_t atom) const { return atomFragmentIds[atom]; }

  /// Number of multi-atom (rigid-body) fragments.
  std::size_t numberOfRigidFragments() const { return rigidFragmentCount; }

  /**
   * \brief Whether the molecule is semi-flexible: it has at least one rigid body and more than one
   *        fragment (so it is neither fully rigid nor fully flexible).
   */
  bool isSemiFlexible() const { return semiFlexible; }

  /// Whether all given atom ids lie inside one and the same rigid-body fragment.
  bool isInsideRigidFragment(std::span<const std::size_t> ids) const;

  /**
   * \brief Returns the atoms of the next fragment the deterministic CBMC growth would place given
   *        'placedBeads', or std::nullopt when growth cannot proceed.
   *
   * The growth is fragment-at-a-time: it finds the first frontier bond (a placed bead connected to an
   * unplaced bead, scanning placed beads in order and neighbors by ascending index), and grows the
   * whole fragment containing that unplaced bead. Returns std::nullopt when there is no frontier bond
   * (nothing left to grow, or a disconnected remainder) or when the target fragment is only partially
   * placed (a fragment must be placed as one unit, e.g. a rigid body cannot be split by partial
   * reinsertion). Grow and retrace call this with the same argument, so they produce the same order.
   */
  std::optional<std::vector<std::size_t>> nextGrowBeads(const std::vector<std::size_t> &placedBeads,
                                                        const ConnectivityTable &connectivity) const;

  /// Indices (into 'fragments') of the rigid-body fragments, in ascending order.
  std::vector<std::size_t> rigidFragmentIndices() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const FragmentGraph &g);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, FragmentGraph &g);
};
