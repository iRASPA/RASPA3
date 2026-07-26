module;

export module voronoi_accessibility;

import std;

import int3;
import double3;
import simulationbox;
import voronoi_network;
import voronoi_channels;

// Classifies arbitrary sample points as solid / accessible-void / inaccessible-void, the
// machinery shared by the accessible surface-area, accessible-volume and blocking-sphere
// analyses.
//
// Following zeo++, the atoms are inflated by the probe radius, a Voronoi network is built
// on the inflated atoms, and its nodes are labelled as belonging to channels (accessible)
// or pockets (inaccessible). A sample point is assigned to the nearest atom; if it lies
// inside that atom's inflated sphere it is solid, otherwise the nearest Voronoi node of
// that atom's cell that is in line of sight decides accessible vs inaccessible.

export struct PointClassification
{
  bool inside{false};      // inside an inflated atom (solid)
  bool accessible{false};  // in an accessible channel (only meaningful when !inside)
  std::int32_t poreId{-1}; // channel/pocket id of the deciding node, or -1
  bool resample{false};    // no line-of-sight node found; caller should resample
};

// The samplers built on this machinery measure a fixed structure, so they are seeded from a constant:
// two runs on the same input have to give the same answer. Only a simulation wants a seed that varies.
export constexpr std::size_t samplingSeed = 1400uz;

// How many times a sampler redraws a point `classify` could not decide before giving up on it. A point
// has to be replaced rather than dropped, or the undecidable ones quietly cost the structure their
// share of its area or volume; the limit is there so that a region where every point is undecidable
// cannot spin.
export constexpr std::size_t resampleLimit = 8uz;

export struct VoronoiAccessibility
{
  // Which distance decides the cell a sample point falls in. It has to be the one the network's cells
  // are cut by, since the point is attributed to a cell before that cell's nodes are consulted: the
  // power distance for the radical diagram, the clearance for the Apollonius diagram.
  enum class Metric
  {
    Power,     // |x - c|^2 - r^2
    Clearance  // |x - c| - r
  };

  VoronoiNetwork network;
  ChannelAnalysis channels;
  Metric metric{Metric::Power};
  std::vector<std::int8_t> nodeAccessible;  // per node: 1 accessible, 0 inaccessible

  // Per node, the room left there for the probe's centre. Measured rather than read off the node,
  // because a radical diagram's node radius is taken against the atoms that generated the vertex,
  // which need not be the nearest ones. A point within this distance of the node is provably in the
  // same channel or pocket as the node, which is what `classify` falls back on.
  std::vector<double> nodeClearance;
  std::vector<double3> atomPositions;       // Cartesian, home cell
  std::vector<double> atomRadii;            // inflated radii used for the network
  SimulationBox simulationBox;

  // Cell list over the fractional unit cube, used to make the nearest-atom and overlap
  // queries O(1) per sample instead of O(number of atoms).
  int3 gridSize{1, 1, 1};
  std::vector<std::vector<std::size_t>> bins;
  double minimumBinWidth{0.0};
  double maximumAtomRadius{0.0};

  // The same cell list over the nodes, which makes the containment test below O(1) rather than a scan.
  // Each bin remembers the largest free ball it holds, so a bin too far away to reach the point with
  // even that ball is skipped without looking inside it, and the search stops once no bin left can
  // reach at all.
  std::vector<std::vector<std::size_t>> nodeBins;
  std::vector<double> nodeBinMaximumClearance;
  double maximumNodeClearance{0.0};

  // Builds the accessibility structure. `radii` are the bare atom radii; atoms are
  // inflated internally by `probeRadius`.
  static VoronoiAccessibility create(const SimulationBox& simulationBox,
                                     const std::vector<double3>& fractionalPositions,
                                     const std::vector<double>& radii, double probeRadius);

  // Builds it from a pore network already constructed on the inflated atoms, for a caller that takes
  // its network from somewhere other than the radical diagram. The network carries the cell, the atom
  // positions and the inflated radii, so it is all that is needed besides the metric its cells are
  // cut by.
  static VoronoiAccessibility createFromNetwork(VoronoiNetwork network, Metric metric);

  PointClassification classify(const double3& point) const;

  // The node whose free ball holds `point`, if there is one, taking the deepest where several do.
  //
  // This is the one statement about a point that is a proof rather than a reading of the diagram. The
  // ball of radius `nodeClearance` about a node holds no inflated atom, so every part of it is room the
  // probe's centre may occupy; being convex it is connected; so a point inside it belongs to the same
  // channel or pocket as the node, whatever the diagram's arcs happen to say. Where several balls hold
  // the point they are all one component -- overlapping free balls could not be otherwise -- so the
  // choice among them cannot change the answer.
  //
  // What it cannot do is decide a point no ball holds, which is any point closer to an atom than the
  // nearest node has room: on a surface-area sweep that is every sample. Those fall through to the
  // heuristics in `classify`.
  std::optional<std::size_t> containingNode(const double3& point, bool accessibleOnly = false) const;

  // Whether `point` is provably in a channel, by the above. Its use is to refuse a blocking sphere over
  // a point that has such a proof: a strongly degenerate diagram can leave a stretch of channel joined
  // to nothing and counted as a pocket, and blocking part of a pore is the one error here a simulation
  // cannot survive.
  bool provablyAccessible(const double3& point) const;

  // Room left for the probe's centre at `point`: min over atoms of |x - x_i| - r_i against the inflated
  // radii, so it is negative inside an atom and is the radius of the largest ball around the point that
  // holds no atom. That ball is connected and free, so it lies wholly within whichever channel or
  // pocket the point belongs to, which is what makes it a safe size for a blocking sphere.
  //
  // Not the same as the distance `classify` finds: that one is measured in the metric the cells are cut
  // by, which for a radical diagram is the power distance.
  double clearance(const double3& point) const;

  // True when the point lies strictly inside the inflated sphere of any atom other than
  // `excludedAtom` (pass a value >= number of atoms to test against all atoms).
  bool overlapsAtom(const double3& point, std::size_t excludedAtom) const;
};
