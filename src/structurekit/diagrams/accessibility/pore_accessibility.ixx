module;

export module pore_accessibility;

import std;

import int3;
import double3;
import unit_cell;
import voronoi_network;
import channel_analysis;

// Classifies arbitrary sample points as solid / accessible-void / inaccessible-void, the
// machinery shared by the accessible surface-area, accessible-volume, pore-size and
// blocking-sphere analyses.
//
// Following zeo++, the atoms are inflated by the probe radius, a pore network is built on
// the inflated atoms, and its nodes are labelled as belonging to channels (accessible) or
// pockets (inaccessible). A sample point is assigned to the cell it falls in; if it lies
// inside that atom's inflated sphere it is solid, otherwise the nearest node of that cell
// that is in line of sight decides accessible vs inaccessible.
//
// Which diagram the network came from is not this classifier's business, and the name says
// so. `create` builds the radical diagram, which is zeo++'s own; `createFromNetwork` takes
// one already built, which is how the Apollonius analyses reach the same classifier. All
// that changes between them is the metric the cells are cut by, and the classifier is told
// which one so that it attributes a point to the cell that diagram would put it in.

export struct PointClassification
{
  bool inside{false};      // inside an inflated atom (solid)
  bool accessible{false};  // in an accessible channel (only meaningful when !inside)
  std::int32_t poreId{-1}; // channel/pocket id of the deciding node, or -1
  bool resample{false};    // no line-of-sight node found; caller should resample

  // The lattice translation carrying `point` into the frame the nodes of `poreId` are assembled in, i.e.
  // into the one lift of that pore rather than whichever periodic copy of it happens to lie next to the
  // point.  Zero when nothing was decided.
  //
  // A caller asking only whether a point is in a channel has no use for this.  It is there for the ones
  // that integrate something over the boundary of a bounded pore, where the pieces are found on atoms of
  // the home cell but belong to different lifts, and add up to a closed surface only once brought
  // together.  Meaningful for a pocket; for a channel the offset depends on the path taken to the node,
  // there being more than one, and it says nothing.
  int3 latticeOffset{0, 0, 0};
};

// One periodic image of one atom, as a neighbour of some query point: where its centre is relative to
// that point, how big it is, which atom of the home cell it is a copy of, and which copy.
//
// `image` is the lattice translation of that copy, so its centre is `atomPositions[index]` shifted by
// `cell * image`. It is what lets a piece of surface found on one atom be recognised as the same piece
// found on the other: the two see it in lifts that differ by exactly this translation.
export struct NeighbourImage
{
  double3 delta;             // from the query point to this image's centre
  double radius{0.0};        // inflated
  std::size_t index{0};      // atom of the home cell this is a copy of
  int3 image{0, 0, 0};       // which copy
};

// The samplers built on this machinery measure a fixed structure, so they are seeded from a constant:
// two runs on the same input have to give the same answer. Only a simulation wants a seed that varies.
export constexpr std::size_t samplingSeed = 1400uz;

// How many times a sampler redraws a point `classify` could not decide before giving up on it. A point
// has to be replaced rather than dropped, or the undecidable ones quietly cost the structure their
// share of its area or volume; the limit is there so that a region where every point is undecidable
// cannot spin.
export constexpr std::size_t resampleLimit = 8uz;

export struct PoreAccessibility
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
  UnitCell unitCell;

  // What the radii above were inflated by, so that the bare atom is recoverable from them. Only the
  // inflated radii bound the room left for the probe's centre, and that is what the network is cut on,
  // but the solvent excluded surface lies on the bare spheres and every caller that wants it would
  // otherwise have to be told the probe radius a second time and be trusted to name the same one.
  double probeRadius{0.0};  // Å

  // The bare radius of an atom, which is what the excluded surface is carried on. Never negative: a
  // radius smaller than the probe would mean the atom was inflated from nothing, which the callers
  // above do not do, and clamping is cheaper than a precondition every caller has to restate.
  double bareRadius(std::size_t atomIndex) const { return std::max(0.0, atomRadii[atomIndex] - probeRadius); }

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

  // Builds the accessibility structure on the radical diagram. `radii` are the bare atom radii; atoms are
  // inflated internally by `probeRadius`.
  static PoreAccessibility create(const UnitCell& unitCell,
                                  const std::vector<double3>& fractionalPositions, const std::vector<double>& radii,
                                  double probeRadius);

  // Builds it from a pore network already constructed on the inflated atoms, for a caller that takes
  // its network from somewhere other than the radical diagram. The network carries the cell, the atom
  // positions and the inflated radii, so it is all that is needed besides the metric its cells are
  // cut by.
  static PoreAccessibility createFromNetwork(VoronoiNetwork network, Metric metric, double probeRadius = 0.0);

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

  // Which lift of `poreId` lies nearest `point`: the lattice translation `t` for which some node of the
  // pore, at `nodes[i].position + cell * (nodeLatticeOffset[i] + t)`, comes closest to it.
  //
  // A pore is one connected set of nodes, but it stands for a periodic family of them, and a point next to
  // one member of that family is not next to the others. Which member matters to anything integrating over
  // a pore's boundary, and the node that `classify` happened to decide the point by need not belong to it:
  // the fallbacks there are chosen to get the pore right, that being all they were ever asked for.
  int3 nearestLift(const double3& point, std::int32_t poreId) const;

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

  // Every inflated atom, in every periodic image, whose centre lies within `reach` of `point`: the
  // vector from the point to that centre, and the inflated radius there.
  //
  // Each image is an entry of its own, including the images of the atom sitting on the point itself,
  // which in a cell smaller than twice an atom's radius cut into its sphere exactly as any other
  // neighbour does. The copy at zero distance is returned too, so a caller measuring the surface of
  // one atom has to drop that one and keep the rest.
  std::vector<std::pair<double3, double>> neighbourAtoms(const double3& point, double reach) const;

  // The same query, carrying which atom each image is of and which image of it. A caller that only
  // measures one sphere needs neither; one that has to recognise the same piece of surface from both of
  // the spheres it lies on needs both, the two descriptions of it being on different atoms in different
  // lifts and having nothing else in common.
  std::vector<NeighbourImage> neighbourAtomImages(const double3& point, double reach) const;
};
