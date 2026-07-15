module;

export module phonon_kpath;

import std;

import double3;
import simulationbox;

/**
 * A labeled node along a phonon band-structure path, expressed in fractional (crystallographic)
 * reciprocal-lattice coordinates. Consecutive nodes define the straight segments of the path; the label
 * (e.g. "G", "X", "L") is used to annotate the high-symmetry point in the output.
 */
export struct PhononPathNode
{
  double3 kPoint{};
  std::string label{};
  /**
   * When true, the segment connecting the previous node to this node is a discontinuity and is not sampled:
   * this node instead starts a new (disconnected) sub-path at the same cumulative path coordinate. This is
   * used for the breaks that appear in crystallographic high-symmetry paths (e.g. "... U | K ..."). The first
   * node in the list ignores this flag.
   */
  bool startsNewSegment{false};
};

/**
 * A single sampled point of the band-structure path: its fractional coordinate, the cumulative distance
 * along the path measured in the Cartesian reciprocal metric (so segments are correctly scaled for
 * plotting), and an optional label inherited from the node when the point coincides with one.
 */
export struct PhononKPoint
{
  double3 kFractional{};
  double pathCoordinate{};
  std::string label{};
};

/**
 * Expand a sequence of high-symmetry nodes into a densely sampled band-structure path.
 *
 * Each segment `[nodes[i], nodes[i+1]]` is sampled with `pointsPerSegment` equal fractional steps
 * (endpoints inclusive; interior nodes are shared between neighbouring segments and appear once). The
 * cumulative path coordinate accumulates the Euclidean distance in Cartesian reciprocal space, computed
 * from the reciprocal-lattice basis of `box`, so that flat and steep regions are spaced physically. An
 * empty node list yields an empty path; a single node yields that node at coordinate zero.
 */
export std::vector<PhononKPoint> buildPhononKPath(std::span<const PhononPathNode> nodes,
                                                  std::size_t pointsPerSegment, const SimulationBox& box);
