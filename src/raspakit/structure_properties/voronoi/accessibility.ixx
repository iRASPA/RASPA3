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

export struct VoronoiAccessibility
{
  VoronoiNetwork network;
  ChannelAnalysis channels;
  std::vector<std::int8_t> nodeAccessible;  // per node: 1 accessible, 0 inaccessible
  std::vector<double3> atomPositions;       // Cartesian, home cell
  std::vector<double> atomRadii;            // inflated radii used for the network
  SimulationBox simulationBox;

  // Cell list over the fractional unit cube, used to make the nearest-atom and overlap
  // queries O(1) per sample instead of O(number of atoms).
  int3 gridSize{1, 1, 1};
  std::vector<std::vector<std::size_t>> bins;
  double minimumBinWidth{0.0};
  double maximumAtomRadius{0.0};

  // Builds the accessibility structure. `radii` are the bare atom radii; atoms are
  // inflated internally by `probeRadius`.
  static VoronoiAccessibility create(const SimulationBox& simulationBox,
                                     const std::vector<double3>& fractionalPositions,
                                     const std::vector<double>& radii, double probeRadius);

  PointClassification classify(const double3& point) const;

  // True when the point lies strictly inside the inflated sphere of any atom other than
  // `excludedAtom` (pass a value >= number of atoms to test against all atoms).
  bool overlapsAtom(const double3& point, std::size_t excludedAtom) const;
};
