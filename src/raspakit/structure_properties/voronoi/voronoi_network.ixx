module;

export module voronoi_network;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;

// Voronoi network in the sense of zeo++: the graph of Voronoi vertices (nodes) and the
// polyhedron edges connecting them, annotated with the geometry needed for pore analysis.
//
// The network is built on top of the metric-aware fractional Voronoi construction
// (module skvoronoi). Each Voronoi vertex becomes a node whose radius `rad_stat_sphere`
// is the radius of the largest sphere centred on the vertex that does not overlap any
// atom (distance from the vertex to its nearest generating atom minus that atom's
// radius). Each polyhedron edge becomes a (periodic) network edge whose radius
// `rad_moving_sphere` is the bottleneck: the largest sphere that can travel along the
// edge without overlapping the atoms whose Voronoi cells share the edge.

export struct VoronoiNode
{
  double3 position;   // Cartesian position in the home unit cell [Å]
  double3 fractional; // fractional position wrapped into [0,1)
  double radius;      // rad_stat_sphere: largest included sphere radius at this node [Å]
  std::vector<std::size_t> atomIndices;  // indices of the atoms whose cells meet here
};

export struct VoronoiEdge
{
  std::size_t from;
  std::size_t to;
  int3 delta;      // lattice image shift added to `to` relative to `from`
  double radius;   // rad_moving_sphere: bottleneck radius along the edge [Å]
  double length;   // Cartesian edge length [Å]
};

export struct VoronoiNetwork
{
  std::vector<VoronoiNode> nodes;
  std::vector<VoronoiEdge> edges;  // directed; every undirected edge is stored both ways

  // Per atom, the nodes of its own Voronoi cell together with the Cartesian vector from
  // the atom centre to that vertex (used by the accessibility line-of-sight test).
  std::vector<std::vector<std::pair<std::size_t, double3>>> atomNodeVectors;

  SimulationBox simulationBox;
  std::vector<double3> atomPositionsFractional;  // wrapped, as used by the construction
  std::vector<double> atomRadii;

  // Diameter of the largest included sphere (Di) = 2 * max node radius.
  double largestIncludedSphereDiameter() const;

  // Builds the network for atoms at the given fractional positions with the given radii.
  static VoronoiNetwork create(const SimulationBox& simulationBox,
                               const std::vector<double3>& fractionalPositions,
                               const std::vector<double>& radii);
};
