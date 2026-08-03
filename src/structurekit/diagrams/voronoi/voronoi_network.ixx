module;

export module voronoi_network;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;

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
  // Radius of the largest empty sphere in this node's own pocket of free space. A radical
  // Voronoi vertex is not where the clearance peaks once the atom radii differ; the peak sits
  // at a nearby Apollonius vertex, so this is at least `radius` and is what Di is read from.
  double maximalRadius{0.0};
  double3 maximalPosition{}; // Cartesian centre of that sphere, the deepest point of the pocket
  std::vector<std::size_t> atomIndices;  // indices of the atoms whose cells meet here
};

// A sphere tangent to, and outside, a set of atoms: a vertex of the Apollonius (additively
// weighted Voronoi) diagram, where the clearance min_j(|x - x_j| - r_j) is locally maximal.
export struct ApolloniusSphere
{
  double3 centre;
  double radius;
};

// The spheres tangent to and outside all four given spheres. Four such constraints determine
// the centre and radius up to a quadratic, so there are at most two; fewer are returned for
// degenerate configurations or when no solution has a non-negative radius.
export std::vector<ApolloniusSphere> apolloniusTangentSpheres(const std::array<double3, 4>& centres,
                                                              const std::array<double, 4>& radii);

export struct VoronoiEdge
{
  std::size_t from;
  std::size_t to;
  int3 delta;      // lattice image shift added to `to` relative to `from`
  double radius;   // rad_moving_sphere: bottleneck radius along the edge [Å]
  double length;   // Cartesian edge length [Å]

  // Where along the edge the bottleneck sits and which way the passage runs there, in the frame in
  // which `from` is the node in the home cell. Together they cut the window a probe has to get
  // through: the plane through the position, perpendicular to the direction.
  //
  // Only the Apollonius builder fills them, and `hasBottleneckGeometry` says whether it did. Its
  // radius is the clearance at a sampled point of the arc, so it has a point to name; the radical
  // builder's radius is a bound taken over the whole segment at once, min over atoms of the distance
  // from that atom to the segment, and belongs to no point of it.
  double3 bottleneckPosition{};
  double3 bottleneckDirection{};
  bool hasBottleneckGeometry{false};
};

export struct VoronoiNetwork
{
  std::vector<VoronoiNode> nodes;
  std::vector<VoronoiEdge> edges;  // directed; every undirected edge is stored both ways

  // Per atom, the nodes of its own Voronoi cell together with the Cartesian vector from
  // the atom centre to that vertex (used by the accessibility line-of-sight test).
  std::vector<std::vector<std::pair<std::size_t, double3>>> atomNodeVectors;

  UnitCell unitCell;
  std::vector<double3> atomPositionsFractional;  // wrapped, as used by the construction
  std::vector<double> atomRadii;

  // Diameter of the largest included sphere (Di) = 2 * max node radius.
  double largestIncludedSphereDiameter() const;

  // Builds the network for atoms at the given fractional positions with the given radii.
  static VoronoiNetwork create(const UnitCell& unitCell,
                               const std::vector<double3>& fractionalPositions,
                               const std::vector<double>& radii);
};
