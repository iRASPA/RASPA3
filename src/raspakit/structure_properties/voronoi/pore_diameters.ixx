module;

export module voronoi_pore_diameters;

import std;

import framework;
import forcefield;
import voronoi_network;

// Largest included / free / included-along-free-path sphere diameters (Di, Df, Dif),
// the classic zeo++ ".res" quantities.
//
//   Di  : diameter of the largest sphere that fits anywhere in the pore space (static);
//         equals twice the largest Voronoi-node radius.
//   Df  : diameter of the largest sphere that can travel through the network along a
//         path that percolates through the periodic boundary (its bottleneck is the
//         smallest edge radius along the path).
//   Dif : diameter of the largest sphere that can be inscribed anywhere along that same
//         percolating (free-sphere) path.
export struct PoreDiameters
{
  double includedSphereDiameter{0.0};           // Di
  double freeSphereDiameter{0.0};               // Df
  double includedAlongFreePathDiameter{0.0};    // Dif

  static PoreDiameters compute(const VoronoiNetwork& network);
};

// The widest path that percolates through the periodic boundary, over the given nodes and the network
// edges between them: how wide its bottleneck is, which edge is that bottleneck, and the nodes the path
// can reach. Passing every node of the network answers the question Df asks; passing the nodes of one
// pore asks it of that pore alone.
//
// Widest is in the max-min sense, the path whose narrowest edge is as wide as possible, so the edges are
// added widest first and the first one that closes a loop with a non-zero net lattice offset is the
// bottleneck of the best percolating path there is.
export struct PercolatingPath
{
  bool percolates{false};
  std::size_t limitingEdge{0};  // index into network.edges of the bottleneck
  double radius{0.0};           // its radius [Å]; twice this is the Df of these nodes
  std::vector<std::size_t> componentNodes;
};

export PercolatingPath widestPercolatingPath(const VoronoiNetwork& network, std::span<const std::size_t> nodes);

export struct VoronoiPoreDiameters
{
  PoreDiameters result;
  VoronoiNetwork network;  // the network the diameters were read from, kept for the analyses that follow

  void run(const ForceField& forceField, const Framework& framework);
};
