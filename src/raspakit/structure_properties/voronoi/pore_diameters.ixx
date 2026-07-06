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

export struct VoronoiPoreDiameters
{
  PoreDiameters result;

  void run(const ForceField& forceField, const Framework& framework);
};
