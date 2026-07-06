module;

export module voronoi_surface_area;

import std;

import framework;
import forcefield;

// Monte-Carlo accessible surface area split into accessible and inaccessible parts,
// using the Voronoi accessibility classifier. Points are sampled on each atom's
// probe-inflated sphere; points overlapping another inflated atom are rejected, and the
// remaining points are classified accessible (channel) or inaccessible (pocket).
export struct VoronoiSurfaceArea
{
  double accessibleSurfaceArea{0.0};     // Å²
  double inaccessibleSurfaceArea{0.0};   // Å²

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> samplesPerAtom = std::nullopt);
};
