module;

export module voronoi_surface_area;

import std;

import framework;
import forcefield;
import voronoi_accessibility;

// Monte-Carlo accessible surface area split into accessible and inaccessible parts,
// using the Voronoi accessibility classifier. Points are sampled on each atom's
// probe-inflated sphere; points overlapping another inflated atom are rejected, and the
// remaining points are classified accessible (channel) or inaccessible (pocket).
export struct SurfaceAreaSample
{
  double accessible{0.0};    // Å²
  double inaccessible{0.0};  // Å²
};

// The sampling itself, over whatever accessibility classifier is handed to it, so that the same
// estimate can be made of a network taken from the radical diagram or from the Apollonius diagram.
// `density` is the number of sample points per Å² of inflated sphere surface.
export SurfaceAreaSample sampleAccessibleSurfaceArea(const VoronoiAccessibility& accessibility, std::size_t density);

export struct VoronoiSurfaceArea
{
  double accessibleSurfaceArea{0.0};     // Å²
  double inaccessibleSurfaceArea{0.0};   // Å²

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> samplesPerAtom = std::nullopt);
};
