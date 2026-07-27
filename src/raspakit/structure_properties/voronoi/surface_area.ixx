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

// Accessible and inaccessible surface area against the radical network, either measured or sampled.
//
// The area itself is a property of the union of the inflated atoms and not of any diagram, so the
// measured route here is the same one the Apollonius analysis uses and returns the same total to the last
// bit. What the network decides is the split: the radical network places the percolation threshold too
// low and reaches further into pockets than a probe can, so its division of the area is the approximate
// one. Running the same measurement against both networks therefore separates the error of the
// classifier from the error of the geometry, which two sampled numbers cannot do.
export struct VoronoiSurfaceArea
{
  enum class Method
  {
    Exact,   // the latitude integral, accurate to round-off
    Sampled  // points on the inflated spheres, as zeo++
  };

  double accessibleSurfaceArea{0.0};     // Å²
  double inaccessibleSurfaceArea{0.0};   // Å²
  double undecidedSurfaceArea{0.0};      // Å², only ever nonzero for the exact method

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           Method method = Method::Exact, std::optional<std::size_t> samplesPerAtom = std::nullopt,
           std::optional<std::size_t> subdivisions = std::nullopt);
};
