module;

export module brute_force_surface_area;

import std;

import double3;
import brute_force_structure;
import brute_force_voxels;

// The surface of the union of the inflated atoms, by throwing points at each atom's sphere.
//
// The route being checked integrates each patch of the surface exactly: it works out where each atom's
// sphere is cut by its neighbours, and integrates the uncut part latitude by latitude with Gauss-Legendre
// quadrature. That is a closed-form answer whose error is the quadrature's, parts in 1e12, and it depends
// on getting every circle of intersection, every crossing of two circles and every arc between them right.
//
// Here there are no circles and no arcs. A direction is drawn at random on each atom, the point at that
// direction is tested against every other atom of every image, and the fraction that survive is what the
// atom contributes of its own 4 pi R squared. Nothing is solved and nothing can be mis-ordered; the price
// is that the answer has a standard error of order one over the root of the number of points, so this can
// confirm the exact area to three or four digits and never to twelve.
//
// The split between the surface a molecule can reach and the surface behind a sealed wall is made the same
// blunt way: step off the surface point along its outward normal and ask which piece of the flooded void
// the point lands in. That shares nothing with the pore classifier the exact route asks.
export struct BruteForceSurfaceArea
{
  double totalArea{0.0};         // Å², the whole of the union's boundary
  double accessibleArea{0.0};    // Å², the part facing a channel
  double inaccessibleArea{0.0};  // Å², the part facing a pocket
  double undecidedArea{0.0};     // Å², where stepping off the surface landed nowhere the flood had labelled

  // The standard error of the total, from the binomial count on each atom. The split has the same error
  // scale; it is reported once, for the total, because the three are not independent.
  double totalAreaError{0.0};  // Å²

  std::vector<double> areaOfAtom;  // Å², one per atom of the unit cell

  std::size_t samplesPerAtom{0};
  std::size_t numberOfExposedSamples{0};

  double seconds{0.0};

  // The exposed sample points, kept for the solvent-excluded check that follows: these are exactly the
  // places a probe's centre can be while touching the framework, which is what generates the surface a
  // molecule actually rolls over. Only filled when asked for, being large.
  std::vector<double3> exposedPoints;

  // `structure` must carry the inflated radii, atom plus probe, so that its union's boundary is the surface
  // being measured. `voxels` decides the accessible split and must have been built from the same structure.
  static BruteForceSurfaceArea compute(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                                       std::size_t samplesPerAtom, bool keepExposedPoints = false);
};
