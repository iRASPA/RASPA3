module;

export module brute_force_solvent_excluded;

import std;

import double3;
import brute_force_structure;
import brute_force_voxels;
import brute_force_surface_area;

// The solvent-excluded surface split into its three kinds of patch, by generating each kind and measuring
// what survives.
//
// The surface a molecule rolls over is not the surface of the atoms. Where the probe touches one atom it
// traces that atom's own sphere, and those pieces are convex; where it rolls along two it sweeps part of a
// torus, and those pieces are saddle-shaped; where it settles into a corner between three it rests, and the
// piece of its own sphere facing the corner is concave. The exact route finds all three by solving for the
// circles, arcs and vertices of the accessible surface and integrating each piece in closed form.
//
// Here each kind is generated from its own definition and sampled:
//
//   Convex, from each atom. A point at direction u on atom i belongs to the surface exactly when the probe
//   resting against the atom there, centred at distance r_i + r_p along u, is free of every other atom.
//   That is one clearance evaluation, and the atom contributes its 4 pi r_i squared times the fraction of
//   directions that pass.
//
//   Saddle, from each pair of atoms whose probes can touch both. The centres that do form a circle, and as
//   the probe goes round it the arc of its surface between the two contact points sweeps a torus. The area
//   element of a surface of revolution is known, so the sweep is integrated over the part of the circle
//   where the probe is free, and over the arc between the contacts. Where the circle is narrower than the
//   probe the torus turns itself inside out; the surface stops where it reaches the axis, and so does this.
//
//   Concave, from each corner where a free probe touches three atoms. The piece is the spherical triangle
//   on the probe's own sphere whose corners are the three contact directions, and it is sampled as such.
//
// None of that is quite enough on its own, because a patch generated correctly can still be buried: another
// position of the probe may reach through and cut part of it away, and the exact route handles that by
// construction. Here it is caught by a test that applies to all three kinds at once. A point belongs to the
// surface only if no probe can be placed with the point strictly inside it, and the centres a probe can
// take are exactly the accessible surface, which the surface-area check has already sampled densely. So
// every generated point is held against that cloud, and any point with a sampled centre nearer than the
// probe's radius is discarded as buried.
export struct BruteForceSolventExcluded
{
  double convexArea{0.0};   // Å², the atoms' own surface where the probe touches it
  double saddleArea{0.0};   // Å², swept by the probe rolling along a pair
  double concaveArea{0.0};  // Å², the probe's own surface resting in a corner

  double totalArea() const { return convexArea + saddleArea + concaveArea; }
  double reentrantArea() const { return saddleArea + concaveArea; }

  double convexAreaError{0.0};
  double saddleAreaError{0.0};
  double concaveAreaError{0.0};

  // How much of each kind was generated and then discarded as buried under a probe that could reach it.
  // Large on any of the three means the structure has the probe overlapping itself and that the sampling
  // below is working near the edge of what it can resolve.
  double buriedConvexArea{0.0};
  double buriedSaddleArea{0.0};
  double buriedConcaveArea{0.0};

  std::size_t numberOfPairs{0};    // pairs whose probes can touch both
  std::size_t numberOfCorners{0};  // free probe positions touching three atoms or more

  // How many corners rest on three atoms, how many on four, and so on. A probe resting on more than three
  // is a corner whose patch is bounded by more edges than a triangle has, and finding the same corner once
  // per triple of the atoms it touches is how a sweep over triples double-counts it.
  std::array<std::size_t, 12> cornersByContacts{};

  // The same split, by area rather than by count, which is what says whether a disagreement about the
  // concave total is about the wedged corners or about the ordinary ones.
  std::array<double, 12> concaveByContacts{};

  // Corners with another corner almost on top of them. Two probe positions a hair apart are one wedge that
  // the arithmetic resolved twice, and each would lay down the same patch again, so a count above zero is
  // the first thing to look at when the concave area comes out high.
  std::size_t crowdedCorners{0};
  double closestCorners{0.0};  // Å, between the two nearest distinct corners found

  double probeRadius{0.0};  // Å

  double seconds{0.0};

  // `structure` carries the atoms' own radii, not inflated. `inflated` carries the same atoms with the
  // probe added, whose union's boundary is where a probe's centre can be. `surface` must be a sample of
  // that boundary with its points kept, which is what the buried test is made of.
  static BruteForceSolventExcluded compute(const BruteForceStructure &structure,
                                           const BruteForceStructure &inflated,
                                           const BruteForceSurfaceArea &surface, double probeRadius,
                                           std::size_t samplesPerAtom, std::size_t creaseSteps = 256,
                                           std::size_t cornerSamples = 4096);
};
