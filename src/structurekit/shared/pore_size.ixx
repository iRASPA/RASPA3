module;

export module grid_pore_size;

import std;

import uint3;
import double3;
import unit_cell;

// The pore-size distribution in the sense of Gelb and Gubbins, as arithmetic on a field rather than as a
// property of a framework.
//
// The size at a point of the void is the diameter of the largest sphere that fits in the void and covers
// that point, max{2 d(c) : |x - c| <= d(c)} over the centres c. It is a two-step measurement: first how much
// room there is at every point, which is the distance from the point to the boundary of the void, and then a
// sphere of that size drawn about every point, with every point a sphere reaches told how big it was and
// keeping the largest it hears.
//
// The second step is the same whatever drew the boundary, and the first is only different in how the
// boundary is found. A clearance field has already done it: the field is the distance, analytically and to
// full precision. A field that is not a distance, an energy above all, has to have the distance measured off
// it, and `distanceToIsosurface` is that measurement.
//
// Everything here reads its field as an openness, in the same sense as the percolation sweep and the
// connected components alongside it: larger means more open, and the void is what lies at or above the level
// given. A clearance is an openness already; an energy is one once negated.

// How far a sphere is allowed to reach past its own radius and still be said to cover a grid point. A grid
// point stands for its voxel, and a sphere that reaches into the voxel has covered what the point stands
// for. The size of it is set by the shell of void within a voxel of the boundary: the largest sphere
// covering such a point touches it from the inside, so its centre lies on the line from the boundary through
// the point, and the nearest grid point to where the centre belongs is off that line by half a voxel
// diagonal, which costs half a diagonal in the sphere's radius and another half in the distance to be
// covered. Without it that shell reports the size of nothing at all and the curve grows a spurious tail.
export double coveringSlack(const UnitCell &unitCell, uint3 gridSize);

// The distance from every grid point to the surface where `openness` falls to `lowestOpenness`, positive
// inside the void and negative outside it, so that the result can be read wherever a clearance field is.
//
// The surface is placed where it crosses the edges of the grid, between the pairs of neighbouring points
// that straddle it, by interpolating the field along the edge. That is where marching cubes puts its
// vertices, and it is what keeps the answer from being quantised to whole voxels: without it every distance
// would be a multiple of the grid spacing and the distribution would come out as a comb rather than a curve.
//
// Interpolating over anything longer than a single edge would be worse than not interpolating at all. An
// energy climbs as the twelfth power of separation, so a point in a wall exceeds a point in a pore by ten
// orders of magnitude rather than by any modest factor, and a straight line drawn between the two puts the
// crossing almost on top of the pore point. Every distance then comes out far too small, by an amount that
// depends on how deep the wall behind it happens to be, and the largest sphere the void holds stops even
// growing with the level.
//
// What is left is still an approximation, in a way the clearance field is not. The nearest point of the
// surface need not lie on a grid edge, so a distance measured this way can be too large by up to about a
// voxel, and unlike the clearance it does not improve by being asked more politely, only by refining the
// grid.
export std::vector<float> distanceToIsosurface(uint3 gridSize, const UnitCell &unitCell,
                                               std::span<const float> openness, double lowestOpenness);

// A sphere about every point of the void, each point keeping the largest sphere that reached it. `distance`
// is the room at each point, as a clearance field holds it or as `distanceToIsosurface` measured it, and
// points where it is negative take no part.
//
// None of the spheres can be dropped for sitting inside another. The field is a distance, so it gains
// exactly one along its own gradient and less in every other direction, and one ball holds another only when
// the two centres and the nearest boundary point are collinear, which on a grid never happens exactly.
export std::vector<float> poreRadiusField(uint3 gridSize, const UnitCell &unitCell,
                                          std::span<const float> distance, double slack);

export struct PoreSizeCurvePoint
{
  double diameter{0.0};

  // The share of the void lying in pores at least this wide, and the derivative of it, which is the
  // distribution itself.
  double cumulativeVoidFraction{0.0};
  double distribution{0.0};

  // The same two, over the part of the void that can be reached from outside the crystal.
  double cumulativeAccessibleFraction{0.0};
  double accessibleDistribution{0.0};
};

export struct PoreSizeCurve
{
  std::vector<PoreSizeCurvePoint> points;

  // The largest sphere that fits in the void anywhere, which is the same quantity as Di.
  double largestDiameter{0.0};

  // The mean size, over whatever the curve was weighted by.
  double meanDiameter{0.0};

  std::size_t numberOfVoidVoxels{0};
  std::size_t numberOfAccessibleVoxels{0};

  // What the curve was normalised by: a count of points when unweighted, and the sum of the weights when
  // not. Only the ratio matters to the curve, so a weight known up to a constant is enough.
  double totalWeight{0.0};
  double totalAccessibleWeight{0.0};

  double binWidth{0.0};
};

// Gathers grid points against the size found at each of them. Nothing is fitted: the cumulative curve is
// exact for the field it was measured on, and only its derivative depends on how finely the sizes are
// divided up. `accessible` marks the points reachable from outside, and may be empty when nothing is.
//
// `weight` says how much each point counts for, and may be empty, which is every point counting for one.
// That is the usual reading of a pore-size distribution: the share of the *space* lying in pores of a given
// width. On an energy landscape there is a better one available. Weighting each point by its Boltzmann
// factor gives the share of the *molecules* in pores of a given width, which is what a measurement would
// see, and it does not treat a point the molecule visits constantly and one it never reaches as the same
// thing merely because both lie below the contour.
//
// It also all but removes the arbitrariness of that contour. Where the level is drawn changes which points
// are counted, but the points it adds and drops are the ones near it, and those are the ones whose weights
// are smallest. A curve that shifts when the level moves by a few kT is a curve reporting the level rather
// than the framework.
export PoreSizeCurve binPoreSizes(std::span<const float> poreRadius, std::span<const float> distance,
                                  std::span<const std::uint8_t> accessible, std::span<const double> weight,
                                  std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins);
