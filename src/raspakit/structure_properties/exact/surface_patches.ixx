module;

export module exact_surface_patches;

import std;

import double3;
import voronoi_accessibility;
import exact_boundary_components;

// The surface of the union of the probe-inflated atoms, measured patch by patch instead of by throwing
// points at it.
//
// The surface in question is the sheet the probe's centre traces when it is rolled over the framework.
// Nothing here belongs to one diagram or another: a patch is the part of one atom's sphere that lies
// outside every other inflated atom, which asks only for the atoms and a neighbour list. Every point of
// the union's boundary is on exactly one such patch, so the patches are the surface, with nothing missed
// and nothing counted twice.
//
// Each patch is bounded by the circles in which its neighbours cut the sphere, and on any circle of
// latitude each neighbour covers one arc, known in closed form. The uncovered length of the latitude is
// the complement of a union of arcs, and the patch area is that length integrated over latitude. The
// integrand is piecewise analytic, breaking only where a bounding circle turns back on itself and at the
// corners where two of them cross uncovered, so integrating between those breaks by Gauss-Legendre
// quadrature returns the area to round-off. No sample count enters the answer.
//
// A classifier is taken all the same, for two things it alone can say: which neighbours there are, and,
// if asked, which pore each patch faces.

// The order of the Gauss-Legendre rule used on each half of a smooth piece of the latitude integral. Ten
// nodes integrate a polynomial of degree nineteen exactly, and what the pieces and the endpoint
// substitution leave to integrate is analytic on each of them, so ten already puts the error at
// round-off. Raising it buys nothing and costs a classifier call per node. It is exported to be reported
// with the result rather than to be chosen by a caller.
export inline constexpr std::size_t exactQuadratureOrder = 10uz;

// What one pore's share of the boundary carries, beyond its area: enough to measure the volume the pore
// encloses, if it encloses one.
//
// A pore that does not percolate is bounded, and its boundary lies wholly on the spheres, so the
// divergence theorem turns its volume into an integral over exactly these arcs. With the field x - c
// about a fixed c, the normal component on the sphere of atom i is R_i on the patch and (p_i + L_a - c).n
// besides, so the two sums below are all that is needed:
//
//   3 |pore|  =  - ( radiusWeightedArea + originWeighted ) .
//
// `vectorArea` is the integral of the normal over the same arcs. Over a closed surface it vanishes, which
// is what makes the choice of c immaterial, so what it comes out as here is a check on the boundary
// having been assembled whole and in one frame rather than a quantity wanted for itself.
export struct PoreBoundaryMoments
{
  double area{0.0};                // Å²
  double radiusWeightedArea{0.0};  // Å³
  double originWeighted{0.0};      // Å³, sum of (p_i + L_a - c).m_a
  double3 vectorArea{0.0, 0.0, 0.0};  // Å², sum of m_a
};

export struct ExactSurfaceAreaSample
{
  double accessible{0.0};    // Å²
  double inaccessible{0.0};  // Å²

  // Area the classifier could not place, which is area on a piece of surface that faces no node of the
  // network with room for the probe. It is reported rather than divided up, so that a structure where
  // the diagram is too degenerate to answer says so instead of quietly assigning the area to a side.
  double undecided{0.0};  // Å²

  // Arcs measured, i.e. connected pieces of latitude, one classifier call each.
  std::size_t numberOfArcs{0};

  // The whole area, and the same area with each atom's patch weighted by its own inflated radius. The
  // second is what the divergence theorem asks of the spheres when this decomposition is used for the
  // volume of the union rather than for its area, and it costs nothing to accumulate here.
  double area{0.0};                // Å²
  double radiusWeightedArea{0.0};  // Å³

  // Per pore of the classifier, the same arcs collected pore by pore. Empty unless the arcs were
  // classified. Only the bounded pores are of any use: a channel's share is not a closed surface and its
  // moments do not mean anything.
  std::vector<PoreBoundaryMoments> pores;

  // Per connected surface, the same arcs collected surface by surface. Filled only by the route that
  // classifies surfaces rather than arcs. This is the collection the divergence theorem wants: a surface
  // is closed or it is not, and which arcs belong to it is a matter of geometry rather than of a
  // classification, so a bounded one's moments close to round-off whatever the classifier makes of it.
  std::vector<PoreBoundaryMoments> components;

  // Arcs whose patch could not be identified, and their area. A patch is looked up from a point of the
  // arc's own edge, which lies on one of the bounding circles and is found there in closed form, so this
  // counts the places where that failed rather than a tolerance being missed by a little.
  std::size_t unplacedArcs{0};
  double unplacedArea{0.0};  // Å²

  // How the surfaces were divided, by the argument each of them was decided on. A surface that closes on a
  // translate of itself runs away through the crystal and has the void running away with it; one that closes
  // on itself around void is sealed; one that closes on itself around solid is a cluster of atoms, and only
  // for those is the network asked anything --- whether the void standing around the cluster is sealed in
  // turn. So `clusterSurfaces` is how much of this division is the network's and the rest is geometry.
  std::size_t numberOfSurfaces{0};
  std::size_t runawaySurfaces{0};
  std::size_t sealedSurfaces{0};
  std::size_t clusterSurfaces{0};

  double total() const { return accessible + inaccessible + undecided; }
};

// The measured area of the boundary of the union of `accessibility`'s inflated atoms, split by that
// classifier.
//
// `subdivisions` cuts each smooth piece of the latitude integral into that many panels. One settles the
// area to some ten digits, which is what the pieces being the right ones buys; raising it is a check on
// that rather than an improvement to it, and costs proportionally, the classifier being asked about
// every arc of every panel.
//
// With `classifyArcs` false the classifier is not consulted at all and the split is left at zero, which
// is what a caller after the area itself, or after the volume, wants: the classifier is the whole cost.
export ExactSurfaceAreaSample exactAccessibleSurfaceArea(const VoronoiAccessibility& accessibility,
                                                         std::size_t subdivisions = 1,
                                                         bool classifyArcs = true);

// The same measurement, split by the connected surfaces of the boundary rather than arc by arc.
//
// Every arc lies on one patch of one sphere, and the patch is found from a point of the arc's own edge: a
// gap in a latitude begins and ends where a bounding circle crosses it, and the arcs of the circles already
// carry the patch they bound. So the lookup is a search along a circle, not a classification, and an arc
// cannot come out on a different side from the surface it is part of.
//
// What that buys is both the cost and the correctness. The classifier is consulted once per connected
// surface instead of once per arc -- some tens of times instead of some hundreds of thousands -- and the
// share of a bounded surface is then a closed surface by construction, which is what the divergence theorem
// needs and what the arc-by-arc route could only be checked for afterwards.
export ExactSurfaceAreaSample exactAccessibleSurfaceAreaByComponent(const VoronoiAccessibility& accessibility,
                                                                    const BoundaryComponents& components,
                                                                    const std::vector<ComponentLabel>& labels,
                                                                    std::size_t subdivisions = 1);

// The same, decomposing the boundary and labelling it itself. For a caller that wants only the area, the
// decomposition costing less than the sweep it replaces the classification in.
export ExactSurfaceAreaSample exactAccessibleSurfaceAreaByComponent(const VoronoiAccessibility& accessibility,
                                                                    std::size_t subdivisions = 1);
