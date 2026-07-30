module;

export module exact_surface_patches;

import std;

import double3;
import pore_accessibility;
import exact_boundary_components;

// The order of the Gauss-Legendre rule used on each half of a smooth piece of the latitude integral, and the
// rest of the machinery for sweeping a sphere. Re-exported because the order is reported with the result.
export import exact_sphere_sweep;

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
// Where the region sits, as well as how large it is, comes from the same arcs and one more moment. The
// divergence theorem with the field |x - c|^2 e_j / 2, whose divergence is the jth component of x - c, turns
// the first moment of the region into a surface integral as well:
//
//   integral over the region of (x - c)  =  - sum_a [ (|d_a|^2 + R_i^2) m_a / 2 + R_i T_a d_a ] ,
//
// with d_a = p_i + L_a - c and T_a the integral of the outer product of the normal with itself over the arc.
// That tensor is as elementary in the latitude sweep as m_a is --- the azimuthal integrals of products of two
// sines and cosines are the same endpoints again --- and only its product with d_a is ever wanted, so what is
// accumulated is a vector and not a tensor. Dividing the sum by the enclosed volume gives the centroid.
export struct BoundaryMoments
{
  double area{0.0};                // Å²
  double radiusWeightedArea{0.0};  // Å³
  double originWeighted{0.0};      // Å³, sum of (p_i + L_a - c).m_a
  double3 vectorArea{0.0, 0.0, 0.0};  // Å², sum of m_a

  // The same patches carried down onto the bare spheres, which is the convex part of the solvent excluded
  // surface behind this one. A point of the accessible surface is the probe's centre touching atom i, and
  // the probe touches the atom along the same radius, so the excluded surface's convex patch is the radial
  // projection of the accessible one and its area is (r_i / R_i)^2 times as large --- exactly, patch by
  // patch, whatever shape the patch is. Accumulated here because the ratio is the atom's and the total is
  // not: a sum over patches cannot be rescaled after the fact.
  double convexArea{0.0};  // Å²

  // Å³, the volume between this share of the accessible surface and the convex patch below it, which is the
  // part of the shell R_s standing over these patches: the area of each of them times R_i/3 (1 - (r_i/R_i)³).
  // With it, the volume a probe of this radius opens up in the region this surface bounds is the volume the
  // surface encloses plus this, so a bounded pore's own share of the pore volume is available surface by
  // surface and not only as a total.
  double shellVolume{0.0};

  // Å⁴, the integral of x - c over the region the surface encloses, so the centroid is c plus this over the
  // enclosed volume. It carries the volume's sign with it, and so is meaningless for a surface that does not
  // close, exactly as the volume is.
  double3 enclosedFirstMoment{0.0, 0.0, 0.0};
};

// How much of a surface's moments a sweep is asked for.
//
// Everything below `enclosedFirstMoment` is wanted by every route: the area, the volume the surface encloses,
// the shell over it, and the vector area that says whether it closes at all. The first moment is wanted by one
// route only, the pocket centres of the void split, and it is much the most expensive part of an arc: the
// volume needs the integral of the normal over the arc, which is three sums, while the first moment needs the
// integral of the outer product of the normal with itself, which is six more and a tensor applied to a vector
// for every one of some hundreds of thousands of arcs. So it is asked for rather than always taken.
export enum class SurfaceMoments : std::uint8_t
{
  volume,     // area, the volume enclosed, the shell over it, and the closure defect
  andCentre,  // and the first moment, which is what a centroid comes from
};

export struct MeasuredPatches
{
  double accessible{0.0};    // Å²
  double inaccessible{0.0};  // Å²

  // Area the classifier could not place, which is area on a piece of surface that faces no node of the
  // network with room for the probe. It is reported rather than divided up, so that a structure where
  // the diagram is too degenerate to answer says so instead of quietly assigning the area to a side.
  double undecided{0.0};  // Å²

  // The whole area, and the same area with each atom's patch weighted by its own inflated radius. The
  // second is what the divergence theorem asks of the spheres when this decomposition is used for the
  // volume of the union rather than for its area, and it costs nothing to accumulate here.
  double area{0.0};                // Å²
  double radiusWeightedArea{0.0};  // Å³

  // Per atom, the area of its own exposed patch. The whole-surface total says nothing about how the area is
  // spread over the spheres, and a weighting that is neither the radius nor a constant --- the shell between
  // the accessible and the excluded surface wants r_i^3 / R_i^2 --- cannot be recovered from the total. It is
  // filled on every route, costing one addition per arc.
  //
  // Whole-structure totals of such weightings are taken from this rather than accumulated arc by arc
  // alongside it. There are some hundreds of atoms and some hundreds of thousands of arcs, so the sum over
  // atoms is both the shorter one and the more accurate, and an arc costs one addition instead of one per
  // weighting anyone might want.
  std::vector<double> atomArea;  // Å²

  // Per pore of the classifier, the same arcs collected pore by pore. Empty unless the arcs were
  // classified. Only the bounded pores are of any use: a channel's share is not a closed surface and its
  // moments do not mean anything.
  std::vector<BoundaryMoments> pores;

  // Per connected surface, the same arcs collected surface by surface. Filled only by the route that
  // classifies surfaces rather than arcs. This is the collection the divergence theorem wants: a surface
  // is closed or it is not, and which arcs belong to it is a matter of geometry rather than of a
  // classification, so a bounded one's moments close to round-off whatever the classifier makes of it.
  std::vector<BoundaryMoments> components;

  // What the sweep did, and where it could not answer. None of this enters any of the numbers above; it is
  // what a caller looks at to decide whether to believe them.
  struct Diagnostics
  {
    // Arcs measured, i.e. connected pieces of latitude, one classifier call each.
    std::size_t numberOfArcs{0};

    // Arcs whose patch could not be identified, and their area. A patch is looked up from a point of the
    // arc's own edge, which lies on one of the bounding circles and is found there in closed form, so this
    // counts the places where that failed rather than a tolerance being missed by a little.
    std::size_t unplacedArcs{0};
    double unplacedArea{0.0};  // Å²
  };
  Diagnostics diagnostics;

  // How the surfaces were divided, by the argument each of them was decided on. A surface that closes on a
  // translate of itself runs away through the crystal and has the void running away with it; one that closes
  // on itself around void is sealed; one that closes on itself around solid is a cluster of atoms, and only
  // for those is the network asked anything --- whether the void standing around the cluster is sealed in
  // turn. So `clusterSurfaces` is how much of this division is the network's and the rest is geometry.
  std::size_t numberOfSurfaces{0};
  std::size_t runawaySurfaces{0};
  std::size_t sealedSurfaces{0};
  std::size_t clusterSurfaces{0};

  // How much surface runs away in each number of independent directions, indexed by that number: a surface
  // walls a pore of its own dimensionality, so this is the area of the framework divided by the kind of pore
  // it faces. Index 0 is everything bounded, whether it seals void or stands round a cluster.
  //
  // The dimensionality of a surface is the rank of the group of lattice translations it closes on itself by,
  // which the decomposition accumulates on its way to deciding whether the surface is bounded at all. It is
  // integer arithmetic on the same walk, so it costs nothing and admits no tolerance, and it consults no pore
  // network: where the network route has to prune its edges at the probe radius and trust that connectivity
  // through a discrete graph is connectivity of the void, this reads the periodicity off the boundary itself.
  std::array<std::size_t, 4> surfacesOfDimension{};
  std::array<double, 4> areaOfDimension{};  // Å²

  // What was asked of the sweep, so that a caller reading `enclosedFirstMoment` can find out whether there is
  // one to read rather than reading a zero as an answer.
  SurfaceMoments moments{SurfaceMoments::volume};

  // The largest of those that carries any surface. It is the dimensionality of the pore system wherever each
  // pore is walled by a single surface, and a lower bound on it otherwise, a pore walled by several being one
  // only a network can put back together.
  int dimensionality() const
  {
    for (int rank = 3; rank > 0; --rank)
    {
      if (surfacesOfDimension[static_cast<std::size_t>(rank)] > 0) return rank;
    }
    return 0;
  }

  double total() const { return accessible + inaccessible + undecided; }
};

// The measured area of the boundary of the union of `accessibility`'s inflated atoms, undivided: the
// classifier is not consulted and `accessible`, `inaccessible` and `undecided` are all left at zero. That
// is what a caller after the area itself, or after the volume of the union, wants, and the classifier is
// the whole of the cost of the two routes that do consult it.
//
// `subdivisions` cuts each smooth piece of the latitude integral into that many panels of equal length. One
// settles the area to some thirteen digits, which is what the pieces being the right ones buys; raising it
// is a check on that rather than an improvement to it, and costs proportionally. A piece ending near a pole
// is cut finer there whatever is asked for, that being the one place where one panel is not enough.
export MeasuredPatches exactSurfaceArea(const PoreAccessibility& accessibility, std::size_t subdivisions = 1);

// The same, divided by asking the classifier about every arc, and collected pore by pore.
//
// This is the reference route and not the one to reach for. It asks the classifier about every arc of every
// panel --- some hundreds of thousands of questions where the route below asks some tens --- and each of
// those is an independent decision that can disagree with the others about the same piece of surface, which
// is what leaves a pocket's boundary open and its volume meaningless. Nothing outside the tests calls it. It
// is kept because it rests on nothing the other route rests on: the two arrive at the same area by different
// arguments, and a structure where they differ is a structure where one of them is wrong.
export MeasuredPatches exactAccessibleSurfaceAreaByPore(const PoreAccessibility& accessibility,
                                                        std::size_t subdivisions = 1);

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
//
// `moments` is what the surfaces are to be measured for; see `SurfaceMoments`. The default leaves out the one
// part no caller here needs, and a caller that does need it says so.
export MeasuredPatches exactAccessibleSurfaceAreaByComponent(const PoreAccessibility& accessibility,
                                                             const BoundaryComponents& components,
                                                             const std::vector<ComponentVerdict>& verdicts,
                                                             std::size_t subdivisions = 1,
                                                             SurfaceMoments moments = SurfaceMoments::volume);

// What settled which side of a connected surface the void that can be reached is on.
export enum class SurfaceClosure : std::uint8_t
{
  empty,    // no area: the surface is not there to be decided about
  runaway,  // closes on a translate of itself, so the void runs away along it
  sealed,   // closes on itself around void, which is then shut in behind it
  cluster,  // closes on itself around solid, and only the network can place the void outside it
};

export struct SurfaceSide
{
  SurfaceClosure closure{SurfaceClosure::empty};
  int side{0};  // plus one reachable, minus one shut in, zero where nothing could say
};

// Which side of each connected surface the reachable void is on.
//
// Three arguments decide it and only the last of them needs the network. A surface that closes on a translate
// of itself runs away through the crystal and the void runs away along it. A surface that closes on itself
// has the void on one side and the solid on the other, and the sign of the volume it encloses says which:
// positive is void, and void enclosed by a closed surface can reach nothing outside it. Only a surface
// enclosing solid --- a cluster of atoms with the void standing around it --- leaves a question, which is
// whether that void is itself sealed, and that is the one thing asked of `verdicts`.
//
// The area, the excluded volume and the void split all divide their totals by this, and they have to divide
// them the same way or the three disagree about the same structure. It is the whole of what the network is
// consulted for on this route, which is why it is one function and not a line repeated in each of them.
export std::vector<SurfaceSide> surfaceSides(const BoundaryComponents& components,
                                             const MeasuredPatches& patches,
                                             const std::vector<ComponentVerdict>& verdicts);

// The point each surface's moments are taken about: a point of the surface itself, carried into the frame the
// surface closes in. The choice cannot change a closed surface's volume, and it cancels out of the centroid
// too, but both are read back against it, so a caller that wants either has to ask for the same point the
// sweep used rather than assume one.
export std::vector<double3> surfaceMomentOrigins(const PoreAccessibility& accessibility,
                                                 const BoundaryComponents& components);

// The same, decomposing the boundary and judging it itself. For a caller that wants only the area, the
// decomposition costing less than the sweep it replaces the classification in.
export MeasuredPatches exactAccessibleSurfaceAreaByComponent(const PoreAccessibility& accessibility,
                                                             std::size_t subdivisions = 1,
                                                             SurfaceMoments moments = SurfaceMoments::volume);
