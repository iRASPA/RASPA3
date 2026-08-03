module;

export module exact_boundary_components;

import std;

import int3;
import double3;
import pore_accessibility;

// The boundary of the union of the probe-inflated atoms, cut into its connected surfaces.
//
// The surface-area sweep measures that boundary arc by arc and asks the classifier about every arc
// separately, which is both the whole cost of the sweep and its one source of error: two arcs of the same
// piece of surface can come back labelled with different pores, and a pocket whose boundary has lost an arc
// to its neighbour is no longer closed. Neither can happen to a labelling that is applied to a connected
// surface as a whole.
//
// What is connected to what is settled here, exactly and without sampling. The boundary is cut along the
// circles in which the spheres cut one another, into patches: a patch is a connected piece of one sphere's
// exposed surface. Patches meet along those circles, so the whole boundary is recovered by joining each
// patch to the ones across its own edges, and the connected surfaces are the components of that graph.
//
// Two things make this exact. The edges are the arcs of the bounding circles that no third sphere covers,
// and the crossings of two circles that no third one covers are the only places where the exposed part of a
// circle can begin or end -- a crossing hidden inside a third sphere has that sphere's surface on both sides
// of it and interrupts nothing. And at an uncovered crossing of two circles the exposed region occupies
// exactly one of the four quadrants, so the two edges bounding it are determined and the walk around a
// patch never has a choice to make.
//
// Each patch also carries the lattice translation into its component's own frame, accumulated along the
// graph in the same way the percolation walk accumulates it over nodes. That is what the divergence theorem
// needs of an arc, and it is what decides whether a surface closes on itself or on a translate of itself:
// a component with a cycle whose translations do not cancel is the wall of a channel, and one without is
// the boundary of something bounded.

// Every sphere reaching the sphere of one atom, as the cap it cuts on it: `take(axis, cosineHalfAngle,
// neighbourIndex, neighbourImage)` for each of them, `axis` being the unit direction from this atom's centre
// towards that neighbour's and `cosineHalfAngle` the circle the two spheres meet in. False where a
// neighbour swallows this sphere whole, in which case the atom carries no exposed surface at all and what
// has been taken so far is to be thrown away rather than used.
//
// There is one of these because the decomposition and the sweep have to agree about which caps a sphere
// has and in what order they come. The sweep is allowed to take the decomposition's crossings rather than
// search for them again, and they are the same crossings only if these are the same caps; two spellings of
// the same law of cosines in two files, agreeing by inspection today, is not something that stays true.
export template <typename Take>
bool eachCapCoveringAtom(const PoreAccessibility& accessibility, std::size_t atomIndex, Take&& take)
{
  const double radius = accessibility.atomRadii[atomIndex];
  const double3 centre = accessibility.atomPositions[atomIndex];

  for (const NeighbourImage& neighbour :
       accessibility.neighbourAtomImages(centre, radius + accessibility.maximumAtomRadius))
  {
    const double distance = neighbour.delta.length();
    if (distance < 1.0e-12)
    {
      // This sphere itself, which the query returns along with the rest, or another sitting exactly on it:
      // only a strictly larger one covers anything, and then it covers all of it.
      if (neighbour.radius > radius) return false;
      continue;
    }

    // The circle of intersection sits at this polar angle from the neighbour's direction. At least one
    // means the neighbour reaches nothing of this sphere; at most minus one means it swallows it whole.
    const double cosineHalfAngle =
        (radius * radius + distance * distance - neighbour.radius * neighbour.radius) / (2.0 * radius * distance);
    if (cosineHalfAngle >= 1.0) continue;
    if (cosineHalfAngle <= -1.0) return false;

    take(neighbour.delta * (1.0 / distance), cosineHalfAngle, neighbour.index, neighbour.image);
  }
  return true;
}

// One of the circles in which a neighbouring sphere cuts the sphere of one atom, with the identity of that
// neighbour: the same circle lies on the neighbour's sphere too, and joining a patch to the one across it
// is the only step here that leaves one atom for another.
export struct SphereCircle
{
  double3 axis{0.0, 0.0, 0.0};  // unit, from this atom's centre toward the neighbour's
  double cosineHalfAngle{0.0};
  double sineHalfAngle{0.0};
  double halfAngle{0.0};

  std::size_t neighbourIndex{0};  // atom of the home cell the neighbour is a copy of
  int3 neighbourImage{0, 0, 0};   // which copy

  // An orthonormal pair perpendicular to the axis, in which a point of the circle is
  // `axis * cosineHalfAngle + (first * cos t + second * sin t) * sineHalfAngle`.
  double3 first{0.0, 0.0, 0.0};
  double3 second{0.0, 0.0, 0.0};

  // The angles, sorted, of the crossings on this circle that no third sphere covers. They cut the circle
  // into arcs: arc k runs from `cornerAngles[k]` to `cornerAngles[k+1]`, the last wrapping round, and a
  // circle with no crossings at all is one arc that closes on itself.
  std::vector<double> cornerAngles;

  // Per arc, the patch of this atom it bounds, or -1 where the arc is covered and bounds nothing.
  std::vector<std::int32_t> arcPatch;

  std::size_t numberOfArcs() const { return std::max<std::size_t>(1, cornerAngles.size()); }

  // The point of the circle at angle `t`, as a unit direction from this atom's centre.
  double3 direction(double t) const
  {
    return axis * cosineHalfAngle + (first * std::cos(t) + second * std::sin(t)) * sineHalfAngle;
  }

  // Where a direction known to lie on this circle sits along it.
  double angleOf(const double3& unitDirection) const
  {
    double angle = std::atan2(double3::dot(unitDirection, second), double3::dot(unitDirection, first));
    return (angle < 0.0) ? angle + 2.0 * std::numbers::pi : angle;
  }
};

// The exposed surface of one atom: the circles that bound it and the patches they cut it into.
export struct SphereBoundary
{
  std::vector<SphereCircle> circles;
  std::size_t numberOfPatches{0};

  // The crossings of two of the circles that no third sphere covers, as unit directions from the centre.
  //
  // These are the corners of the patches, and they are also exactly the latitudes at which the exposed
  // length of a circle of latitude stops being analytic, which is what the surface sweep has to break its
  // quadrature at. Finding them is cubic in the number of circles --- a pair at a time, each tested against
  // all the rest --- and the sweep would otherwise find the same ones over again from the same circles, so
  // they are kept here for it.
  std::vector<double3> crossings;

  // A direction inside each patch, a little way in from one of its edges. It is where the patch is asked
  // which pore it faces, and it is what a point of the patch is joined to when one is looked up.
  std::vector<double3> patchRepresentative;

  // Directions in the exposed region worth trying as way points between two points of it. A single arc of
  // a great circle joins most pairs, but not all: two points can be separated by a cap that the arc between
  // them runs into and a path round would miss, and in a symmetric arrangement they can come out exactly
  // antipodal, where no arc between them is even determined. Both are settled by going through further
  // points, and these are the ones tried.
  std::vector<double3> wayPoints;

  // Which way points are joined to which, by arcs among themselves: the way points fall into groups, and two
  // points of the region that reach the same group are joined through it. One leg of a path is rarely enough
  // and two are not always either --- a region that winds round the sphere needs as many legs as it has
  // turns --- so the groups are what make the test independent of how many that is.
  //
  // Both are found only under `LoopMerge::paths`, which is what wants them. Under the nesting rule they stay
  // empty and `connectedOnSphere` falls back to the single arc, which proves less but proves nothing false.
  std::vector<std::int32_t> wayPointGroup;

  // No exposed surface at all: some neighbour swallows this atom whole.
  bool buried{false};
};

export struct BoundaryComponents
{
  std::vector<SphereBoundary> atoms;

  // Per atom, per patch: which connected surface the patch belongs to, and the lattice translation
  // carrying the patch into that surface's own frame.
  std::vector<std::vector<std::int32_t>> componentOfPatch;
  std::vector<std::vector<int3>> offsetOfPatch;

  std::size_t numberOfComponents{0};

  // Per component: whether the surface closes on a translate of itself rather than on itself, which is what
  // the wall of a channel does and what the boundary of anything bounded cannot do. Read off the
  // translations as they are accumulated, a cycle that does not cancel being exactly this.
  std::vector<std::uint8_t> componentPercolates;

  // Per component: the translations it closes on itself by, and how many independent directions they span.
  //
  // Walking the patch graph and adding up the translations along the way gives each patch a lift of itself.
  // Where the walk returns to a patch it has already lifted and disagrees, the difference is a translation
  // taking the surface onto itself, and the differences generate the whole group of such translations --- the
  // surface is periodic under exactly that subgroup of the lattice and no more. Its rank is how many
  // independent directions the surface runs away in, which is the dimensionality of the pore it walls: a cage
  // gives zero, the wall of a one-dimensional channel closes on itself along the channel and nowhere else and
  // gives one, a wall between sheets gives two, and an intersecting system gives three.
  //
  // Read off the same walk that decides whether the surface is bounded, in integer arithmetic, without a pore
  // network and without a probe path being traced through anything.
  //
  // What it cannot do alone is join two surfaces that wall the same pore. A surface faces one pore, and every
  // translation of the surface is one of the pore, so its rank is never more than the pore's; where the pore
  // is walled by a single surface the two are equal and the rank is the pore's own. Where it is walled by
  // several --- a cluster of atoms standing in open void, two interpenetrating nets, a bundle of parallel rods
  // in a common void --- each surface still gives a rank the pore's cannot fall below, and which of them face
  // the same pore is the one question a pore network is for. It is the same residue the area and the volume
  // leave behind, and no framework of the atlas has it.
  std::vector<std::vector<int3>> componentTranslations;
  std::vector<int> componentDimensionality;

  // Per component, the atom and patch the surface is asked about, chosen to be the piece with the most room
  // in front of it.
  std::vector<std::pair<std::size_t, std::size_t>> componentRepresentative;

  // A few more of them, in the same order, for a surface whose most open patch still does not answer. Any
  // patch of a surface stands for the whole of it, so a second one is not a second opinion: it is the same
  // question asked where it may be easier to answer.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> componentCandidates;

  // Diagnostics. `unjoinedPatches` counts the patches of one atom that could not be shown to be joined to
  // another patch of the same atom nor across an edge, which over-cuts the surface without mislabelling it;
  // `looseEdges` counts edges whose far side could not be identified, which would.
  std::size_t numberOfPatches{0};
  std::size_t unjoinedPatches{0};
  std::size_t looseEdges{0};

  // What joining the loops of a sphere into its patches cost, and how many loops there were to join. It is
  // the one step of the decomposition that can be done in more than one way, so it is the one worth timing.
  double mergeSeconds{0.0};
  std::size_t loopsToMerge{0};

  // Which patch of `atomIndex` holds the point at `unitDirection` from its centre, or -1 if the point could
  // not be joined to any of them. The join is by the arc of a great circle, which lies in the exposed
  // region or does not, in closed form.
  std::int32_t patchOfDirection(std::size_t atomIndex, const double3& unitDirection) const;

  // Which patch a point already known to lie on one of the sphere's bounding circles belongs to. This is
  // the cheap case and the certain one: the circles are cut into arcs that each carry their patch, so the
  // answer is a search for the circle the point is on and then for the arc along it, with no path to find
  // and nothing to fail. Every gap the sweep measures ends on such a point.
  std::int32_t patchOfCirclePoint(std::size_t atomIndex, const double3& unitDirection) const;
};

// What one connected surface faces, asked once for the whole of it.
export struct ComponentVerdict
{
  bool decided{false};
  bool accessible{false};
  std::int32_t poreId{-1};

  // Whether a node's free ball was reached from the surface along a path that stays in the void, which
  // settles the question rather than inferring it. The point a surface is asked about is chosen where the
  // surface is most open, which is what makes this reachable at all: on the boundary itself no free ball
  // ever holds the point, so a per-arc classification can never have a proof and this usually does.
  bool proved{false};

  // How far out the walk had to go, and how many steps it took. Reported to be looked at rather than used.
  double walked{0.0};
  std::size_t steps{0};
};

// Asks each connected surface which pore it faces, at the most open point of it, by walking out along the
// normal in steps no longer than the room at each: every step then stays inside a free ball about the
// previous point, so the whole path is in the void and whatever node ball it reaches is provably on the
// same side of the surface. Falls back to the classifier where no ball is reached.
//
// `reference`, where it is given, is the network the pore is read from, the surfaces and the walk staying
// with `accessibility`. The two differ when a boundary taken at one probe radius is to be divided by what
// another probe can reach: a surface enclosing a pore that a large probe cannot get out of may still stand
// in a channel a small one moves freely along, and it is the smaller probe's network that says so. It has to
// be built from the same atoms at a probe radius no larger, since the walk is clear of the larger atoms only.
export std::vector<ComponentVerdict> boundaryComponentVerdicts(const PoreAccessibility& accessibility,
                                                               const BoundaryComponents& components,
                                                               const PoreAccessibility* reference = nullptr);

// How two loops of edges on one sphere are shown to bound the same patch. Walking the edges settles a patch
// bounded by a single loop; a patch with a hole in it has more than one, and those have to be recognised as
// belonging together or the patch is cut into as many pieces as it has loops.
export enum class LoopMerge : std::uint8_t
{
  // Look for a path from one loop to the other that stays in the exposed region: one arc of a great circle,
  // or several through the way points. Finding one is a proof; not finding one is not a proof of the
  // contrary, so this can over-cut. Kept to check the other against, the two being independent.
  paths,

  // Ask instead whether each loop lies inside the other, in the sense of Quan and Stamm: from a point of the
  // one, the nearest point of the other's own circles is reached without crossing any of them, so the two
  // are on the same side of that loop, and which side is read off from whether that nearest point is on the
  // loop or on a covered stretch of the same circle. This decides rather than searches, and it is the
  // default: over the zeolite atlas it cuts the spheres the same way as the search everywhere except where
  // it merges what the search could not, and it costs a small fraction of what the search costs.
  nesting,
};

// Cuts the boundary of the union of `accessibility`'s inflated atoms into its connected surfaces.
export BoundaryComponents boundaryComponents(const PoreAccessibility& accessibility,
                                             LoopMerge rule = LoopMerge::nesting);

// Whether the arc of the great circle from `from` to `to` stays outside every one of `circles`, which for
// two directions of one sphere's exposed region is a proof that they are on the same patch of it. Exported
// for the tests, the closed form for it being worth checking on its own.
export bool exposedGreatCircleArc(const std::vector<SphereCircle>& circles, const double3& from,
                                  const double3& to);

// Whether two directions of one sphere's exposed region can be shown to be on the same patch: by an arc of a
// great circle between them, or by a path of any number of arcs through the way points. A true answer is a
// proof; a false one is only a failure to find a path, which over-cuts the surface and is counted rather than
// assumed away.
export bool connectedOnSphere(const SphereBoundary& boundary, const double3& from, const double3& to);
