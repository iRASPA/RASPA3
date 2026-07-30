module;

export module exact_solvent_excluded;

import std;

import double3;
import pore_accessibility;
import exact_boundary_components;
import exact_surface_patches;

// The solvent excluded surface of the framework, and with it the pore-size distribution in closed form.
//
// Gelb and Gubbins define the pore size at a point of the void as the diameter of the largest sphere that
// contains the point and fits in the void, and the distribution as the derivative of the volume where that
// diameter exceeds a given value. The set of such points is
//
//   V(r)  =  union over the balls of radius r that fit in the void  =  (E_r) dilated by B_r ,
//
// with E_r the room left for the probe's centre, and that is the morphological opening of the void. Dilation
// and erosion being dual, its complement is the erosion of the solid by the same ball:
//
//   V(r)  =  |cell|  -  |A_r eroded by B_r| ,     A_r = union of the atoms inflated by r .
//
// The eroded set is what a probe of radius r cannot enter, so its boundary is the solvent excluded surface of
// Richards and Connolly, and the whole distribution is the derivative of the volume that surface encloses. No
// point of the void is ever tested, and no trial sphere is ever drawn: the definition asks for the largest
// enclosing sphere at every point at once, and this is it.
//
// The volume comes from the shell decomposition of Quan and Stamm, which subtracts from the union of inflated
// atoms the shell between the accessible and the excluded surface,
//
//   |A_r eroded by B_r|  =  |A_r|  -  |R_s| ,      R_s = { p in A_r : the accessible surface is within r } ,
//
// rather than integrating over the excluded surface directly. The reason is that the excluded surface
// self-intersects where the probe is too large for the crevice it is rolling into, and the pieces of it then
// have to be clipped against one another, whereas the shell splits over the accessible surface's own patches,
// arcs and vertices into pieces that are disjoint because the nearest point of that surface is unique almost
// everywhere. One patch, one arc, one vertex, one closed-form piece each, and nothing counted twice.
//
// The distribution is the derivative of that volume in r, which by the moving-boundary argument is the
// integral of the normal speed of the excluded surface. The convex pieces of it lie on the spheres of the bare
// atoms, which do not move as the probe grows, so they contribute nothing however far the reentrant surface
// eats into them: their extent changes, but a boundary curve sliding along a fixed surface moves no volume.
// So
//
//   P(r) dr  =  integral over the reentrant surface of its normal speed ,
//
// a sum over the toroidal patches and the concave ones alone. That is the precise sense in which the
// pore-size distribution measures the concavity of the void, and it is why a structure of isolated convex
// atoms has P(r) = 0 for every r.

// One of the caps a region of a sphere is cut out by, as the direction of its axis and the cosine of its half
// angle. The cap is the part of the sphere it covers, and the region measured below is what no cap covers.
export struct SphericalCap
{
  double3 axis{0.0, 0.0, 1.0};
  double cosineHalfAngle{1.0};
};

export struct SphericalRegion
{
  double solidAngle{0.0};                // steradian
  double3 moment{0.0, 0.0, 0.0};         // the integral of the outward unit direction over the region
};

// The solid angle of the part of the unit sphere that no cap covers, and the integral of the direction over
// it. Both by the latitude sweep the surface area uses: the uncovered length of a circle of latitude is the
// complement of a union of arcs and is elementary, as are the two azimuthal moments, and the latitude
// integral is analytic between the latitudes where a cap turns back on itself or two caps cross uncovered.
//
// Nothing is assumed about the region. It may be empty, it may be the whole sphere, and it need not be
// connected or simply connected, which is why the sweep is used here and not the Gauss-Bonnet theorem: a
// concave patch clipped against its neighbours can come apart into several pieces, and then the theorem needs
// the loops traced and their number counted while the sweep needs nothing.
export SphericalRegion regionOutsideCaps(const std::vector<SphericalCap>& caps, std::size_t subdivisions = 1);

// The area of the excluded surface, by the kind of patch it is carried on.
//
// The surface has three kinds and no others, which is Richards's decomposition. Where the probe rests against
// one atom it lies on that atom's bare sphere, convex and curving away from the void. Where it rolls along the
// crease between two it traces a piece of a torus, curving away in one direction and towards the void in the
// other. Where it is wedged against three or more it stands still, and the surface is a piece of the probe's
// own sphere, concave and curving towards the void everywhere. So the split is the shape of the wall: convex
// area is atom seen bare, and reentrant area --- saddle and concave together --- is the part of the wall the
// probe cannot follow, which is the corrugation it bridges over.
//
// It is also exactly the split the pore-size distribution is carried by. Only the reentrant patches move as
// the probe grows, so a wall with no reentrant area has no distribution, and a wall that is nearly all
// reentrant is one whose every feature is finer than the probe.
export struct ExcludedSurfaceAreas
{
  double convex{0.0};   // Å², on the bare spheres of the atoms
  double saddle{0.0};   // Å², on the tori of the creases between two atoms
  double concave{0.0};  // Å², on the probe's own sphere where it is wedged against three or more

  double total() const { return convex + saddle + concave; }
  double reentrant() const { return saddle + concave; }

  // The three as fractions of the whole, which add to one wherever there is any area at all.
  double convexFraction() const { return (total() > 0.0) ? convex / total() : 0.0; }
  double saddleFraction() const { return (total() > 0.0) ? saddle / total() : 0.0; }
  double concaveFraction() const { return (total() > 0.0) ? concave / total() : 0.0; }
  double reentrantFraction() const { return (total() > 0.0) ? reentrant() / total() : 0.0; }
};

// The excluded volume of one probe size, and the distribution there.
export struct SolventExcludedGeometry
{
  double probeRadius{0.0};  // Å

  double accessibleVolume{0.0};  // Å³, |A_r|, the union of the inflated atoms
  double shellVolume{0.0};       // Å³, |R_s|, between the accessible surface and the excluded one
  double excludedVolume{0.0};    // Å³, |A_r| - |R_s|, what a probe of this radius cannot enter
  double poreVolume{0.0};        // Å³, the rest of the cell: the volume this probe opens up

  // dV_ses/dr = -dV/dr, in Å², which is the unnormalised distribution. Dividing by the void volume, the
  // pore volume at r = 0, normalises it to integrate to one over all r.
  double distribution{0.0};

  // The same, split by which pore the piece of reentrant surface it is on faces. Every toroidal and concave
  // patch belongs to one connected surface of the boundary, and that surface faces one pore, so the split is
  // geometry rather than a classification of samples.
  double accessibleDistribution{0.0};
  double inaccessibleDistribution{0.0};
  double undecidedDistribution{0.0};

  // The pore volume divided the same way. A surface that closes bounds a region of its own, and the volume
  // this probe opens up inside that region is the volume the surface encloses plus the shell standing over
  // it, both of which the same arcs give. The channels then take the rest, so nothing is measured twice and
  // the three add to `poreVolume`.
  double accessiblePoreVolume{0.0};
  double inaccessiblePoreVolume{0.0};
  double undecidedPoreVolume{0.0};

  // Per connected surface of the boundary: which side it faces, its share of the derivative, and the two parts
  // of the volume opened up in the region it bounds, which is what it encloses plus the shell over it. The
  // enclosed part means nothing for a surface that runs away through the crystal and is left at zero there.
  //
  // These are here for a caller that has to divide the curve by something other than what this probe can
  // reach --- the pore-size distribution divides it by a fixed probe while sweeping the diameter, so at every
  // diameter but one the two questions have different answers --- and such a caller needs the pieces rather
  // than the totals.
  std::vector<int> componentSide;
  std::vector<double> componentDistribution;    // Å²
  std::vector<double> componentEnclosedVolume;  // Å³, room for the probe's centre, zero if it does not close
  std::vector<double> componentShellVolume;     // Å³, the shell standing over it, on the pore's side of it

  // Where it comes from, the two families of reentrant patch.
  double toroidalDistribution{0.0};
  double concaveDistribution{0.0};

  // The area of the excluded surface, whole and by the side it faces. The sides are the connected surfaces of
  // the boundary again, so a patch is on the channel side or the pocket side because of where it sits and not
  // because a sample near it was classified: the convex patch belongs to the surface its accessible patch does,
  // and a crease or a wedge to the surface whose patches it joins.
  ExcludedSurfaceAreas area;
  ExcludedSurfaceAreas accessibleArea;
  ExcludedSurfaceAreas inaccessibleArea;
  ExcludedSurfaceAreas undecidedArea;

  // What the surface was made of, and which of the awkward cases the geometry turned out to hold. None of
  // this enters any of the numbers above; it is what a caller looks at to decide whether to believe them, and
  // it is kept apart from them so that the two are not mistaken for one another.
  struct Diagnostics
  {
    // `cuspedArcs` are the toroidal patches whose rolling circle is smaller than the probe, which are the ones
    // that fold back on themselves and have to be cut at the cusp; `clippedVertices` are the concave patches
    // that a neighbouring probe position reaches into.
    std::size_t numberOfArcs{0};
    std::size_t cuspedArcs{0};
    std::size_t numberOfVertices{0};
    std::size_t clippedVertices{0};

    // Vertices where more than three inflated spheres meet, which a framework's own symmetry produces rather
    // than a coincidence, and vertices whose concave patch a neighbouring probe clipped away entirely. Both
    // are handled and both are counted, a structure that is full of them being one to look at twice.
    std::size_t degenerateVertices{0};
    std::size_t vanishedVertices{0};

    // Corners of the decomposition that turned out not to be vertices of the surface at all. Counted apart
    // from the two above because it is a different thing being reported: those two are vertices the geometry
    // really has and this is the decomposition offering something the geometry does not.
    std::size_t discardedCorners{0};
  };
  Diagnostics diagnostics;
};

// The excluded volume and the distribution at one probe radius, against a boundary already decomposed and
// judged. `probeRadius` has to be the radius the accessibility's atoms were inflated by, the bare radii
// being read back from it, and it is the radius the accessibility itself records: building it from bare radii
// and a probe and passing the same probe here is the only consistent use.
export SolventExcludedGeometry solventExcludedGeometry(const PoreAccessibility& accessibility,
                                                       double probeRadius, const BoundaryComponents& components,
                                                       const std::vector<ComponentVerdict>& verdicts,
                                                       std::size_t subdivisions = 1);

// The same again, for a caller that has already swept the spheres. The sweep is the cost of the whole
// calculation, and a caller reporting the accessible surface area alongside the excluded one has done it
// already: `patches` has to be the one that sweep returned, against these components and verdicts.
export SolventExcludedGeometry solventExcludedGeometry(const PoreAccessibility& accessibility,
                                                       double probeRadius, const BoundaryComponents& components,
                                                       const std::vector<ComponentVerdict>& verdicts,
                                                       const MeasuredPatches& patches,
                                                       std::size_t subdivisions = 1);

// The same, decomposing the boundary and judging it itself.
export SolventExcludedGeometry solventExcludedGeometry(const PoreAccessibility& accessibility,
                                                       double probeRadius, std::size_t subdivisions = 1);
