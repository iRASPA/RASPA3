module;

export module exact_solvent_excluded;

import std;

import double3;
import voronoi_accessibility;
import exact_boundary_components;

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

  // The same, split by which pore the piece of reentrant surface it is carried on faces. Every toroidal and
  // concave patch belongs to one connected surface of the boundary, and that surface faces one pore, so the
  // split is geometry rather than a classification of samples.
  double accessibleDistribution{0.0};
  double inaccessibleDistribution{0.0};
  double undecidedDistribution{0.0};

  // Where it comes from, the two families of reentrant patch.
  double toroidalDistribution{0.0};
  double concaveDistribution{0.0};

  // What the boundary was made of. `cuspedArcs` are the toroidal patches whose rolling circle is smaller than
  // the probe, which are the ones that fold back on themselves and have to be cut at the cusp;
  // `clippedVertices` are the concave patches that a neighbouring probe position reaches into.
  std::size_t numberOfArcs{0};
  std::size_t cuspedArcs{0};
  std::size_t numberOfVertices{0};
  std::size_t clippedVertices{0};

  // Vertices where more than three inflated spheres meet, which a framework's own symmetry produces rather
  // than a coincidence, and vertices whose concave patch was clipped away entirely. Both are handled and both
  // are counted, a structure that is full of them being one to look at twice.
  std::size_t degenerateVertices{0};
  std::size_t vanishedVertices{0};
};

// The excluded volume and the distribution at one probe radius, against a boundary already decomposed and
// labelled. `probeRadius` has to be the radius the accessibility's atoms were inflated by, the bare radii
// being read back from it.
export SolventExcludedGeometry solventExcludedGeometry(const VoronoiAccessibility& accessibility,
                                                       double probeRadius, const BoundaryComponents& components,
                                                       const std::vector<ComponentLabel>& labels,
                                                       std::size_t subdivisions = 1);

// The same, decomposing and labelling the boundary itself.
export SolventExcludedGeometry solventExcludedGeometry(const VoronoiAccessibility& accessibility,
                                                       double probeRadius, std::size_t subdivisions = 1);
