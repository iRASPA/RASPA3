module;

export module surface_curvature;

import std;

import uint3;
import double3;
import double3x3;

// Whether a triangle of a marching-cubes surface sits on a convex, a saddle or a concave piece of the wall.
//
// One triangle is flat, and its own normal says nothing about curvature: curvature is how the normal turns, so
// it takes at least the three normals of the three corners. Marching cubes already has them. It puts a vertex
// where the field crosses the level along a cube edge and stores the field's gradient there, interpolated
// along that same edge, and the gradient of a field is the normal of its level sets. So the mesh already
// carries the derivative a curvature needs, and the area loops have simply been throwing it away.
//
// The estimate is the usual one. On the plane of the triangle pick any orthonormal pair (u, v). The shape
// operator S is the 2 x 2 symmetric matrix carrying a step along the surface to the turn of the normal over
// that step,
//
//     n_j - n_i  =  S (p_j - p_i) ,
//
// and each of the three edges gives that twice, once against u and once against v. Six equations for the three
// entries of S, solved by least squares. The eigenvalues of S are the principal curvatures, and they do not
// depend on which pair (u, v) was picked, which is why no canonical frame has to be chosen.
//
// The sign convention is the whole content of the classification, so it is stated as a choice the caller makes
// rather than left implicit: `FieldSense` below says which way the field grows relative to the void, and the
// normal is turned to point away from the solid either way. With the normal that way round, a bump of solid ---
// a sphere of radius R seen from outside --- has both curvatures at +1/R, and a spherical pocket hollowed out of
// solid has both at -1/R. So positive is convex, negative is concave, and one of each is a saddle.
//
// What this is *not*: it is not Richards's convex/saddle/concave split of the solvent-excluded surface, which
// the exact route reports and which is a statement about how many atoms the probe is touching at once --- one,
// two, or three and more. That split lives on the excluded surface, where the saddle pieces are tori and the
// concave pieces are patches of the probe's own sphere. The two are different decompositions of two different
// surfaces and there is no reason for the numbers to agree.
//
// Whether the split converges as the grid is refined depends entirely on which field it is run on, and the two
// fields in this library fall on opposite sides of that.
//
// On the *clearance* field it does not converge, and the reason is worth spelling out because the columns look
// like a material property and are not one. The clearance field is a minimum over the atoms, and its level set
// is the boundary of a union of balls: every point of that not on a crease between two balls lies on a sphere
// and is convex, and the creases are curves carrying no area. So the honest limit is all convex and nothing
// else. On a grid it is not, because a crease is rounded over about a voxel and the band either side of it,
// where the central differences straddle the kink, comes back saddle or concave. Refining narrows that band and
// the convex share climbs towards one: on MFI from a third at 96^3 to two thirds at 320^3 without settling.
// There the columns are a way of localising creases at a stated spacing and nothing more.
//
// On an *energy* field it does converge, and the split is a real property of the surface. An energy field is a
// sum over the atoms rather than a minimum over them, so there is no crease anywhere: the field is smooth, the
// walls of neighbouring atoms blend into one another over a distance set by how fast the repulsion climbs, and
// the saddle between two atoms is a genuine smooth saddle with a curvature that does not depend on the grid.
// The three classes then mean what they say about the surface the probe sees.
//
// The two curvature integrals below converge in both cases and are the safest thing to quote.

// The two principal curvatures of one triangle, in 1/Å, largest first.
export struct TriangleCurvature
{
  double kappa1{0.0};
  double kappa2{0.0};

  // False when the fit had nothing to work with: a triangle of no area, or three normals so nearly equal that
  // the least-squares matrix is singular. Such a triangle carries area but no shape.
  bool resolved{false};

  double mean() const { return 0.5 * (this->kappa1 + this->kappa2); }
  double gaussian() const { return this->kappa1 * this->kappa2; }
  double sharpest() const { return std::max(std::abs(this->kappa1), std::abs(this->kappa2)); }
};

export enum class CurvatureKind
{
  Flat,
  Convex,
  Saddle,
  Concave,
  Unresolved
};

// The two curvatures that bound the range worth believing, both in 1/Å.
//
// Below `flat` a curvature is a radius larger than anything in the cell and is reported as flat rather than
// assigned a sign it does not really have. Above `sharpest` it is a radius smaller than a voxel, which the
// grid cannot represent, so it is not curvature that has been measured but the crease of the field: the
// clearance field has a crease everywhere two atoms meet, its gradient jumps across it, and the triangles
// straddling it produce enormous curvatures of either sign. Setting the bound to one over the coarsest voxel
// edge puts that debris in its own column instead of letting it flood the saddle one.
//
// Leaving `sharpest` at zero turns the upper bound off and classifies everything.
export struct CurvatureBands
{
  double flat{0.01};
  double sharpest{0.0};
};

export CurvatureKind classifyCurvature(const TriangleCurvature &curvature, const CurvatureBands &bands);

// Area by the shape of the wall it lies on, and how many triangles each column took.
export struct CurvatureAreas
{
  double convex{0.0};      // Å², both curvatures away from the void
  double saddle{0.0};      // Å², one each way
  double concave{0.0};     // Å², both curvatures towards the void
  double flat{0.0};        // Å², neither curvature large enough to have a sign
  double unresolved{0.0};  // Å², sharper than the grid can carry, or no shape at all

  std::size_t numberOfConvex{0};
  std::size_t numberOfSaddle{0};
  std::size_t numberOfConcave{0};
  std::size_t numberOfFlat{0};
  std::size_t numberOfUnresolved{0};

  // Of the unresolved ones, those where the fit itself failed rather than the curvature coming out too sharp.
  // A handful is ordinary; many means degenerate triangles, which is a fault of the extraction and not of the
  // grid.
  std::size_t numberOfDegenerate{0};

  // The integral of the mean curvature H over the whole surface, in Å, and of the Gaussian curvature K, which
  // is dimensionless.
  //
  // These are the part of this that survives refinement, and they are worth more than the columns above. The
  // columns are shares of area and the area near a crease shrinks with the spacing, so they move with the grid;
  // an integral over a crease does not, a crease carrying a curvature of one over the rounding radius over an
  // area proportional to that radius. So both integrals take every triangle whose fit stood up, whether the
  // band held it back from a column or not.
  //
  // The mean one is a Minkowski functional of the solid and can be had exactly for a union of balls, so it is
  // checkable, and on MFI it settles to within a percent by 256^3. The Gaussian one is 2 pi times the Euler
  // characteristic over a closed surface, so it is a topological count rather than a measurement: negative on a
  // network of channels, and the more negative the more connected the network. It is the noisier of the two by
  // a long way, because the fit here is a least-squares one over three normals and not the angle deficit that
  // satisfies the discrete Gauss-Bonnet theorem exactly. On a single sphere it lands within a part in a
  // thousand of 4 pi; on a framework it is worth a few percent at best, so read it for its sign and its size
  // and not past the first digit or two.
  double integratedMeanCurvature{0.0};
  double integratedGaussianCurvature{0.0};

  double classified() const { return this->convex + this->saddle + this->concave + this->flat; }
  double total() const { return this->classified() + this->unresolved; }

  // Fractions of the classified area, so that the first four add to one. The unresolved share is against the
  // whole instead, it being the thing the other four are not accounting for.
  double convexFraction() const { return (classified() > 0.0) ? this->convex / classified() : 0.0; }
  double saddleFraction() const { return (classified() > 0.0) ? this->saddle / classified() : 0.0; }
  double concaveFraction() const { return (classified() > 0.0) ? this->concave / classified() : 0.0; }
  double flatFraction() const { return (classified() > 0.0) ? this->flat / classified() : 0.0; }
  double unresolvedFraction() const { return (total() > 0.0) ? this->unresolved / total() : 0.0; }

  void add(double area, const TriangleCurvature &curvature, const CurvatureBands &bands);
};

// Which way the field grows relative to the void, which is what fixes the outward normal and so the sign of
// every curvature. A clearance field grows into the void, larger meaning more room, so its gradient already
// points away from the solid. An energy field grows into the wall, larger meaning deeper into the repulsion, so
// its gradient points at the solid and has to be turned round. Getting this wrong exchanges convex for concave
// and leaves the saddles where they are, which is why it is a named choice and not a sign buried in a caller.
//
// Every extractor here hands back the gradient in the field's own sense, pointing towards larger values,
// whatever the field is and whichever device found it. So this is the only place a sign is decided, and it is
// decided by what the field means rather than by which code path produced it.
export enum class FieldSense
{
  GrowsIntoVoid,
  GrowsIntoSolid
};

// The bands to use on a grid of this spacing. The upper one is one over the coarsest voxel edge, a curvature
// sharper than that being finer than the grid can carry whatever the field.
export CurvatureBands curvatureBandsForGrid(const double3x3 &unitCell, uint3 gridSize);

// Carries a gradient held per grid step to a gradient with respect to position, in the field's own sense.
//
// Everything here computes the gradient by differencing neighbouring samples, so what it holds is the
// derivative with respect to the grid index and not with respect to a position. A step along an index is 1/n
// of the cell along that axis, so the fractional gradient is n times as large; and a gradient is a covector,
// carried to a position by the inverse transpose of the cell rather than by the cell. Dividing each component
// by its own spacing instead is right only for a cell whose axes are orthogonal, and silently tilts every
// gradient on one whose axes are not.
//
// The magnitude is kept, so this is the one to use where the size of the gradient means something and not only
// its direction --- a directional derivative along a ray, for one.
export double3 cartesianGradientOfGridGradient(const double3x3 &inverseCell, uint3 gridSize, double3 gradient);

// The same, turned to point away from the solid and scaled to unit length: the outward normal of the level set
// through that point. The magnitude is dropped, the classification needing unit normals only, so it does not
// matter whether the caller's gradient was normalised beforehand or scaled by anything.
export double3 cartesianNormalOfGridGradient(const double3x3 &inverseCell, uint3 gridSize, double3 gradient,
                                             FieldSense sense);

// The principal curvatures of one triangle from its Cartesian corners and Cartesian unit normals, in the same
// order. The winding of the corners does not matter: the tangent frame is signed to agree with the normals.
export TriangleCurvature triangleCurvature(const std::array<double3, 3> &corners,
                                           const std::array<double3, 3> &normals);
