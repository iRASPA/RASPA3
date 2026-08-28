#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import double3x3;
import surface_curvature;

// The curvature estimate, on surfaces whose principal curvatures can be written down.
//
// What is being checked is a fit and a sign convention, and the two fail differently. The fit is checked by
// putting a small triangle on a sphere, a cylinder and a saddle and asking for the curvatures back; the
// convention is checked by turning each of those inside out and asking that the classification follows. The
// last group checks the transform from a gradient held per grid step to a Cartesian normal, which is the part
// that is right by accident in a cubic cell and wrong in every other one.

namespace
{

// A triangle on the sphere of radius `radius` about the origin, small enough that the fit sees a sphere and
// large enough that the differences are not rounding. `outward` chooses whether the solid is inside the sphere
// or outside it, which is the whole of the difference between a bump and a pocket.
struct Patch
{
  std::array<double3, 3> corners;
  std::array<double3, 3> normals;
};

Patch spherePatch(double radius, double extent, bool outward)
{
  Patch patch;

  const std::array<double3, 3> directions{double3(0.0, 0.0, 1.0), double3(extent, 0.0, 1.0),
                                          double3(0.0, extent, 1.0)};

  for (std::size_t corner = 0; corner < 3; ++corner)
  {
    double3 unit = double3::normalize(directions[corner]);
    patch.corners[corner] = unit * radius;
    patch.normals[corner] = outward ? unit : -unit;
  }

  return patch;
}

// A triangle on the cylinder of radius `radius` about the z axis.
Patch cylinderPatch(double radius, double extent, bool outward)
{
  Patch patch;

  const std::array<double3, 3> angleAndHeight{double3(0.0, 0.0, 0.0), double3(extent, 0.0, 0.0),
                                              double3(0.5 * extent, radius * extent, 0.0)};

  for (std::size_t corner = 0; corner < 3; ++corner)
  {
    double angle = angleAndHeight[corner].x;
    double3 unit(std::cos(angle), std::sin(angle), 0.0);
    patch.corners[corner] = unit * radius + double3(0.0, 0.0, angleAndHeight[corner].y);
    patch.normals[corner] = outward ? unit : -unit;
  }

  return patch;
}

// A triangle on z = (x^2 - y^2) / (2 R) about the origin, whose principal curvatures at the origin are +1/R
// and -1/R with the normal taken upwards.
Patch saddlePatch(double radius, double extent)
{
  Patch patch;

  const std::array<std::array<double, 2>, 3> plane{
      std::array<double, 2>{-extent, -extent}, std::array<double, 2>{extent, -extent},
      std::array<double, 2>{0.0, extent}};

  for (std::size_t corner = 0; corner < 3; ++corner)
  {
    double x = plane[corner][0];
    double y = plane[corner][1];
    patch.corners[corner] = double3(x, y, (x * x - y * y) / (2.0 * radius));

    // The upward normal of z = f(x, y) is (-df/dx, -df/dy, 1) normalised.
    patch.normals[corner] = double3::normalize(double3(-x / radius, y / radius, 1.0));
  }

  return patch;
}

}  // namespace

// A sphere seen from outside is convex with both curvatures at one over its radius. The patch is put well away
// from the axes so that no accidental alignment of the tangent frame is being relied on.
TEST(grid_surface_curvature, a_sphere_from_outside_is_convex)
{
  for (double radius : {1.5, 3.0, 12.0})
  {
    Patch patch = spherePatch(radius, 0.02, true);
    TriangleCurvature curvature = triangleCurvature(patch.corners, patch.normals);

    ASSERT_TRUE(curvature.resolved) << "radius " << radius;
    EXPECT_NEAR(curvature.kappa1, 1.0 / radius, 1.0e-3 / radius) << "radius " << radius;
    EXPECT_NEAR(curvature.kappa2, 1.0 / radius, 1.0e-3 / radius) << "radius " << radius;
    EXPECT_NEAR(curvature.gaussian(), 1.0 / (radius * radius), 1.0e-3 / (radius * radius));

    EXPECT_EQ(classifyCurvature(curvature, CurvatureBands{}), CurvatureKind::Convex) << "radius " << radius;
  }
}

// The same sphere with the solid on the other side of it, which is a spherical pocket, and every curvature
// changes sign. This is the test that pins the convention: nothing about the triangle has moved.
TEST(grid_surface_curvature, a_sphere_from_inside_is_concave)
{
  for (double radius : {1.5, 3.0, 12.0})
  {
    Patch patch = spherePatch(radius, 0.02, false);
    TriangleCurvature curvature = triangleCurvature(patch.corners, patch.normals);

    ASSERT_TRUE(curvature.resolved) << "radius " << radius;
    EXPECT_NEAR(curvature.kappa1, -1.0 / radius, 1.0e-3 / radius) << "radius " << radius;
    EXPECT_NEAR(curvature.kappa2, -1.0 / radius, 1.0e-3 / radius) << "radius " << radius;

    // Gaussian curvature cannot tell the two apart, which is exactly why the sign of the mean is needed.
    EXPECT_GT(curvature.gaussian(), 0.0);
    EXPECT_LT(curvature.mean(), 0.0);

    EXPECT_EQ(classifyCurvature(curvature, CurvatureBands{}), CurvatureKind::Concave) << "radius " << radius;
  }
}

// A cylinder is curved one way and flat the other, so its Gaussian curvature is zero and the sign has to come
// from the curvature that is not. A rod counts with the convex and the wall of a cylindrical channel with the
// concave, which is what a chemist would say of either.
TEST(grid_surface_curvature, a_cylinder_takes_its_sign_from_the_one_curved_direction)
{
  const double radius = 4.0;

  Patch outside = cylinderPatch(radius, 0.01, true);
  TriangleCurvature rod = triangleCurvature(outside.corners, outside.normals);
  ASSERT_TRUE(rod.resolved);
  EXPECT_NEAR(rod.kappa1, 1.0 / radius, 1.0e-3 / radius);
  EXPECT_NEAR(rod.kappa2, 0.0, 1.0e-3 / radius);
  EXPECT_NEAR(rod.gaussian(), 0.0, 1.0e-3 / (radius * radius));
  EXPECT_EQ(classifyCurvature(rod, CurvatureBands{}), CurvatureKind::Convex);

  Patch inside = cylinderPatch(radius, 0.01, false);
  TriangleCurvature channel = triangleCurvature(inside.corners, inside.normals);
  ASSERT_TRUE(channel.resolved);
  EXPECT_NEAR(channel.kappa1, 0.0, 1.0e-3 / radius);
  EXPECT_NEAR(channel.kappa2, -1.0 / radius, 1.0e-3 / radius);
  EXPECT_EQ(classifyCurvature(channel, CurvatureBands{}), CurvatureKind::Concave);
}

// A saddle, where the two curvatures straddle zero. Its Gaussian curvature is negative, and that alone
// identifies it whichever way the normal points.
TEST(grid_surface_curvature, a_saddle_has_one_curvature_each_way)
{
  const double radius = 5.0;

  Patch patch = saddlePatch(radius, 0.01);
  TriangleCurvature curvature = triangleCurvature(patch.corners, patch.normals);

  ASSERT_TRUE(curvature.resolved);
  EXPECT_NEAR(curvature.kappa1, 1.0 / radius, 1.0e-3 / radius);
  EXPECT_NEAR(curvature.kappa2, -1.0 / radius, 1.0e-3 / radius);
  EXPECT_NEAR(curvature.mean(), 0.0, 1.0e-3 / radius);
  EXPECT_LT(curvature.gaussian(), 0.0);

  EXPECT_EQ(classifyCurvature(curvature, CurvatureBands{}), CurvatureKind::Saddle);
}

// A plane has no curvature at all and must not be given a sign. A large enough sphere is a plane as far as any
// cell here is concerned, and lands in the same column.
TEST(grid_surface_curvature, a_plane_is_flat_and_so_is_a_very_large_sphere)
{
  Patch flat;
  flat.corners = {double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0)};
  flat.normals = {double3(0.0, 0.0, 1.0), double3(0.0, 0.0, 1.0), double3(0.0, 0.0, 1.0)};

  TriangleCurvature curvature = triangleCurvature(flat.corners, flat.normals);
  ASSERT_TRUE(curvature.resolved);
  EXPECT_NEAR(curvature.kappa1, 0.0, 1.0e-12);
  EXPECT_NEAR(curvature.kappa2, 0.0, 1.0e-12);
  EXPECT_EQ(classifyCurvature(curvature, CurvatureBands{}), CurvatureKind::Flat);

  // The default band is one part in a hundred per Ångström, a radius of 100 Å.
  Patch huge = spherePatch(1000.0, 0.001, true);
  EXPECT_EQ(classifyCurvature(triangleCurvature(huge.corners, huge.normals), CurvatureBands{}),
            CurvatureKind::Flat);
}

// A triangle of no area, and one whose corners are collinear, carry area but no shape and have to say so
// rather than return a number.
TEST(grid_surface_curvature, a_degenerate_triangle_is_not_resolved)
{
  std::array<double3, 3> normals{double3(0.0, 0.0, 1.0), double3(0.0, 0.0, 1.0), double3(0.0, 0.0, 1.0)};

  std::array<double3, 3> coincident{double3(1.0, 2.0, 3.0), double3(1.0, 2.0, 3.0), double3(1.0, 2.0, 3.0)};
  EXPECT_FALSE(triangleCurvature(coincident, normals).resolved);

  std::array<double3, 3> collinear{double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(2.0, 0.0, 0.0)};
  EXPECT_FALSE(triangleCurvature(collinear, normals).resolved);
}

// A curvature sharper than the band is held back rather than counted, and a fit that failed goes to the same
// column but is counted separately so that a coarse grid can be told from a broken mesh.
TEST(grid_surface_curvature, the_bands_hold_back_what_the_grid_cannot_carry)
{
  CurvatureBands bands;
  bands.sharpest = 1.0;  // a radius of 1 Å

  Patch sharp = spherePatch(0.1, 0.02, true);
  TriangleCurvature curvature = triangleCurvature(sharp.corners, sharp.normals);
  ASSERT_TRUE(curvature.resolved);
  EXPECT_GT(curvature.kappa1, 1.0);
  EXPECT_EQ(classifyCurvature(curvature, bands), CurvatureKind::Unresolved);

  // With the bound off the same triangle is convex, so it is the bound doing the work and not the fit.
  EXPECT_EQ(classifyCurvature(curvature, CurvatureBands{}), CurvatureKind::Convex);

  Patch gentle = spherePatch(4.0, 0.02, true);

  CurvatureAreas areas;
  areas.add(2.0, curvature, bands);
  areas.add(3.0, triangleCurvature(gentle.corners, gentle.normals), bands);

  EXPECT_EQ(areas.numberOfUnresolved, 1u);
  EXPECT_EQ(areas.numberOfDegenerate, 0u);
  EXPECT_EQ(areas.numberOfConvex, 1u);
  EXPECT_NEAR(areas.unresolved, 2.0, 1.0e-12);
  EXPECT_NEAR(areas.classified(), 3.0, 1.0e-12);
  EXPECT_NEAR(areas.total(), 5.0, 1.0e-12);
  EXPECT_NEAR(areas.convexFraction(), 1.0, 1.0e-12);
  EXPECT_NEAR(areas.unresolvedFraction(), 0.4, 1.0e-12);

  // The integrals take both triangles, the band being about the columns only: 3 Å² of a sphere of radius 4 and
  // 2 Å² of one of radius 0.1.
  EXPECT_NEAR(areas.integratedMeanCurvature, 3.0 * 0.25 + 2.0 * 10.0, 1.0e-2);
  EXPECT_NEAR(areas.integratedGaussianCurvature, 3.0 * 0.0625 + 2.0 * 100.0, 1.0e-1);
}

// The four columns account for the classified area exactly, whatever went into them.
TEST(grid_surface_curvature, the_columns_add_up_to_the_whole)
{
  CurvatureBands bands;
  CurvatureAreas areas;

  Patch bump = spherePatch(3.0, 0.02, true);
  Patch pocket = spherePatch(3.0, 0.02, false);
  Patch saddle = saddlePatch(5.0, 0.01);

  areas.add(1.0, triangleCurvature(bump.corners, bump.normals), bands);
  areas.add(2.0, triangleCurvature(pocket.corners, pocket.normals), bands);
  areas.add(4.0, triangleCurvature(saddle.corners, saddle.normals), bands);

  EXPECT_NEAR(areas.convex, 1.0, 1.0e-12);
  EXPECT_NEAR(areas.concave, 2.0, 1.0e-12);
  EXPECT_NEAR(areas.saddle, 4.0, 1.0e-12);
  EXPECT_NEAR(areas.classified(), 7.0, 1.0e-12);
  EXPECT_NEAR(areas.convexFraction() + areas.saddleFraction() + areas.concaveFraction() + areas.flatFraction(), 1.0,
              1.0e-12);
}

// The gradient the extractors hold is per grid step, and carrying it to a position takes the inverse transpose
// of the cell. Dividing each component by its own spacing gives the same answer in a cubic cell and the wrong
// one in an oblique cell, so both are tested and the oblique one is the point of the test.
TEST(grid_surface_curvature, a_gradient_per_grid_step_becomes_a_cartesian_normal)
{
  const uint3 gridSize{64, 48, 80};

  auto check = [&](const double3x3 &cell, double3 wanted)
  {
    double3x3 inverseCell = cell.inverse();

    // A field whose value is the projection on `wanted`: its fractional gradient is cell^T wanted, and per
    // grid step that is divided by the number of steps along each axis.
    double3 fractionalGradient = transposedMultiply(cell, wanted);
    double3 perStep(fractionalGradient.x / static_cast<double>(gridSize.x),
                    fractionalGradient.y / static_cast<double>(gridSize.y),
                    fractionalGradient.z / static_cast<double>(gridSize.z));

    double3 normal = cartesianNormalOfGridGradient(inverseCell, gridSize, perStep, FieldSense::GrowsIntoVoid);
    double3 expected = double3::normalize(wanted);

    EXPECT_NEAR(normal.x, expected.x, 1.0e-10);
    EXPECT_NEAR(normal.y, expected.y, 1.0e-10);
    EXPECT_NEAR(normal.z, expected.z, 1.0e-10);

    // A field that grows into the wall instead gives the same surface with the normal the other way round.
    double3 turned = cartesianNormalOfGridGradient(inverseCell, gridSize, perStep, FieldSense::GrowsIntoSolid);
    EXPECT_NEAR(turned.x, -expected.x, 1.0e-10);
    EXPECT_NEAR(turned.y, -expected.y, 1.0e-10);
    EXPECT_NEAR(turned.z, -expected.z, 1.0e-10);
  };

  double3x3 cubic(20.0, 20.0, 20.0);
  check(cubic, double3(1.0, 0.0, 0.0));
  check(cubic, double3(0.3, -0.7, 0.5));

  // Monoclinic, so that a component of the fractional gradient mixes into more than one Cartesian direction.
  double3x3 oblique;
  oblique.ax = 18.0;
  oblique.ay = 0.0;
  oblique.az = 0.0;
  oblique.bx = 5.0;
  oblique.by = 14.0;
  oblique.bz = 0.0;
  oblique.cx = -3.0;
  oblique.cy = 2.0;
  oblique.cz = 22.0;

  check(oblique, double3(1.0, 0.0, 0.0));
  check(oblique, double3(0.0, 1.0, 0.0));
  check(oblique, double3(0.3, -0.7, 0.5));

  // The scale of the gradient is irrelevant, only its direction, so a normalised input and an unnormalised one
  // have to agree. Marching cubes normalises its stored gradient and the GPU kernel does not.
  double3 fractionalGradient = transposedMultiply(oblique, double3(0.3, -0.7, 0.5));
  double3 large =
      cartesianNormalOfGridGradient(oblique.inverse(), gridSize, fractionalGradient * 1000.0, FieldSense::GrowsIntoVoid);
  double3 small =
      cartesianNormalOfGridGradient(oblique.inverse(), gridSize, fractionalGradient * 0.001, FieldSense::GrowsIntoVoid);
  EXPECT_NEAR(large.x, small.x, 1.0e-10);
  EXPECT_NEAR(large.y, small.y, 1.0e-10);
  EXPECT_NEAR(large.z, small.z, 1.0e-10);
}


// The sense of the field is the whole of the difference between a bump and a pocket, and getting it wrong
// exchanges the two while leaving the saddles alone. A single atom is the case where it shows: its energy grows
// into the wall, so read with `GrowsIntoSolid` it is convex and read the other way round it is concave.
TEST(grid_surface_curvature, the_sense_of_the_field_exchanges_convex_and_concave)
{
  const uint3 gridSize{100, 100, 100};
  double3x3 cell(20.0, 20.0, 20.0);
  double3x3 inverseCell = cell.inverse();

  const double radius = 3.0;
  const double extent = 0.02;

  // A patch of the sphere of radius `radius`, with the gradient of a field that climbs towards the centre, held
  // per grid step. Pointing at the centre is what an energy field does: it grows into the atom.
  std::array<double3, 3> corners;
  std::array<double3, 3> intoSolid;
  const std::array<double3, 3> directions{double3(0.0, 0.0, 1.0), double3(extent, 0.0, 1.0),
                                          double3(0.0, extent, 1.0)};

  for (std::size_t corner = 0; corner < 3; ++corner)
  {
    double3 unit = double3::normalize(directions[corner]);
    corners[corner] = unit * radius;

    // The Cartesian gradient is -unit; carrying it back to per-grid-step is the inverse of what the function
    // under test does, so that the round trip is what is being checked and not a hand-built normal.
    double3 fractional = transposedMultiply(cell, -unit);
    intoSolid[corner] = double3(fractional.x / static_cast<double>(gridSize.x),
                                fractional.y / static_cast<double>(gridSize.y),
                                fractional.z / static_cast<double>(gridSize.z));
  }

  std::array<double3, 3> outward;
  std::array<double3, 3> inward;
  for (std::size_t corner = 0; corner < 3; ++corner)
  {
    outward[corner] =
        cartesianNormalOfGridGradient(inverseCell, gridSize, intoSolid[corner], FieldSense::GrowsIntoSolid);
    inward[corner] = cartesianNormalOfGridGradient(inverseCell, gridSize, intoSolid[corner], FieldSense::GrowsIntoVoid);
  }

  TriangleCurvature asEnergy = triangleCurvature(corners, outward);
  ASSERT_TRUE(asEnergy.resolved);
  EXPECT_NEAR(asEnergy.kappa1, 1.0 / radius, 1.0e-3 / radius);
  EXPECT_NEAR(asEnergy.kappa2, 1.0 / radius, 1.0e-3 / radius);
  EXPECT_EQ(classifyCurvature(asEnergy, CurvatureBands{}), CurvatureKind::Convex);

  TriangleCurvature theWrongWay = triangleCurvature(corners, inward);
  ASSERT_TRUE(theWrongWay.resolved);
  EXPECT_EQ(classifyCurvature(theWrongWay, CurvatureBands{}), CurvatureKind::Concave);
}


// The bands follow the coarsest voxel edge, so a finer grid believes a sharper curvature.
TEST(grid_surface_curvature, the_bands_follow_the_grid_spacing)
{
  double3x3 cell(20.0, 10.0, 40.0);

  CurvatureBands coarse = curvatureBandsForGrid(cell, uint3{100, 100, 100});
  CurvatureBands fine = curvatureBandsForGrid(cell, uint3{400, 400, 400});

  // The coarsest edge is 40/n, so the bound is n/40.
  EXPECT_NEAR(coarse.sharpest, 100.0 / 40.0, 1.0e-12);
  EXPECT_NEAR(fine.sharpest, 400.0 / 40.0, 1.0e-12);

  // An anisotropic grid is bounded by whichever axis is coarsest, which here is the third.
  CurvatureBands anisotropic = curvatureBandsForGrid(cell, uint3{200, 100, 100});
  EXPECT_NEAR(anisotropic.sharpest, 100.0 / 40.0, 1.0e-12);
}
