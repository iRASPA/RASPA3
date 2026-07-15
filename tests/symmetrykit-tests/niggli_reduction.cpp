#include <gtest/gtest.h>

import std;

import double3;
import double3x3;

import sksymmetrycell;
import sktransformationmatrix;

namespace
{
constexpr double radiansToDegrees = 180.0 / std::numbers::pi;

// An integer, determinant +1 change of basis. Re-expressing a lattice through it yields a different (skewed)
// presentation of the *same* lattice, whose Niggli-reduced cell must therefore be identical.
double3x3 unimodularSkew()
{
  return double3x3(double3(1.0, 0.0, 0.0), double3(1.0, 1.0, 0.0), double3(0.0, 1.0, 1.0));
}

void expectReducedCell(SKSymmetryCell& input, double a, double b, double c, double alphaDegrees,
                       double betaDegrees, double gammaDegrees, double tolerance = 1.0e-6)
{
  auto result = input.computeReducedNiggliCellAndChangeOfBasisMatrix();
  ASSERT_TRUE(result.has_value());
  const SKSymmetryCell& reduced = result->first;

  EXPECT_NEAR(reduced.a(), a, tolerance);
  EXPECT_NEAR(reduced.b(), b, tolerance);
  EXPECT_NEAR(reduced.c(), c, tolerance);
  EXPECT_NEAR(reduced.alpha() * radiansToDegrees, alphaDegrees, tolerance);
  EXPECT_NEAR(reduced.beta() * radiansToDegrees, betaDegrees, tolerance);
  EXPECT_NEAR(reduced.gamma() * radiansToDegrees, gammaDegrees, tolerance);

  // The returned change-of-basis matrix must map the input lattice onto the reduced cell.
  double3x3 reducedLattice = input.unitCell() * result->second;
  SKSymmetryCell fromLattice = SKSymmetryCell::createFromUnitCell(reducedLattice);
  EXPECT_NEAR(fromLattice.a(), reduced.a(), tolerance);
  EXPECT_NEAR(fromLattice.b(), reduced.b(), tolerance);
  EXPECT_NEAR(fromLattice.c(), reduced.c(), tolerance);
  EXPECT_NEAR(fromLattice.alpha(), reduced.alpha(), tolerance);
  EXPECT_NEAR(fromLattice.beta(), reduced.beta(), tolerance);
  EXPECT_NEAR(fromLattice.gamma(), reduced.gamma(), tolerance);
}
}  // namespace

// A simple-cubic lattice presented through a skew unimodular basis must reduce back to the cube.
TEST(NiggliReduction, SkewedSimpleCubicReducesToCube)
{
  double3x3 cubic(double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0));
  SKSymmetryCell skewed = SKSymmetryCell::createFromUnitCell(cubic * unimodularSkew());
  expectReducedCell(skewed, 1.0, 1.0, 1.0, 90.0, 90.0, 90.0);
}

// A face-centered-cubic lattice reduces to the rhombohedral cell with 60 degree angles and edge a/sqrt(2).
TEST(NiggliReduction, SkewedFccReducesToRhombohedral60)
{
  double3x3 fcc(double3(0.0, 0.5, 0.5), double3(0.5, 0.0, 0.5), double3(0.5, 0.5, 0.0));
  const double edge = std::sqrt(0.5);
  SKSymmetryCell skewed = SKSymmetryCell::createFromUnitCell(fcc * unimodularSkew());
  expectReducedCell(skewed, edge, edge, edge, 60.0, 60.0, 60.0);
}

// A body-centered-cubic lattice reduces to a rhombohedral cell with cos(angle) = -1/3 and edge sqrt(3)/2.
TEST(NiggliReduction, SkewedBccReducesToRhombohedral109)
{
  double3x3 bcc(double3(-0.5, 0.5, 0.5), double3(0.5, -0.5, 0.5), double3(0.5, 0.5, -0.5));
  const double edge = std::sqrt(3.0) / 2.0;
  const double angle = std::acos(-1.0 / 3.0) * radiansToDegrees;
  SKSymmetryCell skewed = SKSymmetryCell::createFromUnitCell(bcc * unimodularSkew());
  expectReducedCell(skewed, edge, edge, edge, angle, angle, angle);
}

// An already-reduced triclinic cell must be returned unchanged (idempotence).
TEST(NiggliReduction, AlreadyReducedIsIdempotent)
{
  SKSymmetryCell reducedInput(3.0, 4.0, 5.0, 80.0, 85.0, 88.0);
  expectReducedCell(reducedInput, 3.0, 4.0, 5.0, 80.0, 85.0, 88.0);
}

// The reduction must be scale invariant: scaling every lattice vector by a constant scales the reduced edge
// lengths by the same constant and leaves the reduced angles unchanged. This exercises the relative-tolerance
// behaviour that the fixed absolute epsilon would otherwise break for large or small cells.
TEST(NiggliReduction, IsScaleInvariant)
{
  double3x3 fcc(double3(0.0, 0.5, 0.5), double3(0.5, 0.0, 0.5), double3(0.5, 0.5, 0.0));
  double3x3 skewed = fcc * unimodularSkew();

  for (const double scale : {1.0e-3, 1.0, 1.0e3})
  {
    SKSymmetryCell cell = SKSymmetryCell::createFromUnitCell(scale * skewed);
    auto result = cell.computeReducedNiggliCellAndChangeOfBasisMatrix();
    ASSERT_TRUE(result.has_value()) << "scale " << scale;
    const SKSymmetryCell& reduced = result->first;

    EXPECT_NEAR(reduced.a(), scale * std::sqrt(0.5), 1.0e-9 * scale) << "scale " << scale;
    EXPECT_NEAR(reduced.b(), scale * std::sqrt(0.5), 1.0e-9 * scale) << "scale " << scale;
    EXPECT_NEAR(reduced.c(), scale * std::sqrt(0.5), 1.0e-9 * scale) << "scale " << scale;
    EXPECT_NEAR(reduced.alpha() * radiansToDegrees, 60.0, 1.0e-6) << "scale " << scale;
    EXPECT_NEAR(reduced.beta() * radiansToDegrees, 60.0, 1.0e-6) << "scale " << scale;
    EXPECT_NEAR(reduced.gamma() * radiansToDegrees, 60.0, 1.0e-6) << "scale " << scale;
  }
}

// Two unimodular-equivalent presentations of the same lattice must yield identical reduced cells.
TEST(NiggliReduction, InvariantUnderUnimodularBasisChange)
{
  double3x3 bcc(double3(-0.5, 0.5, 0.5), double3(0.5, -0.5, 0.5), double3(0.5, 0.5, -0.5));

  SKSymmetryCell primitive = SKSymmetryCell::createFromUnitCell(bcc);
  SKSymmetryCell skewed = SKSymmetryCell::createFromUnitCell(bcc * unimodularSkew());

  auto reducedPrimitive = primitive.computeReducedNiggliCellAndChangeOfBasisMatrix();
  auto reducedSkewed = skewed.computeReducedNiggliCellAndChangeOfBasisMatrix();
  ASSERT_TRUE(reducedPrimitive.has_value());
  ASSERT_TRUE(reducedSkewed.has_value());

  EXPECT_NEAR(reducedPrimitive->first.a(), reducedSkewed->first.a(), 1.0e-9);
  EXPECT_NEAR(reducedPrimitive->first.b(), reducedSkewed->first.b(), 1.0e-9);
  EXPECT_NEAR(reducedPrimitive->first.c(), reducedSkewed->first.c(), 1.0e-9);
  EXPECT_NEAR(reducedPrimitive->first.alpha(), reducedSkewed->first.alpha(), 1.0e-9);
  EXPECT_NEAR(reducedPrimitive->first.beta(), reducedSkewed->first.beta(), 1.0e-9);
  EXPECT_NEAR(reducedPrimitive->first.gamma(), reducedSkewed->first.gamma(), 1.0e-9);
}
