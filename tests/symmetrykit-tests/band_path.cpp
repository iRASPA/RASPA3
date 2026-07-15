#include <gtest/gtest.h>

import std;

import double3;
import double3x3;

import skbandpath;

namespace
{
// Build a cell (columns are lattice vectors) from cell parameters (angles in degrees).
double3x3 cell(double a, double b, double c, double alpha, double beta, double gamma)
{
  double d2r = std::numbers::pi / 180.0;
  double ca = std::cos(alpha * d2r);
  double cb = std::cos(beta * d2r);
  double cg = std::cos(gamma * d2r);
  double sg = std::sin(gamma * d2r);
  double vf = std::sqrt(std::max(0.0, 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg));
  return double3x3(double3(a, 0.0, 0.0), double3(b * cg, b * sg, 0.0),
                   double3(c * cb, c * (ca - cb * cg) / sg, c * vf / sg));
}

void expectPoint(const SKBandPath& path, const std::string& label, double x, double y, double z)
{
  std::optional<double3> p = path.coordinatesForLabel(label);
  ASSERT_TRUE(p.has_value()) << "missing k-point " << label;
  EXPECT_NEAR(p->x, x, 1e-9) << label << ".x";
  EXPECT_NEAR(p->y, y, 1e-9) << label << ".y";
  EXPECT_NEAR(p->z, z, 1e-9) << label << ".z";
}

// Every path segment endpoint must be a defined k-point.
void expectConsistentPath(const SKBandPath& path)
{
  for (const auto& [from, to] : path.segments)
  {
    EXPECT_TRUE(path.coordinatesForLabel(from).has_value()) << "segment start undefined: " << from;
    EXPECT_TRUE(path.coordinatesForLabel(to).has_value()) << "segment end undefined: " << to;
  }
  EXPECT_FALSE(path.segments.empty());
  expectPoint(path, "GAMMA", 0.0, 0.0, 0.0);
}
}  // namespace

TEST(BandPath, SimpleCubic)
{
  auto path = SKBrillouinZonePath(cell(4.0, 4.0, 4.0, 90, 90, 90), 221, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->bravaisLattice, "cP");
  EXPECT_EQ(path->extendedBravaisLattice, "cP2");
  expectPoint(*path, "X", 0.0, 0.5, 0.0);
  expectPoint(*path, "M", 0.5, 0.5, 0.0);
  expectPoint(*path, "R", 0.5, 0.5, 0.5);
  expectConsistentPath(*path);
}

TEST(BandPath, FaceCenteredCubic)
{
  auto path = SKBrillouinZonePath(cell(4.0, 4.0, 4.0, 90, 90, 90), 225, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->bravaisLattice, "cF");
  EXPECT_EQ(path->extendedBravaisLattice, "cF2");
  expectPoint(*path, "X", 0.5, 0.0, 0.5);
  expectPoint(*path, "L", 0.5, 0.5, 0.5);
  expectPoint(*path, "W", 0.5, 0.25, 0.75);
  expectPoint(*path, "K", 0.375, 0.375, 0.75);
  expectPoint(*path, "U", 0.625, 0.25, 0.625);
  expectConsistentPath(*path);
}

TEST(BandPath, BodyCenteredCubic)
{
  auto path = SKBrillouinZonePath(cell(4.0, 4.0, 4.0, 90, 90, 90), 229, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->extendedBravaisLattice, "cI1");
  expectPoint(*path, "H", 0.5, -0.5, 0.5);
  expectPoint(*path, "P", 0.25, 0.25, 0.25);
  expectPoint(*path, "N", 0.0, 0.0, 0.5);
  expectConsistentPath(*path);
}

TEST(BandPath, TetragonalPrimitive)
{
  auto path = SKBrillouinZonePath(cell(4.0, 4.0, 6.0, 90, 90, 90), 123, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->extendedBravaisLattice, "tP1");
  expectPoint(*path, "X", 0.0, 0.5, 0.0);
  expectPoint(*path, "M", 0.5, 0.5, 0.0);
  expectPoint(*path, "A", 0.5, 0.5, 0.5);
  expectConsistentPath(*path);
}

TEST(BandPath, TetragonalBodyCentered)
{
  // c <= a  -> tI1 ; c > a -> tI2
  auto flat = SKBrillouinZonePath(cell(5.0, 5.0, 4.0, 90, 90, 90), 139, true);
  ASSERT_TRUE(flat.has_value());
  EXPECT_EQ(flat->extendedBravaisLattice, "tI1");
  expectConsistentPath(*flat);

  auto tall = SKBrillouinZonePath(cell(4.0, 4.0, 6.0, 90, 90, 90), 139, true);
  ASSERT_TRUE(tall.has_value());
  EXPECT_EQ(tall->extendedBravaisLattice, "tI2");
  expectConsistentPath(*tall);
}

TEST(BandPath, Hexagonal)
{
  auto path = SKBrillouinZonePath(cell(3.0, 3.0, 5.0, 90, 90, 120), 194, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->bravaisLattice, "hP");
  EXPECT_EQ(path->extendedBravaisLattice, "hP2");
  expectPoint(*path, "K", 1.0 / 3.0, 1.0 / 3.0, 0.0);
  expectPoint(*path, "M", 0.5, 0.0, 0.0);
  expectPoint(*path, "A", 0.0, 0.0, 0.5);
  expectConsistentPath(*path);
}

TEST(BandPath, Rhombohedral)
{
  // sqrt(3) a <= sqrt(2) c  -> hR1
  auto tall = SKBrillouinZonePath(cell(3.0, 3.0, 8.0, 90, 90, 120), 166, true);
  ASSERT_TRUE(tall.has_value());
  EXPECT_EQ(tall->bravaisLattice, "hR");
  EXPECT_EQ(tall->extendedBravaisLattice, "hR1");
  expectConsistentPath(*tall);

  auto flat = SKBrillouinZonePath(cell(5.0, 5.0, 3.0, 90, 90, 120), 166, true);
  ASSERT_TRUE(flat.has_value());
  EXPECT_EQ(flat->extendedBravaisLattice, "hR2");
  expectConsistentPath(*flat);
}

TEST(BandPath, OrthorhombicBaseCentered)
{
  // a <= b -> oC1, else oC2
  auto oc1 = SKBrillouinZonePath(cell(3.0, 4.0, 5.0, 90, 90, 90), 65, true);
  ASSERT_TRUE(oc1.has_value());
  EXPECT_EQ(oc1->bravaisLattice, "oC");
  EXPECT_EQ(oc1->extendedBravaisLattice, "oC1");
  expectConsistentPath(*oc1);

  auto oc2 = SKBrillouinZonePath(cell(5.0, 4.0, 3.0, 90, 90, 90), 65, true);
  ASSERT_TRUE(oc2.has_value());
  EXPECT_EQ(oc2->extendedBravaisLattice, "oC2");
  expectConsistentPath(*oc2);
}

TEST(BandPath, OrthorhombicFaceCentered)
{
  auto path = SKBrillouinZonePath(cell(3.0, 4.0, 5.0, 90, 90, 90), 69, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->bravaisLattice, "oF");
  EXPECT_TRUE(path->extendedBravaisLattice == "oF1" || path->extendedBravaisLattice == "oF2" ||
              path->extendedBravaisLattice == "oF3");
  expectConsistentPath(*path);
}

TEST(BandPath, MonoclinicBaseCentered)
{
  auto path = SKBrillouinZonePath(cell(4.0, 3.0, 5.0, 90, 105, 90), 15, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->bravaisLattice, "mC");
  EXPECT_TRUE(path->extendedBravaisLattice == "mC1" || path->extendedBravaisLattice == "mC2" ||
              path->extendedBravaisLattice == "mC3");
  expectConsistentPath(*path);
}

TEST(BandPath, Triclinic)
{
  auto path = SKBrillouinZonePath(cell(4.0, 5.0, 6.0, 82, 95, 100), 2, true);
  ASSERT_TRUE(path.has_value());
  EXPECT_EQ(path->bravaisLattice, "aP");
  EXPECT_TRUE(path->extendedBravaisLattice == "aP2" || path->extendedBravaisLattice == "aP3");
  expectConsistentPath(*path);
}

TEST(BandPath, TimeReversalAugmentation)
{
  // Non-centrosymmetric cubic (P432, no. 207) with time reversal disabled should augment the path.
  auto path = SKBrillouinZonePath(cell(4.0, 4.0, 4.0, 90, 90, 90), 207, /*hasInversionSymmetry=*/false,
                                  /*withTimeReversal=*/false);
  ASSERT_TRUE(path.has_value());
  EXPECT_TRUE(path->timeReversalAugmented);
  EXPECT_TRUE(path->coordinatesForLabel("X'").has_value());
}

TEST(BandPath, ReciprocalMappingConsistency)
{
  // The Cartesian wavevector of a labelled point must be identical whether expressed in the primitive
  // reciprocal basis or mapped into the analyzed (here conventional) cell reciprocal basis.
  auto path = SKBrillouinZonePath(cell(4.0, 4.0, 4.0, 90, 90, 90), 225, true);
  ASSERT_TRUE(path.has_value());
  for (const std::string& label : {std::string("X"), std::string("L"), std::string("W"), std::string("K")})
  {
    double3 fPrimitive = path->coordinatesForLabel(label).value();
    double3 kFromPrimitive = path->primitiveReciprocalCell * fPrimitive;

    double3 fConventional = path->reciprocalFractionalInAnalyzedCell(fPrimitive);
    double3x3 conventionalReciprocal = (2.0 * std::numbers::pi) * path->conventionalCell.inverse().transpose();
    double3 kFromConventional = conventionalReciprocal * fConventional;

    EXPECT_NEAR(kFromPrimitive.x, kFromConventional.x, 1e-9) << label;
    EXPECT_NEAR(kFromPrimitive.y, kFromConventional.y, 1e-9) << label;
    EXPECT_NEAR(kFromPrimitive.z, kFromConventional.z, 1e-9) << label;
  }
}

TEST(BandPath, ReciprocalCubicMagnitude)
{
  double a = 4.0;
  auto path = SKBrillouinZonePath(cell(a, a, a, 90, 90, 90), 221, true);
  ASSERT_TRUE(path.has_value());
  double3 f = path->coordinatesForLabel("X").value();
  double3 kCart = path->primitiveReciprocalCell * f;
  // X = (0, 1/2, 0) -> |k| = pi/a
  EXPECT_NEAR(kCart.length(), std::numbers::pi / a, 1e-9);
}
