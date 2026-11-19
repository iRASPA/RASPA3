#include <gtest/gtest.h>

import std;

import double3;
import double3x3;

import skbandpath;
import skspacegroup;

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

namespace
{
// how far outside the first Brillouin zone the worst labelled point lies (relative; 0 when all are inside or on it)
double worstOutside(const SKBandPath& path, std::string* which = nullptr)
{
  double worst = 0.0;
  for (const SKKVector& kv : path.kPoints)
  {
    const double3 k = path.primitiveReciprocalCell * kv.coordinates;
    for (int i = -3; i <= 3; ++i)
      for (int j = -3; j <= 3; ++j)
        for (int l = -3; l <= 3; ++l)
        {
          if (i == 0 && j == 0 && l == 0) continue;
          const double3 G = path.primitiveReciprocalCell * double3(i, j, l);
          const double excess = (k.length_squared() - (k - G).length_squared()) / G.length_squared();
          if (excess > worst)
          {
            worst = excess;
            if (which) *which = kv.label;
          }
        }
  }
  return worst;
}

std::vector<std::tuple<double3, std::size_t, double>> atomsOf(std::size_t hallNumber)
{
  std::vector<std::tuple<double3, std::size_t, double>> atoms;
  std::size_t type = 1;
  for (const double3& position : {double3(0.113, 0.207, 0.341), double3(0.412, 0.069, 0.158)})
  {
    for (const double3& image : SKSpaceGroup(hallNumber).listOfSymmetricPositions(position))
      atoms.emplace_back(double3::fract(image), type, 1.0);
    ++type;
  }
  return atoms;
}
}  // namespace

// The labelled points lie in or on the first Brillouin zone for the standard cells of the lattices with a
// case distinction.
TEST(BandPath, LabelledPointsLieInTheFirstBrillouinZone)
{
  struct Case
  {
    double3x3 cell;
    std::size_t number;
  };
  for (const Case& c : {Case{cell(4, 4, 4, 90, 90, 90), 225}, Case{cell(3, 4, 5, 90, 90, 90), 69},
                        Case{cell(3, 10, 10, 90, 90, 90), 42}, Case{cell(10, 10, 3, 90, 90, 90), 42},
                        Case{cell(4, 3, 5, 90, 105, 90), 15}, Case{cell(4, 7, 9, 90, 100, 90), 14},
                        Case{cell(3, 4, 5, 90, 90, 90), 71}, Case{cell(3, 4, 5, 90, 90, 90), 65},
                        Case{cell(4, 3, 5, 90, 90, 90), 65}, Case{cell(3, 4, 5, 90, 90, 90), 38},
                        Case{cell(3, 5, 4, 90, 90, 90), 38}, Case{cell(3, 3, 5, 90, 90, 90), 139},
                        Case{cell(5, 5, 3, 90, 90, 90), 139}, Case{cell(4, 4, 10, 90, 90, 120), 166},
                        Case{cell(6, 6, 4, 90, 90, 120), 166}, Case{cell(4, 5, 6, 82, 95, 100), 2}})
  {
    const auto path = SKBrillouinZonePath(c.cell, c.number, true);
    ASSERT_TRUE(path.has_value());
    std::string which;
    EXPECT_LT(worstOutside(*path, &which), 1.0e-9) << path->extendedBravaisLattice << " " << which;
  }
}

// A cell without volume has no path (it used to come back with cells that are not numbers).
TEST(BandPath, ACellWithoutVolumeHasNoPath)
{
  EXPECT_FALSE(
      SKBrillouinZonePath(double3x3(double3(1, 0, 0), double3(2, 0, 0), double3(0, 0, 1)), 2, true).has_value());
  EXPECT_FALSE(
      SKBrillouinZonePath(double3x3(double3(1, 0, 0), double3(0, 1, 0), double3(0, 0, 0)), 47, true).has_value());
}

// The path of a structure is given in the Cartesian frame of the cell that was passed in: the conventional cell
// times the reported transformation is that cell, also for a crystal whose conventional axes are a permutation of
// the given ones, or that is turned in space, and for cells in a setting the tables do not take as it is.
TEST(BandPath, ThePathOfAStructureIsInTheFrameOfItsCell)
{
  struct Case
  {
    const char* name;
    double3x3 cell;
    std::size_t hallNumber;
  };
  const double3x3 rotation(double3(0.36, 0.48, -0.8), double3(-0.8, 0.6, 0.0), double3(0.48, 0.64, 0.6));
  for (const Case& c :
       {Case{"P 21/c, beta 118", cell(4, 7, 9, 90, 118, 90), 81},
        Case{"C 2/c, beta 123.3", cell(10.35, 6.42, 11.81, 90, 123.3, 90), 90},
        Case{"F m m 2, b shortest", cell(10, 3, 10, 90, 90, 90), 212}, Case{"F m m m, c < b < a", cell(5, 4, 3, 90, 90, 90), 334},
        Case{"F m m m, turned", rotation * cell(5, 4, 3, 90, 90, 90), 334},
        Case{"P 21/c, turned", rotation * cell(4, 7, 9, 90, 118, 90), 81}})
  {
    const auto path = SKBrillouinZonePath(c.cell, atomsOf(c.hallNumber), false, 1.0e-4);
    ASSERT_TRUE(path.has_value()) << c.name;

    const double3x3 rebuilt = path->conventionalCell * path->conventionalToAnalyzedTransformation;
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_LT((rebuilt[i] - c.cell[i]).length(), 1.0e-6) << c.name << ", cell vector " << i;
    }
    std::string which;
    EXPECT_LT(worstOutside(*path, &which), 1.0e-9) << c.name << " " << which;
  }
}
