#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import unit_cell;

// The periodic cell: the shortest image of a displacement, whatever the shape of the cell and however it is laid
// against the axes.
//
// The shortest image is checked against a search that owes nothing to how the cell finds it: every lattice
// translation out to as many cells as could hold a shorter vector, each one tried. A fractional step of t along
// an axis is at least |t| times the perpendicular width across it, so a vector no longer than the one in hand has
// no step of more than that length over the width, plus the half the rounding left.

namespace
{

double3 shortestBySearch(const double3x3 &cell, const double3 &dr)
{
  const UnitCell measured(cell, UnitCell::Type::Triclinic);
  const double3 widths = measured.perpendicularWidths();
  double3 s = cell.inverse() * dr;
  s = double3(s.x - std::round(s.x), s.y - std::round(s.y), s.z - std::round(s.z));
  const double3 base = cell * s;
  const double reach = base.length();
  const int nx = static_cast<int>(std::floor(reach / widths.x + 0.5)) + 1;
  const int ny = static_cast<int>(std::floor(reach / widths.y + 0.5)) + 1;
  const int nz = static_cast<int>(std::floor(reach / widths.z + 0.5)) + 1;

  double3 best = base;
  for (int i = -nx; i <= nx; ++i)
    for (int j = -ny; j <= ny; ++j)
      for (int k = -nz; k <= nz; ++k)
      {
        const double3 image = cell * double3(s.x + i, s.y + j, s.z + k);
        if (image.length_squared() < best.length_squared()) best = image;
      }
  return best;
}

// A displacement is an image of another when they differ by whole cell vectors.
bool differByALatticeVector(const double3x3 &cell, const double3 &a, const double3 &b)
{
  const double3 s = cell.inverse() * (a - b);
  return std::abs(s.x - std::round(s.x)) < 1.0e-8 && std::abs(s.y - std::round(s.y)) < 1.0e-8 &&
         std::abs(s.z - std::round(s.z)) < 1.0e-8;
}

struct Draw
{
  std::uint64_t state;
  double operator()()
  {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>(state >> 11) / 9007199254740992.0;
  }
};

void expectShortestEverywhere(const UnitCell &unitCell, const double3x3 &cell, std::uint64_t seed, const char *what)
{
  Draw draw{seed};
  for (int sample = 0; sample < 400; ++sample)
  {
    const double3 dr = cell * double3(6.0 * draw() - 3.0, 6.0 * draw() - 3.0, 6.0 * draw() - 3.0);
    const double3 image = unitCell.applyPeriodicBoundaryConditions(dr);
    const double3 shortest = shortestBySearch(cell, dr);
    ASSERT_NEAR(image.length(), shortest.length(), 1.0e-9 * (1.0 + shortest.length()))
        << what << ": the image of (" << dr.x << ", " << dr.y << ", " << dr.z << ") is not the shortest";
    ASSERT_TRUE(differByALatticeVector(cell, image, dr)) << what << ": what came back is not an image at all";
  }
}

double3x3 rotationAboutZ(double angle)
{
  return double3x3(double3(std::cos(angle), std::sin(angle), 0.0), double3(-std::sin(angle), std::cos(angle), 0.0),
                   double3(0.0, 0.0, 1.0));
}

constexpr double degrees = std::numbers::pi / 180.0;

}  // namespace

// The plainest case of it: a hexagonal net, with the two axes at sixty degrees. Rounding the fractional
// coordinates keeps (0.45, 0.45), which is 7.8 long, where the image one cell back along a is 5.1.
TEST(unit_cell, the_image_nearest_in_fractions_is_not_the_shortest_in_a_skewed_cell)
{
  const UnitCell unitCell(10.0, 10.0, 10.0, 90.0 * degrees, 90.0 * degrees, 60.0 * degrees);
  const double3 dr = unitCell.cell * double3(0.45, 0.45, 0.0);
  const double3 image = unitCell.applyPeriodicBoundaryConditions(dr);
  EXPECT_NEAR(image.length(), (unitCell.cell * double3(-0.55, 0.45, 0.0)).length(), 1.0e-12);
  EXPECT_LT(image.length(), 5.2);
}

TEST(unit_cell, the_image_is_the_shortest_in_cells_of_every_shape)
{
  Draw draw{0x1234567ull};
  for (int shape = 0; shape < 150; ++shape)
  {
    const double a = 5.0 + 25.0 * draw();
    const double b = 5.0 + 25.0 * draw();
    const double c = 5.0 + 25.0 * draw();
    const double alpha = (35.0 + 110.0 * draw()) * degrees;
    const double beta = (35.0 + 110.0 * draw()) * degrees;
    const double gamma = (35.0 + 110.0 * draw()) * degrees;
    // three angles are those of a cell only if the cell they make has volume
    const double inside = 1.0 - std::cos(alpha) * std::cos(alpha) - std::cos(beta) * std::cos(beta) -
                          std::cos(gamma) * std::cos(gamma) + 2.0 * std::cos(alpha) * std::cos(beta) * std::cos(gamma);
    if (inside < 0.05) continue;

    const UnitCell unitCell(a, b, c, alpha, beta, gamma);
    expectShortestEverywhere(unitCell, unitCell.cell, 77ull + static_cast<std::uint64_t>(shape), "a triclinic cell");
  }
}

// A lattice has many cells. The same cubic lattice described by three long slanted vectors is the cell a reduction
// was never run on, and the images in it are the cubic ones.
TEST(unit_cell, the_image_does_not_depend_on_which_cell_of_the_lattice_was_given)
{
  const double3x3 cubic(double3(8.0, 0.0, 0.0), double3(0.0, 8.0, 0.0), double3(0.0, 0.0, 8.0));
  const double3x3 slanted(cubic * double3(1.0, 0.0, 0.0), cubic * double3(3.0, 1.0, 0.0),
                          cubic * double3(-2.0, 4.0, 1.0));
  const UnitCell box(8.0, 8.0, 8.0);
  const UnitCell unitCell(slanted);
  EXPECT_NEAR(std::abs(unitCell.volume), 512.0, 1.0e-9);

  Draw draw{99ull};
  for (int sample = 0; sample < 400; ++sample)
  {
    const double3 dr(40.0 * draw() - 20.0, 40.0 * draw() - 20.0, 40.0 * draw() - 20.0);
    const double3 image = unitCell.applyPeriodicBoundaryConditions(dr);
    const double3 inTheBox = box.applyPeriodicBoundaryConditions(dr);
    ASSERT_NEAR(image.length(), inTheBox.length(), 1.0e-9);
    ASSERT_TRUE(differByALatticeVector(slanted, image, dr));
  }
}

// A box is a box however it lies. One turned about an axis, one with its axes named in another order, and one
// with an axis pointing backwards all have right angles, and none of them is the matrix with the three lengths
// down its diagonal that the quick form of the convention reads.
TEST(unit_cell, a_box_that_is_not_along_the_axes_is_still_wrapped_as_the_box_it_is)
{
  const double3x3 box(double3(10.0, 0.0, 0.0), double3(0.0, 12.0, 0.0), double3(0.0, 0.0, 15.0));

  const double3x3 turned = rotationAboutZ(30.0 * degrees) * box;
  expectShortestEverywhere(UnitCell(turned), turned, 5ull, "a box turned about z");
  expectShortestEverywhere(UnitCell(turned, UnitCell::Type::Rectangular), turned, 6ull,
                           "a box turned about z and called rectangular");

  const double3x3 renamed(double3(0.0, 12.0, 0.0), double3(0.0, 0.0, 15.0), double3(10.0, 0.0, 0.0));
  expectShortestEverywhere(UnitCell(renamed), renamed, 7ull, "a box with its axes in another order");

  const double3x3 backwards(double3(-10.0, 0.0, 0.0), double3(0.0, 12.0, 0.0), double3(0.0, 0.0, 15.0));
  expectShortestEverywhere(UnitCell(backwards), backwards, 8ull, "a box with an axis pointing backwards");

  // and one that is along the axes is wrapped the quick way, as it always was
  const UnitCell plain(box);
  EXPECT_EQ(plain.type, UnitCell::Type::Rectangular);
  expectShortestEverywhere(plain, box, 9ull, "a box along the axes");
}

// Which image it is: the whole cell vectors taken off, for a caller that has to know where the other end went.
TEST(unit_cell, the_translation_taken_off_is_the_one_that_leaves_the_shortest_image)
{
  for (const UnitCell &unitCell : {UnitCell(10.0, 12.0, 15.0),
                                   UnitCell(10.0, 10.0, 10.0, 90.0 * degrees, 90.0 * degrees, 60.0 * degrees),
                                   UnitCell(7.0, 9.0, 11.0, 50.0 * degrees, 70.0 * degrees, 115.0 * degrees)})
  {
    Draw draw{31ull};
    for (int sample = 0; sample < 200; ++sample)
    {
      const double3 dr = unitCell.cell * double3(6.0 * draw() - 3.0, 6.0 * draw() - 3.0, 6.0 * draw() - 3.0);
      const int3 translation = unitCell.translationToShortestImage(dr);
      const double3 left = dr - unitCell.cell * double3(translation.x, translation.y, translation.z);
      ASSERT_NEAR(left.length(), shortestBySearch(unitCell.cell, dr).length(), 1.0e-9);
      ASSERT_NEAR((left - unitCell.applyPeriodicBoundaryConditions(dr)).length(), 0.0, 1.0e-9);
    }
  }
}

// The angles are read off the vectors with an arc cosine, which has no value a rounding error past one or minus
// one, and none at all for a vector of no length.
TEST(unit_cell, the_angles_of_a_cell_are_angles_whatever_the_cell)
{
  auto areAngles = [](const UnitCell &unitCell)
  {
    for (double angle : {unitCell.angleAlpha, unitCell.angleBeta, unitCell.angleGamma})
    {
      if (!(angle >= 0.0 && angle <= std::numbers::pi)) return false;
    }
    return true;
  };

  // two vectors along one line, which in floating point can come to a cosine past one
  const double3 along(0.1, 0.2, 0.3);
  EXPECT_TRUE(areAngles(UnitCell(double3x3(along, 3.0 * along, double3(0.0, 0.0, 1.0)))));
  EXPECT_TRUE(areAngles(UnitCell(double3x3(along, -7.0 * along, double3(0.0, 0.0, 1.0)))));

  // a vector of no length
  EXPECT_TRUE(areAngles(UnitCell(double3x3(double3(0.0, 0.0, 0.0), double3(0.0, 5.0, 0.0), double3(0.0, 0.0, 5.0)))));
  EXPECT_TRUE(areAngles(UnitCell(double3x3())));

  // and an honest cell keeps the angles it has
  const UnitCell unitCell(7.0, 9.0, 11.0, 50.0 * degrees, 70.0 * degrees, 115.0 * degrees);
  const UnitCell again(unitCell.cell);
  EXPECT_NEAR(again.angleAlpha, 50.0 * degrees, 1.0e-12);
  EXPECT_NEAR(again.angleBeta, 70.0 * degrees, 1.0e-12);
  EXPECT_NEAR(again.angleGamma, 115.0 * degrees, 1.0e-12);
}

// A cell repeated is the same lattice coarser, and wraps as one.
TEST(unit_cell, a_repeated_cell_wraps_as_the_larger_cell_it_is)
{
  const UnitCell unitCell(7.0, 9.0, 11.0, 50.0 * degrees, 70.0 * degrees, 115.0 * degrees);
  const UnitCell repeated = unitCell.scaled(int3(2, 1, 3));
  EXPECT_NEAR(repeated.volume, 6.0 * unitCell.volume, 1.0e-9 * unitCell.volume);
  expectShortestEverywhere(repeated, repeated.cell, 41ull, "a repeated cell");
}

// A cell repeated stands where the cell stood. Its vectors are whole multiples of the cell's own, so an atom at
// some fraction of the cell is at that fraction, divided by the count, of the repeated one, and at the same
// place in space; and it is the same hand. The repeated cell used to be rebuilt from the lengths and the angles
// alone, which put it with its first vector along x and its second in the xy-plane whatever way the cell itself
// was turned: the same shape, somewhere else, and right-handed whether or not the cell was.
TEST(unit_cell, a_repeated_cell_stands_where_the_cell_stood)
{
  const int3 counts(2, 1, 3);

  // The cell as its lengths and angles put it, the same cell turned, and the same with two vectors changed over.
  const UnitCell canonical(7.0, 9.0, 11.0, 50.0 * degrees, 70.0 * degrees, 115.0 * degrees);
  auto turned = [](const double3 &v)
  {
    const double aboutZ = 0.7;
    const double aboutX = 0.4;
    const double3 first(std::cos(aboutZ) * v.x - std::sin(aboutZ) * v.y,
                        std::sin(aboutZ) * v.x + std::cos(aboutZ) * v.y, v.z);
    return double3(first.x, std::cos(aboutX) * first.y - std::sin(aboutX) * first.z,
                   std::sin(aboutX) * first.y + std::cos(aboutX) * first.z);
  };
  const std::vector<std::pair<std::string, UnitCell>> cells{
      {"as its lengths and angles put it", canonical},
      {"turned", UnitCell(double3x3(turned(canonical.cell[0]), turned(canonical.cell[1]), turned(canonical.cell[2])))},
      {"left-handed", UnitCell(double3x3(canonical.cell[1], canonical.cell[0], canonical.cell[2]))},
      {"a box along the axes", UnitCell(7.0, 9.0, 11.0)},
  };

  for (const auto &[name, unitCell] : cells)
  {
    const UnitCell repeated = unitCell.scaled(counts);

    // Each vector the cell's own, so many times over.
    EXPECT_LT((repeated.cell[0] - 2.0 * unitCell.cell[0]).length(), 1.0e-12) << name;
    EXPECT_LT((repeated.cell[1] - 1.0 * unitCell.cell[1]).length(), 1.0e-12) << name;
    EXPECT_LT((repeated.cell[2] - 3.0 * unitCell.cell[2]).length(), 1.0e-12) << name;

    // So a point of the cell is where it was.
    for (const double3 &fractional : {double3(0.1, 0.2, 0.3), double3(0.9, 0.5, 0.05), double3(1.7, -0.4, 2.2)})
    {
      const double3 inRepeated(fractional.x / 2.0, fractional.y / 1.0, fractional.z / 3.0);
      EXPECT_LT((repeated.cell * inRepeated - unitCell.cell * fractional).length(), 1.0e-12) << name;
    }

    // The same hand, six times the room, and the lengths and angles it had, to the last digit.
    EXPECT_GT(repeated.cell.determinant() * unitCell.cell.determinant(), 0.0) << name;
    EXPECT_NEAR(repeated.volume, 6.0 * unitCell.volume, 1.0e-9 * unitCell.volume) << name;
    EXPECT_DOUBLE_EQ(repeated.lengthA, 2.0 * unitCell.lengthA) << name;
    EXPECT_DOUBLE_EQ(repeated.lengthB, 1.0 * unitCell.lengthB) << name;
    EXPECT_DOUBLE_EQ(repeated.lengthC, 3.0 * unitCell.lengthC) << name;
    EXPECT_DOUBLE_EQ(repeated.angleAlpha, unitCell.angleAlpha) << name;
    EXPECT_DOUBLE_EQ(repeated.angleBeta, unitCell.angleBeta) << name;
    EXPECT_DOUBLE_EQ(repeated.angleGamma, unitCell.angleGamma) << name;
    EXPECT_EQ(repeated.type, unitCell.type) << name;

    expectShortestEverywhere(repeated, repeated.cell, 41ull, name.c_str());
  }
}

// The same in cells nobody would choose: lengths from two to sixty, angles from twenty degrees to a hundred and
// sixty, and each then sheared by whole numbers of its own vectors a few times over, which leaves the lattice
// what it was and the cell long and thin. The search it is held against is run in the cell before the shearing.
TEST(unit_cell, the_image_is_the_shortest_in_cells_nobody_would_choose)
{
  Draw draw{0xABCDEFull};
  std::size_t cells = 0;
  for (int shape = 0; shape < 600; ++shape)
  {
    const double a = 2.0 + 58.0 * draw(), b = 2.0 + 58.0 * draw(), c = 2.0 + 58.0 * draw();
    const double alpha = (20.0 + 140.0 * draw()) * degrees, beta = (20.0 + 140.0 * draw()) * degrees,
                 gamma = (20.0 + 140.0 * draw()) * degrees;
    const double inside = 1.0 - std::cos(alpha) * std::cos(alpha) - std::cos(beta) * std::cos(beta) -
                          std::cos(gamma) * std::cos(gamma) + 2.0 * std::cos(alpha) * std::cos(beta) * std::cos(gamma);
    if (inside < 0.02) continue;
    double3x3 cell = UnitCell(a, b, c, alpha, beta, gamma).cell;
    // shear it by whole numbers a few times: the same lattice, a worse cell
    for (int shear = 0; shear < static_cast<int>(4.0 * draw()); ++shear)
    {
      const int from = static_cast<int>(3.0 * draw()), to = (from + 1 + static_cast<int>(2.0 * draw())) % 3;
      const double times = std::round(4.0 * draw() - 2.0);
      std::array<double3, 3> columns{cell[0], cell[1], cell[2]};
      columns[static_cast<std::size_t>(to)] += times * columns[static_cast<std::size_t>(from)];
      cell = double3x3(columns[0], columns[1], columns[2]);
    }
    const UnitCell unitCell(cell);
    const UnitCell reference(UnitCell(a, b, c, alpha, beta, gamma).cell, UnitCell::Type::Triclinic);
    ++cells;
    for (int sample = 0; sample < 60; ++sample)
    {
      const double3 dr = reference.cell * double3(4.0 * draw() - 2.0, 4.0 * draw() - 2.0, 4.0 * draw() - 2.0);
      const double3 image = unitCell.applyPeriodicBoundaryConditions(dr);
      const double3 shortest = shortestBySearch(reference.cell, dr);
      ASSERT_NEAR(image.length(), shortest.length(), 1.0e-9 * (1.0 + shortest.length()))
          << "cell " << shape << " lengths " << a << " " << b << " " << c << " angles " << alpha / degrees << " "
          << beta / degrees << " " << gamma / degrees;
      ASSERT_TRUE(differByALatticeVector(reference.cell, image, dr));
      const int3 t = unitCell.translationToShortestImage(dr);
      ASSERT_NEAR((dr - unitCell.cell * double3(t.x, t.y, t.z) - image).length(), 0.0, 1.0e-8);
    }
  }
  EXPECT_GT(cells, 200uz);
}
