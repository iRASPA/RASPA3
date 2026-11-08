#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import skdefinitions;
import skcell;
import skrotationmatrix;
import skseitzintegermatrix;
import skspacegroupdatabase;
import skspacegroup;
import skstructure;
import skparser;
import skposcarparser;
import blocking_spheres;
import energy_shared_energy_backend;
import energy_shared_blocking_mask;

TEST(SYMMETRY_STRUCTURE_FIXES, set_box_divides_the_box_over_the_replicas)
{
  SKCell cell(10.0, 10.0, 10.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0);
  cell.setMaximumReplicas(int3(1, 0, 0));  // two replicas along a
  cell.setBox(double3x3(double3(30.0, 0.0, 0.0), double3(0.0, 12.0, 0.0), double3(0.0, 0.0, 14.0)));
  EXPECT_NEAR(cell.unitCell().ax, 15.0, 1e-12);
  EXPECT_NEAR(cell.unitCell().by, 12.0, 1e-12);
  EXPECT_NEAR(cell.inverseUnitCell().ax, 1.0 / 15.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, angle_setters_keep_the_inverse_unit_cell)
{
  SKCell cell(10.0, 11.0, 12.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0);
  cell.setAlphaAngle(80.0 * std::numbers::pi / 180.0);
  const double3x3 product = cell.unitCell() * cell.inverseUnitCell();
  EXPECT_NEAR(product.ax, 1.0, 1e-12);
  EXPECT_NEAR(product.bz, 0.0, 1e-12);
  EXPECT_NEAR(product.cz, 1.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, equal_seitz_matrices_hash_equally)
{
  const SKSeitzIntegerMatrix a(SKRotationMatrix::identity, int3(1, 2, 3));
  const SKSeitzIntegerMatrix b(SKRotationMatrix::identity, int3(25, 26, 27));
  ASSERT_TRUE(a == b);
  EXPECT_EQ(SKSeitzIntegerMatrix::hashFunction{}(a), SKSeitzIntegerMatrix::hashFunction{}(b));
}

TEST(SYMMETRY_STRUCTURE_FIXES, b_centred_settings_are_b_face_centred)
{
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[13].HallNumber(), 13uz);
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[13].centring(), Centring::b_face);
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[15].HallNumber(), 15uz);
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[15].centring(), Centring::b_face);
}

TEST(SYMMETRY_STRUCTURE_FIXES, space_group_search_without_atoms_finds_nothing)
{
  const double3x3 unitCell(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  EXPECT_FALSE(SKSpaceGroup::findSpaceGroup(unitCell, {}, false, 1e-2).has_value());
}

TEST(SYMMETRY_STRUCTURE_FIXES, poscar_blank_comment_negative_volume_and_scaled_cartesian_coordinates)
{
  // blank comment line, a negative scaling factor (a cell volume of 1000 A^3 from a 2x2x2 cube), Cartesian
  // coordinates scaled like the lattice
  const std::string content =
      "\n"
      "-1000.0\n"
      "2.0 0.0 0.0\n"
      "0.0 2.0 0.0\n"
      "0.0 0.0 2.0\n"
      "Si\n"
      "1\n"
      "Cartesian\n"
      "0.5 1.0 1.5\n";
  SKPOSCARParser parser(content);
  parser.startParsing();
  const double3x3 cell = parser.movies().front().front()->cell->unitCell();
  EXPECT_NEAR(cell.ax, 10.0, 1e-9);
  EXPECT_NEAR(cell.determinant(), 1000.0, 1e-6);

  const std::vector<std::tuple<double3, std::size_t, double>> atoms = parser.firstTestFrame();
  ASSERT_EQ(atoms.size(), 1uz);
  // the Cartesian position (0.5, 1, 1.5) scaled by 5 is (2.5, 5, 7.5): fractional (0.25, 0.5, 0.75)
  EXPECT_NEAR(std::get<0>(atoms[0]).x, 0.25, 1e-12);
  EXPECT_NEAR(std::get<0>(atoms[0]).y, 0.5, 1e-12);
  EXPECT_NEAR(std::get<0>(atoms[0]).z, 0.75, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, a_backend_with_blocked_pockets_is_another_backend_for_the_field_cache)
{
  EnergyBackend backend;
  backend.name = "cpu";
  const std::vector<BlockingSphere> first{BlockingSphere{double3(0.1, 0.2, 0.3), 1.5}};
  const std::vector<BlockingSphere> second{BlockingSphere{double3(0.1, 0.2, 0.3), 2.0}};
  const EnergyBackend blocked = withBlockedPockets(backend, first, 1.0, 2.0);
  const EnergyBackend otherBlocked = withBlockedPockets(backend, second, 1.0, 2.0);
  EXPECT_TRUE(backend.blocking.empty());
  EXPECT_FALSE(blocked.blocking.empty());
  EXPECT_NE(blocked.blocking, otherBlocked.blocking);
  EXPECT_EQ(blocked.name, backend.name);
}
