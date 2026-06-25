#include <gtest/gtest.h>

import std;

import double3;
import double3x3;

import skposcarlegacyparser;
import sksymmetrycell;
import skspacegroup;
import skpointgroup;
import skrotationmatrix;

std::pair<double3x3, std::vector<std::tuple<double3, std::size_t, double>>> loadPOSCAR(const std::string& fileName)
{
  std::ifstream t(fileName.c_str());
  std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

  SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
  parser.startParsing();

  return {parser.movies().front().front()->cell->unitCell(), parser.firstTestFrame()};
}

std::optional<std::size_t> findSpaceGroupNumber(double3x3 unitCell,
                                                const std::vector<std::tuple<double3, std::size_t, double>>& atoms,
                                                double symmetryPrecision)
{
  std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
      SKSpaceGroup::findSpaceGroup(unitCell, atoms, false, symmetryPrecision);
  if (!spaceGroup)
  {
    return std::nullopt;
  }
  return SKSpaceGroup(spaceGroup->HallNumber).spaceGroupSetting().number();
}

// Conventional cells with Z > 1 in the POSCAR must use the primitive Delaunay metric when matching symmetry.
TEST(FindSpacegroupSymmetryMetric, FaceCenteredCubicConventionalCell)
{
  const auto& [unitCell, atoms] = loadPOSCAR("spglibtestdata/cubic/POSCAR-225");

  std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms{};
  std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms),
               [](const std::tuple<double3, std::size_t, double>& atom) { return std::get<1>(atom) == 1; });

  double3x3 smallestUnitCell =
      SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, atoms, unitCell, false, 1e-5);
  std::optional<double3x3> primitiveDelaunayUnitCell =
      SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, 1e-5);

  ASSERT_TRUE(primitiveDelaunayUnitCell.has_value());
  EXPECT_LT(std::fabs(primitiveDelaunayUnitCell->determinant()), std::fabs(unitCell.determinant()) - 1.0);

  std::optional<std::size_t> spaceGroupNumber = findSpaceGroupNumber(unitCell, atoms, 1e-5);
  ASSERT_TRUE(spaceGroupNumber.has_value());
  EXPECT_EQ(*spaceGroupNumber, 225u);
}

TEST(FindSpacegroupSymmetryMetric, BodyCenteredTetragonalConventionalCell)
{
  const auto& [unitCell, atoms] = loadPOSCAR("spglibtestdata/tetragonal/POSCAR-083-3");

  std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms{};
  std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms),
               [](const std::tuple<double3, std::size_t, double>& atom) { return std::get<1>(atom) == 1; });

  double3x3 smallestUnitCell =
      SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, atoms, unitCell, false, 1e-5);
  std::optional<double3x3> primitiveDelaunayUnitCell =
      SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, 1e-5);

  ASSERT_TRUE(primitiveDelaunayUnitCell.has_value());
  EXPECT_LT(std::fabs(primitiveDelaunayUnitCell->determinant()), 0.26 * std::fabs(unitCell.determinant()));

  std::optional<std::size_t> spaceGroupNumber = findSpaceGroupNumber(unitCell, atoms, 1e-5);
  ASSERT_TRUE(spaceGroupNumber.has_value());
  EXPECT_EQ(*spaceGroupNumber, 83u);
}

TEST(FindSpacegroupSymmetryMetric, TestSymmetryUsesLatticeOfFractionalCoordinates)
{
  const auto& [unitCell, atoms] = loadPOSCAR("spglibtestdata/cubic/POSCAR-225");

  std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms{};
  std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms),
               [](const std::tuple<double3, std::size_t, double>& atom) { return std::get<1>(atom) == 1; });

  double3x3 smallestUnitCell =
      SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, atoms, unitCell, false, 1e-5);
  std::optional<double3x3> primitiveDelaunayUnitCell =
      SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, 1e-5);
  ASSERT_TRUE(primitiveDelaunayUnitCell.has_value());

  double3 direction(-0.5, -0.01282051, 0.5);
  double primitiveLength = (*primitiveDelaunayUnitCell * direction).length();
  ASSERT_GT(primitiveLength, 0.0);
  const double precision = 7e-6;
  double3 dr = direction * (0.8 * precision / primitiveLength);

  std::vector<std::tuple<double3, std::size_t, double>> singleAtom = {std::make_tuple(double3(0.0, 0.0, 0.0), 1, 1.0)};

  EXPECT_TRUE(SKSymmetryCell::testSymmetry(dr, SKRotationMatrix::identity, singleAtom, *primitiveDelaunayUnitCell,
                                         false, precision));
  EXPECT_FALSE(SKSymmetryCell::testSymmetry(dr, SKRotationMatrix::identity, singleAtom, unitCell, false, precision));
}

TEST(FindPointgroupSymmetryMetric, FaceCenteredCubicConventionalCell)
{
  const auto& [unitCell, atoms] = loadPOSCAR("spglibtestdata/cubic/POSCAR-225");

  std::optional<SKPointGroup> pointGroup = SKSpaceGroup::findPointGroup(unitCell, atoms, false, 1e-5);
  ASSERT_TRUE(pointGroup.has_value());
  EXPECT_EQ(pointGroup->number(), 32u);
}
