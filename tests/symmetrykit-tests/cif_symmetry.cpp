#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import skdefinitions;
import skseitzmatrix;
import skspacegroup;

std::string readFile(const std::string& path)
{
  std::ifstream input(path);
  return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

std::vector<std::string> readSymmetryOperations(const std::string& content)
{
  std::vector<std::string> operations;
  std::istringstream stream(content);
  std::string line;
  bool inLoop = false;

  while (std::getline(stream, line))
  {
    if (line.starts_with("loop_"))
    {
      inLoop = false;
      continue;
    }
    if (line == "_symmetry_equiv_pos_as_xyz")
    {
      inLoop = true;
      continue;
    }
    if (!inLoop || line.empty() || line.starts_with("#"))
    {
      continue;
    }
    if (line.starts_with("_"))
    {
      break;
    }
    operations.push_back(line);
  }

  return operations;
}

double3x3 cubicLattice(double length)
{
  return double3x3(double3(length, 0.0, 0.0), double3(0.0, length, 0.0), double3(0.0, 0.0, length));
}

TEST(CifSymmetry, ParsesBasicSymmetryOperation)
{
  const std::optional<SKSeitzMatrix> identity = SKSpaceGroup::parseCifSymmetryOperation("'+x,+y,+z'");
  ASSERT_TRUE(identity.has_value());
  EXPECT_EQ(identity->rotation.int3x3_m.m11, 1);
  EXPECT_EQ(identity->rotation.int3x3_m.m22, 1);
  EXPECT_EQ(identity->rotation.int3x3_m.m33, 1);
  EXPECT_DOUBLE_EQ(identity->translation.x, 0.0);
  EXPECT_DOUBLE_EQ(identity->translation.y, 0.0);
  EXPECT_DOUBLE_EQ(identity->translation.z, 0.0);

  const std::optional<SKSeitzMatrix> translated = SKSpaceGroup::parseCifSymmetryOperation("'+x,1/2+y,1/2+z'");
  ASSERT_TRUE(translated.has_value());
  EXPECT_DOUBLE_EQ(translated->translation.x, 0.0);
  EXPECT_DOUBLE_EQ(translated->translation.y, 0.5);
  EXPECT_DOUBLE_EQ(translated->translation.z, 0.5);
}

TEST(CifSymmetry, IdentifiesLtnFd3mFromSymmetryOperations)
{
  const std::string content = readFile("testdata/LTN_symmetry.cif");
  const std::vector<std::string> operationStrings = readSymmetryOperations(content);
  ASSERT_EQ(operationStrings.size(), 192U);

  std::vector<SKSeitzMatrix> operations;
  operations.reserve(operationStrings.size());
  for (const std::string& operationString : operationStrings)
  {
    const std::optional<SKSeitzMatrix> operation = SKSpaceGroup::parseCifSymmetryOperation(operationString);
    ASSERT_TRUE(operation.has_value()) << operationString;
    operations.push_back(*operation);
  }

  const std::optional<std::size_t> hallNumber = SKSpaceGroup::HallNumberFromSymmetryOperations(
      operations, cubicLattice(35.622), Centring::face, 227U);
  ASSERT_TRUE(hallNumber.has_value());
  EXPECT_EQ(*hallNumber, 526U);
  EXPECT_TRUE(SKSpaceGroup(*hallNumber).spaceGroupSetting().standardSetting());
}

TEST(CifSymmetry, CentringFromHMSymbol)
{
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("'F d -3 m'"), Centring::face);
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("P 1"), Centring::primitive);
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("Ia-3d"), Centring::body);
}
