#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import skdefinitions;
import skseitzmatrix;
import skseitzintegermatrix;
import skspacegroupsetting;
import skrotationmatrix;
import sksymmetrycell;
import skintegersymmetryoperationset;
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

namespace
{
std::vector<SKSeitzMatrix> operationsOfSetting(std::size_t hallNumber, double noise)
{
  std::vector<SKSeitzMatrix> operations;
  for (const SKSeitzIntegerMatrix& operation :
       SKSpaceGroup(hallNumber).spaceGroupSetting().fullSeitzMatrices().operations)
  {
    operations.emplace_back(operation.rotation,
                            double3(operation.translation.x / 24.0 + noise, operation.translation.y / 24.0 + noise,
                                    operation.translation.z / 24.0 + noise));
  }
  return operations;
}

double3x3 latticeOfSetting(std::size_t number)
{
  const double pi = std::numbers::pi;
  auto cell = [](double a, double b, double c, double alpha, double beta, double gamma)
  {
    const double ca = std::cos(alpha), cb = std::cos(beta), cg = std::cos(gamma), sg = std::sin(gamma);
    const double v = std::sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg);
    return double3x3(double3(a, 0.0, 0.0), double3(b * cg, b * sg, 0.0),
                     double3(c * cb, c * (ca - cb * cg) / sg, c * v / sg));
  };
  if (number <= 2) return cell(5.0, 6.0, 7.0, 0.4 * pi, 0.45 * pi, 0.55 * pi);
  if (number <= 15) return cell(5.0, 6.0, 7.0, 0.5 * pi, 0.6 * pi, 0.5 * pi);
  if (number <= 74) return cell(5.0, 6.0, 7.0, 0.5 * pi, 0.5 * pi, 0.5 * pi);
  if (number <= 142) return cell(5.0, 5.0, 7.0, 0.5 * pi, 0.5 * pi, 0.5 * pi);
  if (number <= 194) return cell(5.0, 5.0, 7.0, 0.5 * pi, 0.5 * pi, 2.0 * pi / 3.0);
  return cell(5.0, 5.0, 5.0, 0.5 * pi, 0.5 * pi, 0.5 * pi);
}
}  // namespace

// The operations of every setting, as a CIF lists them with translations rounded a little up or down, are read as
// that setting and not as another setting of the group: the origin shift of the setting itself is zero, and it used
// to come out as 1/2 or 1/3 (and as a centring translation) when a translation difference fell just below zero.
TEST(CifSymmetry, RoundedTranslationsKeepTheSettingOfTheOperations)
{
  for (const double noise : {0.0, 1.0e-5, -1.0e-5})
  {
    for (std::size_t hallNumber = 1; hallNumber <= 530; ++hallNumber)
    {
      const SKSpaceGroupSetting& setting = SKSpaceGroup(hallNumber).spaceGroupSetting();
      const std::optional<std::size_t> found = SKSpaceGroup::HallNumberFromSymmetryOperations(
          operationsOfSetting(hallNumber, noise), latticeOfSetting(setting.number()), setting.centring(),
          setting.number(), 1.0e-3);
      ASSERT_TRUE(found.has_value()) << hallNumber;
      if (*found != hallNumber)
      {
        // only another name for the same operations
        EXPECT_TRUE(SKSpaceGroup(*found).spaceGroupSetting().fullSeitzMatrices().operations ==
                    setting.fullSeitzMatrices().operations)
            << "Hall number " << hallNumber << " read as " << *found << " with noise " << noise;
      }
    }
  }
}

// One operation is P 1 only when it is the identity.
TEST(CifSymmetry, ASingleOperationThatIsNotTheIdentityIsNoGroup)
{
  const auto hallNumberOf = [](const char* operation)
  {
    return SKSpaceGroup::HallNumberFromSymmetryOperations({*SKSpaceGroup::parseCifSymmetryOperation(operation)},
                                                          latticeOfSetting(1), Centring::primitive, std::nullopt);
  };
  EXPECT_EQ(hallNumberOf("x,y,z"), std::optional<std::size_t>(1));
  EXPECT_EQ(hallNumberOf("x+1,y,z"), std::optional<std::size_t>(1));
  EXPECT_FALSE(hallNumberOf("-x,-y,-z").has_value());
  EXPECT_FALSE(hallNumberOf("x,-y,z+1/2").has_value());
  EXPECT_FALSE(hallNumberOf("x+1/2,y,z").has_value());
}

// The H-prefixed symbols of the hexagonal-axes settings are read without spaces as they are with them.
TEST(CifSymmetry, CompactHexagonalAxesSymbols)
{
  for (const auto& [compact, spaced] :
       {std::pair{"H-3m", "H -3 m"}, std::pair{"H3m", "H 3 m"}, std::pair{"H-3c", "H -3 c"}, std::pair{"H3c", "H 3 c"},
        std::pair{"H32", "H 3 2"}, std::pair{"H-3", "H -3"}, std::pair{"H3", "H 3"}})
  {
    const std::optional<std::size_t> found = SKSpaceGroup::HallNumberFromHMString(compact);
    ASSERT_TRUE(found.has_value()) << compact;
    EXPECT_EQ(found, SKSpaceGroup::HallNumberFromHMString(spaced)) << compact;
    EXPECT_EQ(SKSpaceGroup(*found).spaceGroupSetting().qualifier(), "H") << compact;
  }
}

// What is no symmetry operation is refused, and nothing is thrown.
TEST(CifSymmetry, MalformedSymmetryOperationsAreRefused)
{
  for (const char* expression : {"xx,y,z", "x1/2,y,z", "x,x,x", "x,y,z,", "99999999999/2+x,y,z", "1/99999999999+x,y,z",
                                 "x,y", "x,,z", "2x,y,z", "", ",,"})
  {
    std::optional<SKSeitzMatrix> operation;
    EXPECT_NO_THROW(operation = SKSpaceGroup::parseCifSymmetryOperation(expression)) << expression;
    EXPECT_FALSE(operation.has_value()) << expression;
  }
  for (const char* expression : {"x,y,z", "-x,1/2+y,-z+1/2", "x-y,x,z+0.33333", "'-y, x-y, z+2/3'", "y+3/4,x+1/4,-z+.75"})
  {
    EXPECT_TRUE(SKSpaceGroup::parseCifSymmetryOperation(expression).has_value()) << expression;
  }
}

// The refinement of a translation on no atoms is the translation (it was 0 / 0).
TEST(CifSymmetry, RefinedTranslationOnNoAtoms)
{
  const std::optional<double3> refined = SKSymmetryCell::refinedTranslation(
      double3(0.25, 0.5, 0.0), SKRotationMatrix::identity, {}, latticeOfSetting(200), false, 1.0e-2);
  ASSERT_TRUE(refined.has_value());
  EXPECT_DOUBLE_EQ(refined->x, 0.25);
  EXPECT_DOUBLE_EQ(refined->y, 0.5);
  EXPECT_DOUBLE_EQ(refined->z, 0.0);
}

// An atom that lies just inside the far face of the cell, on a site at the origin, is symmetrized to the origin
// and not to 1.
TEST(CifSymmetry, SymmetrizedPositionsStayInTheCell)
{
  SKIntegerSymmetryOperationSet operations = SKSpaceGroup(2).spaceGroupSetting().fullSeitzMatrices();  // P -1
  const auto symmetrized = operations.symmetrize(
      latticeOfSetting(1), {{double3(0.9999, 0.0, 0.0), 1, 1.0}, {double3(0.5, 0.9999, 0.5), 1, 1.0}}, 1.0e-2);
  ASSERT_EQ(symmetrized.size(), 2U);
  for (const auto& atom : symmetrized)
  {
    const double3 position = std::get<0>(atom);
    for (const double coordinate : {position.x, position.y, position.z})
    {
      EXPECT_GE(coordinate, 0.0);
      EXPECT_LT(coordinate, 1.0);
    }
  }
  EXPECT_DOUBLE_EQ(std::get<0>(symmetrized[0]).x, 0.0);
  EXPECT_DOUBLE_EQ(std::get<0>(symmetrized[1]).y, 0.0);
}

// Two sites closer than the precision give every operation twice, a translation apart that the precision allows;
// kept apart by a fixed fractional tolerance in a small cell, the two identities left no point group.
TEST(CifSymmetry, SplitSitesInASmallCellHaveASpaceGroup)
{
  const double3x3 lattice(double3(3.0, 0.0, 0.0), double3(0.0, 3.0, 0.0), double3(0.0, 0.0, 3.0));
  const std::vector<std::tuple<double3, std::size_t, double>> atoms = {{double3(0.0, 0.0, 0.0), 1, 0.5},
                                                                       {double3(0.007, 0.0, 0.0), 2, 0.5}};
  const auto found = SKSpaceGroup::findSpaceGroupAtPrecision(lattice, atoms, true, 0.05);
  EXPECT_TRUE(found.has_value());
}
