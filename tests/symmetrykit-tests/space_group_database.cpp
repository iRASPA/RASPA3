#include <gtest/gtest.h>

import std;

import int3;
import int3x3;
import double3;
import skdefinitions;
import skrotationmatrix;
import skseitzintegermatrix;
import skspacegroup;
import skspacegroupdatabase;
import skspacegroupsetting;

// The symmetry operations of every space-group setting are held in the database as a
// string, three characters per operation and one character per row of the Seitz matrix.
// Nothing in the reader validates that string: SKSeitzIntegerMatrix::SeitzMatrices takes
// its length divided by three as the number of operations and indexes the row table
// without a bounds check. A string of the wrong length therefore yields the wrong number
// of operations, and, because the grouping shifts, garbage for every operation after the
// damage. The atoms that follow are then at wrong positions rather than on top of each
// other, so the coincident-image merge in the CIF reader cannot repair or even notice it.
//
// The tests below check the invariants that make such a string detectable.

namespace
{
// Translations are stored with a denominator of 24. Composition can drive a component
// negative, and the built-in comparison uses the sign-preserving %, so reduce into
// [0, 24) before comparing.
int3 reducedTranslation(int3 translation)
{
  const auto wrap = [](int value) { return ((value % 24) + 24) % 24; };
  return int3(wrap(translation.x), wrap(translation.y), wrap(translation.z));
}

using OperationKey = std::array<int, 12>;

OperationKey keyOf(const SKSeitzIntegerMatrix& matrix)
{
  const int3x3& r = matrix.rotation.int3x3_m;
  const int3 t = reducedTranslation(matrix.translation);
  return OperationKey{r.m11, r.m12, r.m13, r.m21, r.m22, r.m23, r.m31, r.m32, r.m33, t.x, t.y, t.z};
}

// spaceGroupData[0] is a placeholder; the settings proper are 1 through 530.
constexpr std::size_t firstHallNumber = 1;
constexpr std::size_t lastHallNumber = 530;

std::string describe(std::size_t hallNumber)
{
  const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
  return std::format("Hall {} (space group {}, '{}', Hall symbol '{}', qualifier '{}')", hallNumber,
                     setting.number(), setting.HMString(), setting.HallString(), setting.qualifier());
}
}  // namespace

TEST(SpaceGroupDatabase, EncodedSeitzStringHoldsThreeCharactersPerOperation)
{
  for (std::size_t hallNumber = firstHallNumber; hallNumber <= lastHallNumber; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    EXPECT_EQ(setting.encodedSeitz().size(), 3 * setting.order()) << describe(hallNumber);
  }
}

TEST(SpaceGroupDatabase, EncodedSeitzStringUsesOnlyKnownRows)
{
  const std::size_t rows = SKSeitzIntegerMatrix::SeitzData.size();
  for (std::size_t hallNumber = firstHallNumber; hallNumber <= lastHallNumber; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    for (const char character : setting.encodedSeitz())
    {
      const int index = character - '0';
      EXPECT_GE(index, 0) << describe(hallNumber) << ", character '" << character << "'";
      EXPECT_LT(static_cast<std::size_t>(index), rows)
          << describe(hallNumber) << ", character '" << character << "'";
    }
  }
}

TEST(SpaceGroupDatabase, EverySettingExpandsToAGroupOfItsDeclaredOrder)
{
  for (std::size_t hallNumber = firstHallNumber; hallNumber <= lastHallNumber; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const std::string where = describe(hallNumber);

    const auto operations = setting.fullSeitzMatrices().operations;
    ASSERT_EQ(operations.size(), setting.order()) << where;

    std::vector<SKSeitzIntegerMatrix> elements(operations.begin(), operations.end());
    std::set<OperationKey> present;
    for (const SKSeitzIntegerMatrix& element : elements)
    {
      present.insert(keyOf(element));
    }
    ASSERT_EQ(present.size(), setting.order()) << where << ": operations are not distinct modulo a lattice vector";

    // Every element is a proper crystallographic rotation.
    for (const SKSeitzIntegerMatrix& element : elements)
    {
      SKRotationMatrix rotation = element.rotation;
      EXPECT_EQ(std::abs(rotation.determinant()), 1) << where << ": rotation with a determinant that is not +-1";
    }

    EXPECT_TRUE(present.contains(keyOf(SKSeitzIntegerMatrix(SKRotationMatrix::identity, int3(0, 0, 0)))))
        << where << ": the identity is missing";

    // Closed under composition. Report the first failure rather than one per pair.
    std::optional<std::pair<std::size_t, std::size_t>> notClosed;
    for (std::size_t i = 0; i < elements.size() && !notClosed; ++i)
    {
      for (std::size_t j = 0; j < elements.size(); ++j)
      {
        if (!present.contains(keyOf(elements[i] * elements[j])))
        {
          notClosed = std::make_pair(i, j);
          break;
        }
      }
    }
    EXPECT_FALSE(notClosed.has_value()) << where << ": the operations are not closed under composition";
  }
}

TEST(SpaceGroupDatabase, HexagonalSettingsOfRhombohedralGroupsCarryTheirCentringTranslations)
{
  // These seven are the only settings whose encoded string needs the row '-x+y+1/3',
  // which is a single backslash, so they are the ones an escaping slip in the C++ source
  // silently corrupts. Reading a CIF in any of them then gives a structure that is wrong
  // in both the atom count and the atom positions, with no diagnostic.
  const std::vector<std::pair<std::size_t, std::size_t>> hexagonalSettings = {
      {433, 9}, {436, 18}, {444, 18}, {450, 18}, {452, 18}, {458, 36}, {460, 36}};

  // (2/3, 1/3, 1/3) and (1/3, 2/3, 2/3), with a denominator of 24.
  const std::vector<int3> centringTranslations = {int3(16, 8, 8), int3(8, 16, 16)};

  for (const auto& [hallNumber, expectedOrder] : hexagonalSettings)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const std::string where = describe(hallNumber);

    ASSERT_EQ(setting.qualifier(), "H") << where;
    ASSERT_EQ(setting.centring(), Centring::r) << where;
    ASSERT_EQ(setting.order(), expectedOrder) << where;

    const auto operations = setting.fullSeitzMatrices().operations;
    ASSERT_EQ(operations.size(), expectedOrder) << where;

    std::set<OperationKey> present;
    for (const SKSeitzIntegerMatrix& element : operations)
    {
      present.insert(keyOf(element));
    }

    for (const int3& translation : centringTranslations)
    {
      EXPECT_TRUE(present.contains(keyOf(SKSeitzIntegerMatrix(SKRotationMatrix::identity, translation))))
          << where << ": the R-centring translation (" << translation.x << ", " << translation.y << ", "
          << translation.z << ")/24 is not among the operations";
    }
  }
}

TEST(SpaceGroupDatabase, ExpandingAGeneralPositionGivesOneImagePerOperation)
{
  // A point in general position has as many distinct images as the group has operations.
  // Rhombohedral R -3 m, the setting DDR is given in, is the case that matters most here.
  const double3 generalPosition = double3(0.7267, 0.0511, 0.07);

  for (std::size_t hallNumber = firstHallNumber; hallNumber <= lastHallNumber; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const std::vector<double3> positions = SKSpaceGroup(hallNumber).listOfSymmetricPositions(generalPosition);
    EXPECT_EQ(positions.size(), setting.order()) << describe(hallNumber);
  }
}
