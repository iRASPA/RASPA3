#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <optional>
#include <print>
#include <random>
#include <string>
#include <tuple>
#include <vector>

import double3;
import double3x3;
import randomnumbers;

import skposcarlegacyparser;
import sksymmetrycell;
import skspacegroup;

TEST(TransformationMatrix, Triclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {"spglibtestdata/triclinic/POSCAR-001",
                                             "spglibtestdata/triclinic/POSCAR-002"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Monoclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {
      "spglibtestdata/monoclinic/POSCAR-003",   "spglibtestdata/monoclinic/POSCAR-004",
      "spglibtestdata/monoclinic/POSCAR-004-2", "spglibtestdata/monoclinic/POSCAR-005",
      "spglibtestdata/monoclinic/POSCAR-005-2", "spglibtestdata/monoclinic/POSCAR-006",
      "spglibtestdata/monoclinic/POSCAR-006-2", "spglibtestdata/monoclinic/POSCAR-007",
      "spglibtestdata/monoclinic/POSCAR-007-2", "spglibtestdata/monoclinic/POSCAR-008",
      "spglibtestdata/monoclinic/POSCAR-008-2", "spglibtestdata/monoclinic/POSCAR-009",
      "spglibtestdata/monoclinic/POSCAR-009-2", "spglibtestdata/monoclinic/POSCAR-010",
      "spglibtestdata/monoclinic/POSCAR-010-2", "spglibtestdata/monoclinic/POSCAR-011",
      "spglibtestdata/monoclinic/POSCAR-011-2", "spglibtestdata/monoclinic/POSCAR-012",
      "spglibtestdata/monoclinic/POSCAR-012-2", "spglibtestdata/monoclinic/POSCAR-012-3",
      "spglibtestdata/monoclinic/POSCAR-013",   "spglibtestdata/monoclinic/POSCAR-013-2",
      "spglibtestdata/monoclinic/POSCAR-013-3", "spglibtestdata/monoclinic/POSCAR-014",
      "spglibtestdata/monoclinic/POSCAR-014-2", "spglibtestdata/monoclinic/POSCAR-015",
      "spglibtestdata/monoclinic/POSCAR-015-2", "spglibtestdata/monoclinic/POSCAR-015-3"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Orthorhombic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {
      "spglibtestdata/orthorhombic/POSCAR-016",   "spglibtestdata/orthorhombic/POSCAR-016-2",
      "spglibtestdata/orthorhombic/POSCAR-017-2", "spglibtestdata/orthorhombic/POSCAR-018",
      "spglibtestdata/orthorhombic/POSCAR-018-2", "spglibtestdata/orthorhombic/POSCAR-019",
      "spglibtestdata/orthorhombic/POSCAR-019-2", "spglibtestdata/orthorhombic/POSCAR-020",
      "spglibtestdata/orthorhombic/POSCAR-021",   "spglibtestdata/orthorhombic/POSCAR-021-2",
      "spglibtestdata/orthorhombic/POSCAR-022",   "spglibtestdata/orthorhombic/POSCAR-023",
      "spglibtestdata/orthorhombic/POSCAR-023-2", "spglibtestdata/orthorhombic/POSCAR-024",
      "spglibtestdata/orthorhombic/POSCAR-024-2", "spglibtestdata/orthorhombic/POSCAR-025",
      "spglibtestdata/orthorhombic/POSCAR-025-2", "spglibtestdata/orthorhombic/POSCAR-026",
      "spglibtestdata/orthorhombic/POSCAR-026-2", "spglibtestdata/orthorhombic/POSCAR-027",
      "spglibtestdata/orthorhombic/POSCAR-027-2", "spglibtestdata/orthorhombic/POSCAR-028",
      "spglibtestdata/orthorhombic/POSCAR-028-2", "spglibtestdata/orthorhombic/POSCAR-029",
      "spglibtestdata/orthorhombic/POSCAR-029-2", "spglibtestdata/orthorhombic/POSCAR-030",
      "spglibtestdata/orthorhombic/POSCAR-030-2", "spglibtestdata/orthorhombic/POSCAR-031",
      "spglibtestdata/orthorhombic/POSCAR-031-2", "spglibtestdata/orthorhombic/POSCAR-032",
      "spglibtestdata/orthorhombic/POSCAR-032-2", "spglibtestdata/orthorhombic/POSCAR-033",
      "spglibtestdata/orthorhombic/POSCAR-033-2", "spglibtestdata/orthorhombic/POSCAR-033-3",
      "spglibtestdata/orthorhombic/POSCAR-034",   "spglibtestdata/orthorhombic/POSCAR-034-2",
      "spglibtestdata/orthorhombic/POSCAR-035",   "spglibtestdata/orthorhombic/POSCAR-035-2",
      "spglibtestdata/orthorhombic/POSCAR-036",   "spglibtestdata/orthorhombic/POSCAR-036-2",
      "spglibtestdata/orthorhombic/POSCAR-037",   "spglibtestdata/orthorhombic/POSCAR-037-2",
      "spglibtestdata/orthorhombic/POSCAR-038",   "spglibtestdata/orthorhombic/POSCAR-038-2",
      "spglibtestdata/orthorhombic/POSCAR-039",   "spglibtestdata/orthorhombic/POSCAR-039-2",
      "spglibtestdata/orthorhombic/POSCAR-040",   "spglibtestdata/orthorhombic/POSCAR-040-2",
      "spglibtestdata/orthorhombic/POSCAR-041",   "spglibtestdata/orthorhombic/POSCAR-041-2",
      "spglibtestdata/orthorhombic/POSCAR-042",   "spglibtestdata/orthorhombic/POSCAR-043",
      "spglibtestdata/orthorhombic/POSCAR-043-2", "spglibtestdata/orthorhombic/POSCAR-044",
      "spglibtestdata/orthorhombic/POSCAR-044-2", "spglibtestdata/orthorhombic/POSCAR-045",
      "spglibtestdata/orthorhombic/POSCAR-045-2", "spglibtestdata/orthorhombic/POSCAR-046",
      "spglibtestdata/orthorhombic/POSCAR-046-2", "spglibtestdata/orthorhombic/POSCAR-047",
      "spglibtestdata/orthorhombic/POSCAR-047-2", "spglibtestdata/orthorhombic/POSCAR-048",
      "spglibtestdata/orthorhombic/POSCAR-048-2", "spglibtestdata/orthorhombic/POSCAR-049",
      "spglibtestdata/orthorhombic/POSCAR-049-2", "spglibtestdata/orthorhombic/POSCAR-050",
      "spglibtestdata/orthorhombic/POSCAR-050-2", "spglibtestdata/orthorhombic/POSCAR-051",
      "spglibtestdata/orthorhombic/POSCAR-051-2", "spglibtestdata/orthorhombic/POSCAR-051-3",
      "spglibtestdata/orthorhombic/POSCAR-052",   "spglibtestdata/orthorhombic/POSCAR-052-2",
      "spglibtestdata/orthorhombic/POSCAR-053",   "spglibtestdata/orthorhombic/POSCAR-053-2",
      "spglibtestdata/orthorhombic/POSCAR-054",   "spglibtestdata/orthorhombic/POSCAR-054-2",
      "spglibtestdata/orthorhombic/POSCAR-055",   "spglibtestdata/orthorhombic/POSCAR-055-2",
      "spglibtestdata/orthorhombic/POSCAR-056",   "spglibtestdata/orthorhombic/POSCAR-056-2",
      "spglibtestdata/orthorhombic/POSCAR-057",   "spglibtestdata/orthorhombic/POSCAR-057-2",
      "spglibtestdata/orthorhombic/POSCAR-058",   "spglibtestdata/orthorhombic/POSCAR-058-2",
      "spglibtestdata/orthorhombic/POSCAR-058-3", "spglibtestdata/orthorhombic/POSCAR-059",
      "spglibtestdata/orthorhombic/POSCAR-059-2", "spglibtestdata/orthorhombic/POSCAR-060",
      "spglibtestdata/orthorhombic/POSCAR-060-2", "spglibtestdata/orthorhombic/POSCAR-060-3",
      "spglibtestdata/orthorhombic/POSCAR-061",   "spglibtestdata/orthorhombic/POSCAR-061-2",
      "spglibtestdata/orthorhombic/POSCAR-062",   "spglibtestdata/orthorhombic/POSCAR-062-2",
      "spglibtestdata/orthorhombic/POSCAR-063",   "spglibtestdata/orthorhombic/POSCAR-063-2",
      "spglibtestdata/orthorhombic/POSCAR-063-3", "spglibtestdata/orthorhombic/POSCAR-064",
      "spglibtestdata/orthorhombic/POSCAR-064-2", "spglibtestdata/orthorhombic/POSCAR-064-3",
      "spglibtestdata/orthorhombic/POSCAR-065",   "spglibtestdata/orthorhombic/POSCAR-065-2",
      "spglibtestdata/orthorhombic/POSCAR-065-3", "spglibtestdata/orthorhombic/POSCAR-066",
      "spglibtestdata/orthorhombic/POSCAR-066-2", "spglibtestdata/orthorhombic/POSCAR-067",
      "spglibtestdata/orthorhombic/POSCAR-067-2", "spglibtestdata/orthorhombic/POSCAR-067-3",
      "spglibtestdata/orthorhombic/POSCAR-068",   "spglibtestdata/orthorhombic/POSCAR-068-2",
      "spglibtestdata/orthorhombic/POSCAR-069",   "spglibtestdata/orthorhombic/POSCAR-069-2",
      "spglibtestdata/orthorhombic/POSCAR-070",   "spglibtestdata/orthorhombic/POSCAR-070-2",
      "spglibtestdata/orthorhombic/POSCAR-071",   "spglibtestdata/orthorhombic/POSCAR-071-2",
      "spglibtestdata/orthorhombic/POSCAR-072",   "spglibtestdata/orthorhombic/POSCAR-072-2",
      "spglibtestdata/orthorhombic/POSCAR-073",   "spglibtestdata/orthorhombic/POSCAR-073-2",
      "spglibtestdata/orthorhombic/POSCAR-074",   "spglibtestdata/orthorhombic/POSCAR-074-2"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Tetragonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {
      "spglibtestdata/tetragonal/POSCAR-075",   "spglibtestdata/tetragonal/POSCAR-075-2",
      "spglibtestdata/tetragonal/POSCAR-076",   "spglibtestdata/tetragonal/POSCAR-076-2",
      "spglibtestdata/tetragonal/POSCAR-077",   "spglibtestdata/tetragonal/POSCAR-077-2",
      "spglibtestdata/tetragonal/POSCAR-077-3", "spglibtestdata/tetragonal/POSCAR-078",
      "spglibtestdata/tetragonal/POSCAR-078-2", "spglibtestdata/tetragonal/POSCAR-079",
      "spglibtestdata/tetragonal/POSCAR-079-2", "spglibtestdata/tetragonal/POSCAR-080",
      "spglibtestdata/tetragonal/POSCAR-080-2", "spglibtestdata/tetragonal/POSCAR-081",
      "spglibtestdata/tetragonal/POSCAR-081-2", "spglibtestdata/tetragonal/POSCAR-082",
      "spglibtestdata/tetragonal/POSCAR-082-2", "spglibtestdata/tetragonal/POSCAR-083",
      "spglibtestdata/tetragonal/POSCAR-083-2", "spglibtestdata/tetragonal/POSCAR-083-3",
      "spglibtestdata/tetragonal/POSCAR-084",   "spglibtestdata/tetragonal/POSCAR-084-2",
      "spglibtestdata/tetragonal/POSCAR-085",   "spglibtestdata/tetragonal/POSCAR-085-2",
      "spglibtestdata/tetragonal/POSCAR-086",   "spglibtestdata/tetragonal/POSCAR-086-2",
      "spglibtestdata/tetragonal/POSCAR-087",   "spglibtestdata/tetragonal/POSCAR-087-2",
      "spglibtestdata/tetragonal/POSCAR-088",   "spglibtestdata/tetragonal/POSCAR-088-2",
      "spglibtestdata/tetragonal/POSCAR-090",   "spglibtestdata/tetragonal/POSCAR-090-2",
      "spglibtestdata/tetragonal/POSCAR-091",   "spglibtestdata/tetragonal/POSCAR-091-2",
      "spglibtestdata/tetragonal/POSCAR-092",   "spglibtestdata/tetragonal/POSCAR-092-2",
      "spglibtestdata/tetragonal/POSCAR-092-3", "spglibtestdata/tetragonal/POSCAR-094",
      "spglibtestdata/tetragonal/POSCAR-094-2", "spglibtestdata/tetragonal/POSCAR-094-3",
      "spglibtestdata/tetragonal/POSCAR-095",   "spglibtestdata/tetragonal/POSCAR-095-2",
      "spglibtestdata/tetragonal/POSCAR-096",   "spglibtestdata/tetragonal/POSCAR-096-2",
      "spglibtestdata/tetragonal/POSCAR-097",   "spglibtestdata/tetragonal/POSCAR-097-2",
      "spglibtestdata/tetragonal/POSCAR-098",   "spglibtestdata/tetragonal/POSCAR-098-2",
      "spglibtestdata/tetragonal/POSCAR-099",   "spglibtestdata/tetragonal/POSCAR-099-2",
      "spglibtestdata/tetragonal/POSCAR-100",   "spglibtestdata/tetragonal/POSCAR-100-2",
      "spglibtestdata/tetragonal/POSCAR-102",   "spglibtestdata/tetragonal/POSCAR-102-2",
      "spglibtestdata/tetragonal/POSCAR-103",   "spglibtestdata/tetragonal/POSCAR-103-2",
      "spglibtestdata/tetragonal/POSCAR-104",   "spglibtestdata/tetragonal/POSCAR-104-2",
      "spglibtestdata/tetragonal/POSCAR-105",   "spglibtestdata/tetragonal/POSCAR-105-2",
      "spglibtestdata/tetragonal/POSCAR-106",   "spglibtestdata/tetragonal/POSCAR-107",
      "spglibtestdata/tetragonal/POSCAR-107-2", "spglibtestdata/tetragonal/POSCAR-107-3",
      "spglibtestdata/tetragonal/POSCAR-108",   "spglibtestdata/tetragonal/POSCAR-108-2",
      "spglibtestdata/tetragonal/POSCAR-109",   "spglibtestdata/tetragonal/POSCAR-109-2",
      "spglibtestdata/tetragonal/POSCAR-110",   "spglibtestdata/tetragonal/POSCAR-110-2",
      "spglibtestdata/tetragonal/POSCAR-111",   "spglibtestdata/tetragonal/POSCAR-111-2",
      "spglibtestdata/tetragonal/POSCAR-112",   "spglibtestdata/tetragonal/POSCAR-112-2",
      "spglibtestdata/tetragonal/POSCAR-113",   "spglibtestdata/tetragonal/POSCAR-113-2",
      "spglibtestdata/tetragonal/POSCAR-114",   "spglibtestdata/tetragonal/POSCAR-114-2",
      "spglibtestdata/tetragonal/POSCAR-115",   "spglibtestdata/tetragonal/POSCAR-115-2",
      "spglibtestdata/tetragonal/POSCAR-115-3", "spglibtestdata/tetragonal/POSCAR-115-4",
      "spglibtestdata/tetragonal/POSCAR-115-5", "spglibtestdata/tetragonal/POSCAR-116",
      "spglibtestdata/tetragonal/POSCAR-116-2", "spglibtestdata/tetragonal/POSCAR-117",
      "spglibtestdata/tetragonal/POSCAR-117-2", "spglibtestdata/tetragonal/POSCAR-118",
      "spglibtestdata/tetragonal/POSCAR-118-2", "spglibtestdata/tetragonal/POSCAR-119",
      "spglibtestdata/tetragonal/POSCAR-119-2", "spglibtestdata/tetragonal/POSCAR-120",
      "spglibtestdata/tetragonal/POSCAR-120-2", "spglibtestdata/tetragonal/POSCAR-121",
      "spglibtestdata/tetragonal/POSCAR-121-2", "spglibtestdata/tetragonal/POSCAR-122",
      "spglibtestdata/tetragonal/POSCAR-122-2", "spglibtestdata/tetragonal/POSCAR-122-3",
      "spglibtestdata/tetragonal/POSCAR-123",   "spglibtestdata/tetragonal/POSCAR-123-2",
      "spglibtestdata/tetragonal/POSCAR-123-3", "spglibtestdata/tetragonal/POSCAR-124",
      "spglibtestdata/tetragonal/POSCAR-124-2", "spglibtestdata/tetragonal/POSCAR-125",
      "spglibtestdata/tetragonal/POSCAR-125-2", "spglibtestdata/tetragonal/POSCAR-126",
      "spglibtestdata/tetragonal/POSCAR-126-2", "spglibtestdata/tetragonal/POSCAR-127",
      "spglibtestdata/tetragonal/POSCAR-127-2", "spglibtestdata/tetragonal/POSCAR-128",
      "spglibtestdata/tetragonal/POSCAR-128-2", "spglibtestdata/tetragonal/POSCAR-129",
      "spglibtestdata/tetragonal/POSCAR-129-2", "spglibtestdata/tetragonal/POSCAR-129-3",
      "spglibtestdata/tetragonal/POSCAR-130",   "spglibtestdata/tetragonal/POSCAR-130-2",
      "spglibtestdata/tetragonal/POSCAR-131",   "spglibtestdata/tetragonal/POSCAR-131-2",
      "spglibtestdata/tetragonal/POSCAR-132",   "spglibtestdata/tetragonal/POSCAR-132-2",
      "spglibtestdata/tetragonal/POSCAR-133",   "spglibtestdata/tetragonal/POSCAR-133-2",
      "spglibtestdata/tetragonal/POSCAR-134",   "spglibtestdata/tetragonal/POSCAR-134-2",
      "spglibtestdata/tetragonal/POSCAR-135",   "spglibtestdata/tetragonal/POSCAR-135-2",
      "spglibtestdata/tetragonal/POSCAR-136",   "spglibtestdata/tetragonal/POSCAR-136-2",
      "spglibtestdata/tetragonal/POSCAR-136-3", "spglibtestdata/tetragonal/POSCAR-136-4",
      "spglibtestdata/tetragonal/POSCAR-136-5", "spglibtestdata/tetragonal/POSCAR-137",
      "spglibtestdata/tetragonal/POSCAR-137-2", "spglibtestdata/tetragonal/POSCAR-137-3",
      "spglibtestdata/tetragonal/POSCAR-138",   "spglibtestdata/tetragonal/POSCAR-138-2",
      "spglibtestdata/tetragonal/POSCAR-139",   "spglibtestdata/tetragonal/POSCAR-139-2",
      "spglibtestdata/tetragonal/POSCAR-140",   "spglibtestdata/tetragonal/POSCAR-140-2",
      "spglibtestdata/tetragonal/POSCAR-141",   "spglibtestdata/tetragonal/POSCAR-141-2",
      "spglibtestdata/tetragonal/POSCAR-142",   "spglibtestdata/tetragonal/POSCAR-142-2",
      "spglibtestdata/tetragonal/POSCAR-142-3"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Trigonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {
      "spglibtestdata/trigonal/POSCAR-143",   "spglibtestdata/trigonal/POSCAR-143-2",
      "spglibtestdata/trigonal/POSCAR-144",   "spglibtestdata/trigonal/POSCAR-144-2",
      "spglibtestdata/trigonal/POSCAR-145",   "spglibtestdata/trigonal/POSCAR-145-2",
      "spglibtestdata/trigonal/POSCAR-146",   "spglibtestdata/trigonal/POSCAR-146-2",
      "spglibtestdata/trigonal/POSCAR-147",   "spglibtestdata/trigonal/POSCAR-147-2",
      "spglibtestdata/trigonal/POSCAR-148",   "spglibtestdata/trigonal/POSCAR-148-2",
      "spglibtestdata/trigonal/POSCAR-149",   "spglibtestdata/trigonal/POSCAR-149-2",
      "spglibtestdata/trigonal/POSCAR-150",   "spglibtestdata/trigonal/POSCAR-150-2",
      "spglibtestdata/trigonal/POSCAR-151",   "spglibtestdata/trigonal/POSCAR-151-2",
      "spglibtestdata/trigonal/POSCAR-152",   "spglibtestdata/trigonal/POSCAR-152-2",
      "spglibtestdata/trigonal/POSCAR-153",   "spglibtestdata/trigonal/POSCAR-154",
      "spglibtestdata/trigonal/POSCAR-154-2", "spglibtestdata/trigonal/POSCAR-154-3",
      "spglibtestdata/trigonal/POSCAR-155",   "spglibtestdata/trigonal/POSCAR-155-2",
      "spglibtestdata/trigonal/POSCAR-156",   "spglibtestdata/trigonal/POSCAR-156-2",
      "spglibtestdata/trigonal/POSCAR-157",   "spglibtestdata/trigonal/POSCAR-157-2",
      "spglibtestdata/trigonal/POSCAR-158",   "spglibtestdata/trigonal/POSCAR-158-2",
      "spglibtestdata/trigonal/POSCAR-159",   "spglibtestdata/trigonal/POSCAR-159-2",
      "spglibtestdata/trigonal/POSCAR-160",   "spglibtestdata/trigonal/POSCAR-160-2",
      "spglibtestdata/trigonal/POSCAR-161",   "spglibtestdata/trigonal/POSCAR-161-2",
      "spglibtestdata/trigonal/POSCAR-162",   "spglibtestdata/trigonal/POSCAR-162-2",
      "spglibtestdata/trigonal/POSCAR-163",   "spglibtestdata/trigonal/POSCAR-163-2",
      "spglibtestdata/trigonal/POSCAR-164",   "spglibtestdata/trigonal/POSCAR-164-2",
      "spglibtestdata/trigonal/POSCAR-165",   "spglibtestdata/trigonal/POSCAR-165-2",
      "spglibtestdata/trigonal/POSCAR-166",   "spglibtestdata/trigonal/POSCAR-166-2",
      "spglibtestdata/trigonal/POSCAR-167",   "spglibtestdata/trigonal/POSCAR-167-2",
      "spglibtestdata/trigonal/POSCAR-167-3"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Hexagonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {
      "spglibtestdata/hexagonal/POSCAR-168",   "spglibtestdata/hexagonal/POSCAR-169",
      "spglibtestdata/hexagonal/POSCAR-169-2", "spglibtestdata/hexagonal/POSCAR-170",
      "spglibtestdata/hexagonal/POSCAR-170-2", "spglibtestdata/hexagonal/POSCAR-171",
      "spglibtestdata/hexagonal/POSCAR-171-2", "spglibtestdata/hexagonal/POSCAR-172",
      "spglibtestdata/hexagonal/POSCAR-173",   "spglibtestdata/hexagonal/POSCAR-173-2",
      "spglibtestdata/hexagonal/POSCAR-174",   "spglibtestdata/hexagonal/POSCAR-174-2",
      "spglibtestdata/hexagonal/POSCAR-175",   "spglibtestdata/hexagonal/POSCAR-175-2",
      "spglibtestdata/hexagonal/POSCAR-176",   "spglibtestdata/hexagonal/POSCAR-176-2",
      "spglibtestdata/hexagonal/POSCAR-177",   "spglibtestdata/hexagonal/POSCAR-179",
      "spglibtestdata/hexagonal/POSCAR-179-2", "spglibtestdata/hexagonal/POSCAR-180",
      "spglibtestdata/hexagonal/POSCAR-180-2", "spglibtestdata/hexagonal/POSCAR-181",
      "spglibtestdata/hexagonal/POSCAR-181-2", "spglibtestdata/hexagonal/POSCAR-182",
      "spglibtestdata/hexagonal/POSCAR-182-2", "spglibtestdata/hexagonal/POSCAR-183",
      "spglibtestdata/hexagonal/POSCAR-183-2", "spglibtestdata/hexagonal/POSCAR-184",
      "spglibtestdata/hexagonal/POSCAR-184-2", "spglibtestdata/hexagonal/POSCAR-185",
      "spglibtestdata/hexagonal/POSCAR-185-2", "spglibtestdata/hexagonal/POSCAR-186",
      "spglibtestdata/hexagonal/POSCAR-186-2", "spglibtestdata/hexagonal/POSCAR-187",
      "spglibtestdata/hexagonal/POSCAR-187-2", "spglibtestdata/hexagonal/POSCAR-188",
      "spglibtestdata/hexagonal/POSCAR-188-2", "spglibtestdata/hexagonal/POSCAR-189",
      "spglibtestdata/hexagonal/POSCAR-189-2", "spglibtestdata/hexagonal/POSCAR-190",
      "spglibtestdata/hexagonal/POSCAR-190-2", "spglibtestdata/hexagonal/POSCAR-191",
      "spglibtestdata/hexagonal/POSCAR-191-2", "spglibtestdata/hexagonal/POSCAR-192",
      "spglibtestdata/hexagonal/POSCAR-192-2", "spglibtestdata/hexagonal/POSCAR-193",
      "spglibtestdata/hexagonal/POSCAR-193-2", "spglibtestdata/hexagonal/POSCAR-194",
      "spglibtestdata/hexagonal/POSCAR-194-2"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Cubic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {
      "spglibtestdata/cubic/POSCAR-195",   "spglibtestdata/cubic/POSCAR-195-2", "spglibtestdata/cubic/POSCAR-196",
      "spglibtestdata/cubic/POSCAR-196-2", "spglibtestdata/cubic/POSCAR-197",   "spglibtestdata/cubic/POSCAR-197-2",
      "spglibtestdata/cubic/POSCAR-198",   "spglibtestdata/cubic/POSCAR-198-2", "spglibtestdata/cubic/POSCAR-199",
      "spglibtestdata/cubic/POSCAR-199-2", "spglibtestdata/cubic/POSCAR-200",   "spglibtestdata/cubic/POSCAR-200-2",
      "spglibtestdata/cubic/POSCAR-205",   "spglibtestdata/cubic/POSCAR-205-3", "spglibtestdata/cubic/POSCAR-206",
      "spglibtestdata/cubic/POSCAR-206-2", "spglibtestdata/cubic/POSCAR-207",   "spglibtestdata/cubic/POSCAR-208",
      "spglibtestdata/cubic/POSCAR-208-2", "spglibtestdata/cubic/POSCAR-209",   "spglibtestdata/cubic/POSCAR-210",
      "spglibtestdata/cubic/POSCAR-210-2", "spglibtestdata/cubic/POSCAR-211",   "spglibtestdata/cubic/POSCAR-212",
      "spglibtestdata/cubic/POSCAR-212-2", "spglibtestdata/cubic/POSCAR-213",   "spglibtestdata/cubic/POSCAR-213-2",
      "spglibtestdata/cubic/POSCAR-214",   "spglibtestdata/cubic/POSCAR-214-2", "spglibtestdata/cubic/POSCAR-215",
      "spglibtestdata/cubic/POSCAR-215-2", "spglibtestdata/cubic/POSCAR-216",   "spglibtestdata/cubic/POSCAR-216-2",
      "spglibtestdata/cubic/POSCAR-217",   "spglibtestdata/cubic/POSCAR-217-2", "spglibtestdata/cubic/POSCAR-218",
      "spglibtestdata/cubic/POSCAR-218-2", "spglibtestdata/cubic/POSCAR-219",   "spglibtestdata/cubic/POSCAR-219-2",
      "spglibtestdata/cubic/POSCAR-220",   "spglibtestdata/cubic/POSCAR-220-2", "spglibtestdata/cubic/POSCAR-221",
      "spglibtestdata/cubic/POSCAR-221-2", "spglibtestdata/cubic/POSCAR-222",   "spglibtestdata/cubic/POSCAR-222-2",
      "spglibtestdata/cubic/POSCAR-223",   "spglibtestdata/cubic/POSCAR-223-2", "spglibtestdata/cubic/POSCAR-224",
      "spglibtestdata/cubic/POSCAR-224-2", "spglibtestdata/cubic/POSCAR-225",   "spglibtestdata/cubic/POSCAR-225-2",
      "spglibtestdata/cubic/POSCAR-226",   "spglibtestdata/cubic/POSCAR-226-2", "spglibtestdata/cubic/POSCAR-227",
      "spglibtestdata/cubic/POSCAR-227-2", "spglibtestdata/cubic/POSCAR-228",   "spglibtestdata/cubic/POSCAR-228-2",
      "spglibtestdata/cubic/POSCAR-229",   "spglibtestdata/cubic/POSCAR-229-2", "spglibtestdata/cubic/POSCAR-230",
      "spglibtestdata/cubic/POSCAR-230-2", "spglibtestdata/cubic/POSCAR-230-3", "spglibtestdata/cubic/POSCAR-230-4"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(TransformationMatrix, Virtual)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::vector<std::string> testData = {"spglibtestdata/virtual_structure/POSCAR-1-221-33",
                                             "spglibtestdata/virtual_structure/POSCAR-1-222-33",
                                             "spglibtestdata/virtual_structure/POSCAR-1-223-33",
                                             "spglibtestdata/virtual_structure/POSCAR-1-224-33",
                                             "spglibtestdata/virtual_structure/POSCAR-1-227-73",
                                             "spglibtestdata/virtual_structure/POSCAR-1-227-93",
                                             "spglibtestdata/virtual_structure/POSCAR-1-227-99",
                                             "spglibtestdata/virtual_structure/POSCAR-1-230-conv-56",
                                             "spglibtestdata/virtual_structure/POSCAR-1-230-prim-33",
                                             "spglibtestdata/virtual_structure/POSCAR-1-bcc-33",
                                             "spglibtestdata/virtual_structure/POSCAR-10-221-18",
                                             "spglibtestdata/virtual_structure/POSCAR-10-223-18",
                                             "spglibtestdata/virtual_structure/POSCAR-10-227-50",
                                             "spglibtestdata/virtual_structure/POSCAR-102-224-13",
                                             "spglibtestdata/virtual_structure/POSCAR-104-222-13",
                                             "spglibtestdata/virtual_structure/POSCAR-105-223-13",
                                             "spglibtestdata/virtual_structure/POSCAR-109-227-13",
                                             "spglibtestdata/virtual_structure/POSCAR-11-227-48",
                                             "spglibtestdata/virtual_structure/POSCAR-110-230-conv-15",
                                             "spglibtestdata/virtual_structure/POSCAR-110-230-prim-13",
                                             "spglibtestdata/virtual_structure/POSCAR-111-221-11",
                                             "spglibtestdata/virtual_structure/POSCAR-111-224-11",
                                             "spglibtestdata/virtual_structure/POSCAR-111-227-66",
                                             "spglibtestdata/virtual_structure/POSCAR-112-222-11",
                                             "spglibtestdata/virtual_structure/POSCAR-112-223-11",
                                             "spglibtestdata/virtual_structure/POSCAR-113-227-68",
                                             "spglibtestdata/virtual_structure/POSCAR-115-221-14",
                                             "spglibtestdata/virtual_structure/POSCAR-115-223-14",
                                             "spglibtestdata/virtual_structure/POSCAR-115-227-33",
                                             "spglibtestdata/virtual_structure/POSCAR-116-230-conv-34",
                                             "spglibtestdata/virtual_structure/POSCAR-117-230-conv-33",
                                             "spglibtestdata/virtual_structure/POSCAR-118-222-14",
                                             "spglibtestdata/virtual_structure/POSCAR-118-224-14",
                                             "spglibtestdata/virtual_structure/POSCAR-12-221-19",
                                             "spglibtestdata/virtual_structure/POSCAR-12-224-19",
                                             "spglibtestdata/virtual_structure/POSCAR-12-227-21",
                                             "spglibtestdata/virtual_structure/POSCAR-12-227-83",
                                             "spglibtestdata/virtual_structure/POSCAR-120-230-conv-16",
                                             "spglibtestdata/virtual_structure/POSCAR-120-230-prim-14",
                                             "spglibtestdata/virtual_structure/POSCAR-122-230-conv-13",
                                             "spglibtestdata/virtual_structure/POSCAR-122-230-prim-11",
                                             "spglibtestdata/virtual_structure/POSCAR-123-221-05",
                                             "spglibtestdata/virtual_structure/POSCAR-126-222-05",
                                             "spglibtestdata/virtual_structure/POSCAR-13-222-18",
                                             "spglibtestdata/virtual_structure/POSCAR-13-224-18",
                                             "spglibtestdata/virtual_structure/POSCAR-13-227-49",
                                             "spglibtestdata/virtual_structure/POSCAR-13-230-conv-44",
                                             "spglibtestdata/virtual_structure/POSCAR-131-223-05",
                                             "spglibtestdata/virtual_structure/POSCAR-134-224-05",
                                             "spglibtestdata/virtual_structure/POSCAR-14-227-47",
                                             "spglibtestdata/virtual_structure/POSCAR-14-227-51",
                                             "spglibtestdata/virtual_structure/POSCAR-14-230-conv-45",
                                             "spglibtestdata/virtual_structure/POSCAR-142-230-conv-05",
                                             "spglibtestdata/virtual_structure/POSCAR-142-230-prim-05",
                                             "spglibtestdata/virtual_structure/POSCAR-146-221-27",
                                             "spglibtestdata/virtual_structure/POSCAR-146-222-27",
                                             "spglibtestdata/virtual_structure/POSCAR-146-223-27",
                                             "spglibtestdata/virtual_structure/POSCAR-146-224-27",
                                             "spglibtestdata/virtual_structure/POSCAR-146-227-92",
                                             "spglibtestdata/virtual_structure/POSCAR-146-230-conv-36",
                                             "spglibtestdata/virtual_structure/POSCAR-146-230-conv-55",
                                             "spglibtestdata/virtual_structure/POSCAR-146-230-prim-27",
                                             "spglibtestdata/virtual_structure/POSCAR-146-bcc-27",
                                             "spglibtestdata/virtual_structure/POSCAR-148-221-15",
                                             "spglibtestdata/virtual_structure/POSCAR-148-222-15",
                                             "spglibtestdata/virtual_structure/POSCAR-148-223-15",
                                             "spglibtestdata/virtual_structure/POSCAR-148-224-15",
                                             "spglibtestdata/virtual_structure/POSCAR-148-227-70",
                                             "spglibtestdata/virtual_structure/POSCAR-148-230-conv-17",
                                             "spglibtestdata/virtual_structure/POSCAR-148-230-conv-37",
                                             "spglibtestdata/virtual_structure/POSCAR-148-230-prim-15",
                                             "spglibtestdata/virtual_structure/POSCAR-148-bcc-15",
                                             "spglibtestdata/virtual_structure/POSCAR-15-222-19",
                                             "spglibtestdata/virtual_structure/POSCAR-15-223-19",
                                             "spglibtestdata/virtual_structure/POSCAR-15-230-conv-21",
                                             "spglibtestdata/virtual_structure/POSCAR-15-230-conv-22",
                                             "spglibtestdata/virtual_structure/POSCAR-15-230-prim-18",
                                             "spglibtestdata/virtual_structure/POSCAR-15-230-prim-19",
                                             "spglibtestdata/virtual_structure/POSCAR-15-bcc-18",
                                             "spglibtestdata/virtual_structure/POSCAR-15-bcc-19",
                                             "spglibtestdata/virtual_structure/POSCAR-155-221-17",
                                             "spglibtestdata/virtual_structure/POSCAR-155-222-17",
                                             "spglibtestdata/virtual_structure/POSCAR-155-223-17",
                                             "spglibtestdata/virtual_structure/POSCAR-155-224-17",
                                             "spglibtestdata/virtual_structure/POSCAR-155-227-72",
                                             "spglibtestdata/virtual_structure/POSCAR-155-230-conv-19",
                                             "spglibtestdata/virtual_structure/POSCAR-155-230-conv-38",
                                             "spglibtestdata/virtual_structure/POSCAR-155-230-prim-17",
                                             "spglibtestdata/virtual_structure/POSCAR-155-bcc-17",
                                             "spglibtestdata/virtual_structure/POSCAR-16-221-20",
                                             "spglibtestdata/virtual_structure/POSCAR-16-222-20",
                                             "spglibtestdata/virtual_structure/POSCAR-16-223-20",
                                             "spglibtestdata/virtual_structure/POSCAR-16-224-20",
                                             "spglibtestdata/virtual_structure/POSCAR-16-227-84",
                                             "spglibtestdata/virtual_structure/POSCAR-160-221-16",
                                             "spglibtestdata/virtual_structure/POSCAR-160-224-16",
                                             "spglibtestdata/virtual_structure/POSCAR-160-227-16",
                                             "spglibtestdata/virtual_structure/POSCAR-160-227-71",
                                             "spglibtestdata/virtual_structure/POSCAR-160-fcc",
                                             "spglibtestdata/virtual_structure/POSCAR-161-222-16",
                                             "spglibtestdata/virtual_structure/POSCAR-161-223-16",
                                             "spglibtestdata/virtual_structure/POSCAR-161-230-conv-18",
                                             "spglibtestdata/virtual_structure/POSCAR-161-230-prim-16",
                                             "spglibtestdata/virtual_structure/POSCAR-161-bcc-16",
                                             "spglibtestdata/virtual_structure/POSCAR-166-221-06",
                                             "spglibtestdata/virtual_structure/POSCAR-166-224-06",
                                             "spglibtestdata/virtual_structure/POSCAR-166-227-06",
                                             "spglibtestdata/virtual_structure/POSCAR-166-227-38",
                                             "spglibtestdata/virtual_structure/POSCAR-167-222-06",
                                             "spglibtestdata/virtual_structure/POSCAR-167-223-06",
                                             "spglibtestdata/virtual_structure/POSCAR-167-230-conv-06",
                                             "spglibtestdata/virtual_structure/POSCAR-167-230-prim-06",
                                             "spglibtestdata/virtual_structure/POSCAR-167-bcc-6",
                                             "spglibtestdata/virtual_structure/POSCAR-17-227-60",
                                             "spglibtestdata/virtual_structure/POSCAR-17-227-85",
                                             "spglibtestdata/virtual_structure/POSCAR-17-230-conv-46",
                                             "spglibtestdata/virtual_structure/POSCAR-18-227-86",
                                             "spglibtestdata/virtual_structure/POSCAR-19-227-59",
                                             "spglibtestdata/virtual_structure/POSCAR-19-227-89",
                                             "spglibtestdata/virtual_structure/POSCAR-19-230-conv-51",
                                             "spglibtestdata/virtual_structure/POSCAR-195-221-07",
                                             "spglibtestdata/virtual_structure/POSCAR-195-222-07",
                                             "spglibtestdata/virtual_structure/POSCAR-195-223-07",
                                             "spglibtestdata/virtual_structure/POSCAR-195-224-07",
                                             "spglibtestdata/virtual_structure/POSCAR-198-227-40",
                                             "spglibtestdata/virtual_structure/POSCAR-198-230-conv-20",
                                             "spglibtestdata/virtual_structure/POSCAR-199-230-conv-07",
                                             "spglibtestdata/virtual_structure/POSCAR-199-230-prim-07",
                                             "spglibtestdata/virtual_structure/POSCAR-2-221-28",
                                             "spglibtestdata/virtual_structure/POSCAR-2-222-28",
                                             "spglibtestdata/virtual_structure/POSCAR-2-223-28",
                                             "spglibtestdata/virtual_structure/POSCAR-2-224-28",
                                             "spglibtestdata/virtual_structure/POSCAR-2-227-41",
                                             "spglibtestdata/virtual_structure/POSCAR-2-227-74",
                                             "spglibtestdata/virtual_structure/POSCAR-2-227-94",
                                             "spglibtestdata/virtual_structure/POSCAR-2-230-conv-39",
                                             "spglibtestdata/virtual_structure/POSCAR-2-230-conv-57",
                                             "spglibtestdata/virtual_structure/POSCAR-2-230-prim-28",
                                             "spglibtestdata/virtual_structure/POSCAR-2-bcc-28",
                                             "spglibtestdata/virtual_structure/POSCAR-20-227-53",
                                             "spglibtestdata/virtual_structure/POSCAR-20-227-90",
                                             "spglibtestdata/virtual_structure/POSCAR-20-230-conv-53",
                                             "spglibtestdata/virtual_structure/POSCAR-200-221-02",
                                             "spglibtestdata/virtual_structure/POSCAR-200-223-02",
                                             "spglibtestdata/virtual_structure/POSCAR-201-222-02",
                                             "spglibtestdata/virtual_structure/POSCAR-201-224-02",
                                             "spglibtestdata/virtual_structure/POSCAR-205-230-conv-08",
                                             "spglibtestdata/virtual_structure/POSCAR-206-230-conv-02",
                                             "spglibtestdata/virtual_structure/POSCAR-206-230-prim-02",
                                             "spglibtestdata/virtual_structure/POSCAR-207-221-04",
                                             "spglibtestdata/virtual_structure/POSCAR-207-222-04",
                                             "spglibtestdata/virtual_structure/POSCAR-208-223-04",
                                             "spglibtestdata/virtual_structure/POSCAR-208-224-04",
                                             "spglibtestdata/virtual_structure/POSCAR-21-221-23",
                                             "spglibtestdata/virtual_structure/POSCAR-21-222-23",
                                             "spglibtestdata/virtual_structure/POSCAR-21-223-23",
                                             "spglibtestdata/virtual_structure/POSCAR-21-224-23",
                                             "spglibtestdata/virtual_structure/POSCAR-21-230-conv-49",
                                             "spglibtestdata/virtual_structure/POSCAR-212-227-19",
                                             "spglibtestdata/virtual_structure/POSCAR-213-230-conv-09",
                                             "spglibtestdata/virtual_structure/POSCAR-214-230-conv-04",
                                             "spglibtestdata/virtual_structure/POSCAR-214-230-prim-04",
                                             "spglibtestdata/virtual_structure/POSCAR-215-221-03",
                                             "spglibtestdata/virtual_structure/POSCAR-215-224-03",
                                             "spglibtestdata/virtual_structure/POSCAR-215-227-18",
                                             "spglibtestdata/virtual_structure/POSCAR-216-227-03",
                                             "spglibtestdata/virtual_structure/POSCAR-218-222-03",
                                             "spglibtestdata/virtual_structure/POSCAR-218-223-03",
                                             "spglibtestdata/virtual_structure/POSCAR-22-230-conv-26",
                                             "spglibtestdata/virtual_structure/POSCAR-22-230-prim-23",
                                             "spglibtestdata/virtual_structure/POSCAR-220-230-conv-03",
                                             "spglibtestdata/virtual_structure/POSCAR-220-230-prim-03",
                                             "spglibtestdata/virtual_structure/POSCAR-221-221-01",
                                             "spglibtestdata/virtual_structure/POSCAR-222-222-01",
                                             "spglibtestdata/virtual_structure/POSCAR-223-223-01",
                                             "spglibtestdata/virtual_structure/POSCAR-224-224-01",
                                             "spglibtestdata/virtual_structure/POSCAR-227-227-01",
                                             "spglibtestdata/virtual_structure/POSCAR-230-230-conv-01",
                                             "spglibtestdata/virtual_structure/POSCAR-230-230-conv-62",
                                             "spglibtestdata/virtual_structure/POSCAR-230-230-prim-01",
                                             "spglibtestdata/virtual_structure/POSCAR-24-230-conv-23",
                                             "spglibtestdata/virtual_structure/POSCAR-24-230-prim-20",
                                             "spglibtestdata/virtual_structure/POSCAR-25-221-21",
                                             "spglibtestdata/virtual_structure/POSCAR-25-223-21",
                                             "spglibtestdata/virtual_structure/POSCAR-25-227-54",
                                             "spglibtestdata/virtual_structure/POSCAR-26-227-64",
                                             "spglibtestdata/virtual_structure/POSCAR-27-230-conv-48",
                                             "spglibtestdata/virtual_structure/POSCAR-28-227-62",
                                             "spglibtestdata/virtual_structure/POSCAR-29-230-conv-52",
                                             "spglibtestdata/virtual_structure/POSCAR-3-221-29",
                                             "spglibtestdata/virtual_structure/POSCAR-3-222-29",
                                             "spglibtestdata/virtual_structure/POSCAR-3-223-29",
                                             "spglibtestdata/virtual_structure/POSCAR-3-224-29",
                                             "spglibtestdata/virtual_structure/POSCAR-3-227-82",
                                             "spglibtestdata/virtual_structure/POSCAR-3-227-95",
                                             "spglibtestdata/virtual_structure/POSCAR-3-230-conv-58",
                                             "spglibtestdata/virtual_structure/POSCAR-30-227-65",
                                             "spglibtestdata/virtual_structure/POSCAR-31-227-58",
                                             "spglibtestdata/virtual_structure/POSCAR-32-230-conv-47",
                                             "spglibtestdata/virtual_structure/POSCAR-33-227-63",
                                             "spglibtestdata/virtual_structure/POSCAR-34-222-21",
                                             "spglibtestdata/virtual_structure/POSCAR-34-224-21",
                                             "spglibtestdata/virtual_structure/POSCAR-35-221-22",
                                             "spglibtestdata/virtual_structure/POSCAR-35-224-22",
                                             "spglibtestdata/virtual_structure/POSCAR-35-227-87",
                                             "spglibtestdata/virtual_structure/POSCAR-37-222-22",
                                             "spglibtestdata/virtual_structure/POSCAR-37-223-22",
                                             "spglibtestdata/virtual_structure/POSCAR-38-221-26",
                                             "spglibtestdata/virtual_structure/POSCAR-39-224-26",
                                             "spglibtestdata/virtual_structure/POSCAR-4-227-77",
                                             "spglibtestdata/virtual_structure/POSCAR-4-227-81",
                                             "spglibtestdata/virtual_structure/POSCAR-4-227-96",
                                             "spglibtestdata/virtual_structure/POSCAR-4-230-conv-59",
                                             "spglibtestdata/virtual_structure/POSCAR-40-223-26",
                                             "spglibtestdata/virtual_structure/POSCAR-41-222-26",
                                             "spglibtestdata/virtual_structure/POSCAR-43-230-conv-25",
                                             "spglibtestdata/virtual_structure/POSCAR-43-230-conv-29",
                                             "spglibtestdata/virtual_structure/POSCAR-43-230-prim-22",
                                             "spglibtestdata/virtual_structure/POSCAR-43-230-prim-26",
                                             "spglibtestdata/virtual_structure/POSCAR-43-bcc-22",
                                             "spglibtestdata/virtual_structure/POSCAR-43-bcc-26",
                                             "spglibtestdata/virtual_structure/POSCAR-44-227-24",
                                             "spglibtestdata/virtual_structure/POSCAR-45-230-conv-24",
                                             "spglibtestdata/virtual_structure/POSCAR-45-230-prim-21",
                                             "spglibtestdata/virtual_structure/POSCAR-46-227-28",
                                             "spglibtestdata/virtual_structure/POSCAR-47-221-08",
                                             "spglibtestdata/virtual_structure/POSCAR-47-223-08",
                                             "spglibtestdata/virtual_structure/POSCAR-48-222-08",
                                             "spglibtestdata/virtual_structure/POSCAR-48-224-08",
                                             "spglibtestdata/virtual_structure/POSCAR-5-221-32",
                                             "spglibtestdata/virtual_structure/POSCAR-5-222-32",
                                             "spglibtestdata/virtual_structure/POSCAR-5-223-32",
                                             "spglibtestdata/virtual_structure/POSCAR-5-224-32",
                                             "spglibtestdata/virtual_structure/POSCAR-5-227-45",
                                             "spglibtestdata/virtual_structure/POSCAR-5-227-75",
                                             "spglibtestdata/virtual_structure/POSCAR-5-227-98",
                                             "spglibtestdata/virtual_structure/POSCAR-5-230-conv-40",
                                             "spglibtestdata/virtual_structure/POSCAR-5-230-conv-43",
                                             "spglibtestdata/virtual_structure/POSCAR-5-230-conv-61",
                                             "spglibtestdata/virtual_structure/POSCAR-5-230-prim-29",
                                             "spglibtestdata/virtual_structure/POSCAR-5-230-prim-32",
                                             "spglibtestdata/virtual_structure/POSCAR-5-bcc-29",
                                             "spglibtestdata/virtual_structure/POSCAR-5-bcc-32",
                                             "spglibtestdata/virtual_structure/POSCAR-51-227-29",
                                             "spglibtestdata/virtual_structure/POSCAR-53-227-32",
                                             "spglibtestdata/virtual_structure/POSCAR-54-230-conv-30",
                                             "spglibtestdata/virtual_structure/POSCAR-6-221-30",
                                             "spglibtestdata/virtual_structure/POSCAR-6-223-30",
                                             "spglibtestdata/virtual_structure/POSCAR-6-227-79",
                                             "spglibtestdata/virtual_structure/POSCAR-61-230-conv-31",
                                             "spglibtestdata/virtual_structure/POSCAR-62-227-31",
                                             "spglibtestdata/virtual_structure/POSCAR-65-221-09",
                                             "spglibtestdata/virtual_structure/POSCAR-66-223-09",
                                             "spglibtestdata/virtual_structure/POSCAR-67-224-09",
                                             "spglibtestdata/virtual_structure/POSCAR-68-222-09",
                                             "spglibtestdata/virtual_structure/POSCAR-7-222-30",
                                             "spglibtestdata/virtual_structure/POSCAR-7-224-30",
                                             "spglibtestdata/virtual_structure/POSCAR-7-227-78",
                                             "spglibtestdata/virtual_structure/POSCAR-7-227-80",
                                             "spglibtestdata/virtual_structure/POSCAR-7-230-conv-60",
                                             "spglibtestdata/virtual_structure/POSCAR-70-230-conv-11",
                                             "spglibtestdata/virtual_structure/POSCAR-70-230-prim-09",
                                             "spglibtestdata/virtual_structure/POSCAR-70-bcc-9",
                                             "spglibtestdata/virtual_structure/POSCAR-73-230-conv-10",
                                             "spglibtestdata/virtual_structure/POSCAR-73-230-prim-08",
                                             "spglibtestdata/virtual_structure/POSCAR-74-227-09",
                                             "spglibtestdata/virtual_structure/POSCAR-75-221-25",
                                             "spglibtestdata/virtual_structure/POSCAR-75-222-25",
                                             "spglibtestdata/virtual_structure/POSCAR-76-227-61",
                                             "spglibtestdata/virtual_structure/POSCAR-77-223-25",
                                             "spglibtestdata/virtual_structure/POSCAR-77-224-25",
                                             "spglibtestdata/virtual_structure/POSCAR-78-227-91",
                                             "spglibtestdata/virtual_structure/POSCAR-78-230-conv-54",
                                             "spglibtestdata/virtual_structure/POSCAR-8-221-31",
                                             "spglibtestdata/virtual_structure/POSCAR-8-224-31",
                                             "spglibtestdata/virtual_structure/POSCAR-8-227-44",
                                             "spglibtestdata/virtual_structure/POSCAR-8-227-97",
                                             "spglibtestdata/virtual_structure/POSCAR-80-230-conv-28",
                                             "spglibtestdata/virtual_structure/POSCAR-80-230-prim-25",
                                             "spglibtestdata/virtual_structure/POSCAR-81-221-24",
                                             "spglibtestdata/virtual_structure/POSCAR-81-222-24",
                                             "spglibtestdata/virtual_structure/POSCAR-81-223-24",
                                             "spglibtestdata/virtual_structure/POSCAR-81-224-24",
                                             "spglibtestdata/virtual_structure/POSCAR-81-227-88",
                                             "spglibtestdata/virtual_structure/POSCAR-81-230-conv-50",
                                             "spglibtestdata/virtual_structure/POSCAR-82-230-conv-27",
                                             "spglibtestdata/virtual_structure/POSCAR-82-230-prim-24",
                                             "spglibtestdata/virtual_structure/POSCAR-83-221-10",
                                             "spglibtestdata/virtual_structure/POSCAR-84-223-10",
                                             "spglibtestdata/virtual_structure/POSCAR-85-222-10",
                                             "spglibtestdata/virtual_structure/POSCAR-86-224-10",
                                             "spglibtestdata/virtual_structure/POSCAR-88-230-conv-12",
                                             "spglibtestdata/virtual_structure/POSCAR-88-230-prim-10",
                                             "spglibtestdata/virtual_structure/POSCAR-89-221-12",
                                             "spglibtestdata/virtual_structure/POSCAR-89-222-12",
                                             "spglibtestdata/virtual_structure/POSCAR-9-222-31",
                                             "spglibtestdata/virtual_structure/POSCAR-9-223-31",
                                             "spglibtestdata/virtual_structure/POSCAR-9-227-43",
                                             "spglibtestdata/virtual_structure/POSCAR-9-230-conv-41",
                                             "spglibtestdata/virtual_structure/POSCAR-9-230-conv-42",
                                             "spglibtestdata/virtual_structure/POSCAR-9-230-prim-30",
                                             "spglibtestdata/virtual_structure/POSCAR-9-230-prim-31",
                                             "spglibtestdata/virtual_structure/POSCAR-9-bcc-30",
                                             "spglibtestdata/virtual_structure/POSCAR-9-bcc-31",
                                             "spglibtestdata/virtual_structure/POSCAR-91-227-67",
                                             "spglibtestdata/virtual_structure/POSCAR-92-227-35",
                                             "spglibtestdata/virtual_structure/POSCAR-92-230-conv-35",
                                             "spglibtestdata/virtual_structure/POSCAR-93-223-12",
                                             "spglibtestdata/virtual_structure/POSCAR-93-224-12",
                                             "spglibtestdata/virtual_structure/POSCAR-95-227-36",
                                             "spglibtestdata/virtual_structure/POSCAR-95-230-conv-32",
                                             "spglibtestdata/virtual_structure/POSCAR-96-227-69",
                                             "spglibtestdata/virtual_structure/POSCAR-98-230-conv-14",
                                             "spglibtestdata/virtual_structure/POSCAR-98-230-prim-12",
                                             "spglibtestdata/virtual_structure/POSCAR-99-221-13"};

  for (const std::string& fileName : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double>> atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                   [randomShift](const std::tuple<double3, int, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup =
        SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      SKSymmetryCell cell = spaceGroup->cell;
      double3x3 transformationMatrix = spaceGroup->transformationMatrix;
      double3x3 rotationMatrix = spaceGroup->rotationMatrix;

      double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

      EXPECT_EQ(originalUnitCell, unitCell);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}
