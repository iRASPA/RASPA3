#include <gtest/gtest.h>
#include <filesystem>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <tuple>
#include <random>
#include <optional>
#include <format>
#include <print>

import double3;
import double3x3;
import randomnumbers;

import skposcarlegacyparser;
import sksymmetrycell;
import skspacegroup;
import skpointgroup;

TEST(FindSpacegroup, Triclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/triclinic/POSCAR-001" , 1),
    std::make_pair("spglibtestdata/triclinic/POSCAR-002" , 2)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Monoclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/monoclinic/POSCAR-003", 3),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-004", 4),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-004-2",  4),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-005",  5),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-005-2",  5),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-006",  6),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-006-2",  6),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-007",  7),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-007-2",  7),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-008",  8),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-008-2",  8),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-009",  9),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-009-2",  9),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-010",  10),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-010-2",  10),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-011",  11),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-011-2",  11),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012",  12),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012-2",  12),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012-3",  12),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013",  13),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013-2",  13),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013-3",  13),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-014",  14),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-014-2",  14),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015",  15),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015-2",  15),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015-3",  15)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Orthorhombic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-016",  16),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-016-2",  16),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-017-2",  17),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-018",  18),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-018-2",  18),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-019",  19),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-019-2",  19),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-020",  20),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-021",  21),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-021-2",  21),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-022",  22),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-023",  23),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-023-2",  23),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-024",  24),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-024-2",  24),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-025",  51),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-025-2",  25),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-026",  26),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-026-2",  26),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-027",  27),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-027-2",  27),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-028",  28),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-028-2",  28),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-029",  29),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-029-2",  29),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-030",  30),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-030-2",  30),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-031",  31),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-031-2",  31),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-032",  32),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-032-2",  32),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033",  33),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-2",  33),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-3",  33),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-034",  34),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-034-2",  34),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-035",  35),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-035-2",  35),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-036",  36),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-036-2",  36),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-037",  37),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-037-2",  37),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-038",  38),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-038-2",  38),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-039",  39),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-039-2",  39),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-040",  40),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-040-2",  40),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-041",  41),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-041-2",  41),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-042",  42),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-043",  43),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-043-2",  43),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-044",  44),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-044-2",  44),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-045",  45),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-045-2",  45),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-046",  46),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-046-2",  46),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-047",  47),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-047-2",  65),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-048",  48),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-048-2",  48),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-049",  49),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-049-2",  49),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-050",  50),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-050-2",  50),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051",  51),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-2",  51),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-3",  51),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-052",  52),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-052-2",  52),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-053",  53),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-053-2",  53),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-054",  54),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-054-2",  54),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-055",  55),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-055-2",  55),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-056",  56),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-056-2",  56),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-057",  57),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-057-2",  57),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058",  58),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-2",  58),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-3",  58),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-059",  59),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-059-2",  59),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060",  60),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-2",  60),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-3",  60),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-061",  61),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-061-2",  61),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-062",  62),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-062-2",  62),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063",  63),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-2",  63),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-3",  63),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064",  64),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-2",  64),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-3",  64),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065",  65),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-2",  65),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-3",  65),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-066",  66),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-066-2",  66),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067",  67),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-2",  67),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-3",  67),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-068",  68),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-068-2",  68),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-069",  69),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-069-2",  69),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-070",  70),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-070-2",  70),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-071",  71),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-071-2",  71),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-072",  72),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-072-2",  72),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-073",  73),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-073-2",  73),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-074",  74),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-074-2",  74)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Tetragonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/tetragonal/POSCAR-075",  75),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-075-2",  75),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-076",  76),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-076-2",  76),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077",  77),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077-2",  77),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077-3",  77),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-078",  78),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-078-2",  78),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-079",  79),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-079-2",  79),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-080",  80),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-080-2",  80),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-081",  81),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-081-2",  81),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-082",  121),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-082-2",  82),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083",  83),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083-2",  139),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083-3",  83),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-084",  84),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-084-2",  84),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-085",  85),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-085-2",  85),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-086",  86),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-086-2",  86),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-087",  87),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-087-2",  87),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-088",  88),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-088-2",  88),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-090",  90),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-090-2",  90),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-091",  91),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-091-2",  91),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092",  92),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092-2",  92),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092-3",  92),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094",  94),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094-2",  94),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094-3",  94),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-095",  95),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-095-2",  95),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-096",  96),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-096-2",  96),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-097",  97),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-097-2",  97),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-098",  98),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-098-2",  98),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-099",  99),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-099-2",  99),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-100",  100),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-100-2",  100),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-102",  102),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-102-2",  102),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-103",  103),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-103-2",  103),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-104",  104),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-104-2",  104),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-105",  105),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-105-2",  105),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-106",  106),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107",  107),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107-2",  107),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107-3",  107),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-108",  108),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-108-2",  108),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-109",  141),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-109-2",  109),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-110",  110),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-110-2",  110),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-111",  111),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-111-2",  111),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-112",  112),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-112-2",  112),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-113",  113),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-113-2",  113),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-114",  114),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-114-2",  114),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115",  115),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-2",  115),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-3",  115),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-4",  115),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-5",  115),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-116",  116),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-116-2",  116),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-117",  117),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-117-2",  117),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-118",  118),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-118-2",  118),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-119",  119),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-119-2",  119),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-120",  120),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-120-2",  120),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-121",  121),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-121-2",  121),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122",  122),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122-2",  122),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122-3",  122),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123",  139),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123-2",  123),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123-3",  123),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-124",  140),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-124-2",  124),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-125",  125),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-125-2",  125),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-126",  126),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-126-2",  126),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-127",  127),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-127-2",  127),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-128",  128),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-128-2",  128),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129",  129),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129-2",  129),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129-3",  129),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-130",  130),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-130-2",  130),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-131",  131),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-131-2",  131),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-132",  132),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-132-2",  132),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-133",  133),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-133-2",  133),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-134",  134),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-134-2",  139),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-135",  135),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-135-2",  135),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136",  136),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-2",  136),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-3",  136),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-4",  136),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-5",  136),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137",  137),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137-2",  137),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137-3",  137),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-138",  138),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-138-2",  138),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-139",  139),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-139-2",  139),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-140",  140),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-140-2",  140),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-141",  141),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-141-2",  141),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142",  142),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142-2",  142),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142-3",  142)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Trigonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/trigonal/POSCAR-143",  143),
    std::make_pair("spglibtestdata/trigonal/POSCAR-143-2",  143),
    std::make_pair("spglibtestdata/trigonal/POSCAR-144",  144),
    std::make_pair("spglibtestdata/trigonal/POSCAR-144-2",  144),
    std::make_pair("spglibtestdata/trigonal/POSCAR-145",  145),
    std::make_pair("spglibtestdata/trigonal/POSCAR-145-2",  145),
    std::make_pair("spglibtestdata/trigonal/POSCAR-146",  146),
    std::make_pair("spglibtestdata/trigonal/POSCAR-146-2",  146),
    std::make_pair("spglibtestdata/trigonal/POSCAR-147",  147),
    std::make_pair("spglibtestdata/trigonal/POSCAR-147-2",  147),
    std::make_pair("spglibtestdata/trigonal/POSCAR-148",  148),
    std::make_pair("spglibtestdata/trigonal/POSCAR-148-2",  148),
    std::make_pair("spglibtestdata/trigonal/POSCAR-149",  149),
    std::make_pair("spglibtestdata/trigonal/POSCAR-149-2",  149),
    std::make_pair("spglibtestdata/trigonal/POSCAR-150",  150),
    std::make_pair("spglibtestdata/trigonal/POSCAR-150-2",  150),
    std::make_pair("spglibtestdata/trigonal/POSCAR-151",  151),
    std::make_pair("spglibtestdata/trigonal/POSCAR-151-2",  151),
    std::make_pair("spglibtestdata/trigonal/POSCAR-152",  152),
    std::make_pair("spglibtestdata/trigonal/POSCAR-152-2",  152),
    std::make_pair("spglibtestdata/trigonal/POSCAR-153",  153),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154",  154),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154-2",  154),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154-3",  154),
    std::make_pair("spglibtestdata/trigonal/POSCAR-155",  155),
    std::make_pair("spglibtestdata/trigonal/POSCAR-155-2",  155),
    std::make_pair("spglibtestdata/trigonal/POSCAR-156",  156),
    std::make_pair("spglibtestdata/trigonal/POSCAR-156-2",  156),
    std::make_pair("spglibtestdata/trigonal/POSCAR-157",  157),
    std::make_pair("spglibtestdata/trigonal/POSCAR-157-2",  157),
    std::make_pair("spglibtestdata/trigonal/POSCAR-158",  158),
    std::make_pair("spglibtestdata/trigonal/POSCAR-158-2",  158),
    std::make_pair("spglibtestdata/trigonal/POSCAR-159",  159),
    std::make_pair("spglibtestdata/trigonal/POSCAR-159-2",  159),
    std::make_pair("spglibtestdata/trigonal/POSCAR-160",  160),
    std::make_pair("spglibtestdata/trigonal/POSCAR-160-2",  160),
    std::make_pair("spglibtestdata/trigonal/POSCAR-161",  161),
    std::make_pair("spglibtestdata/trigonal/POSCAR-161-2",  161),
    std::make_pair("spglibtestdata/trigonal/POSCAR-162",  162),
    std::make_pair("spglibtestdata/trigonal/POSCAR-162-2",  162),
    std::make_pair("spglibtestdata/trigonal/POSCAR-163",  163),
    std::make_pair("spglibtestdata/trigonal/POSCAR-163-2",  163),
    std::make_pair("spglibtestdata/trigonal/POSCAR-164",  164),
    std::make_pair("spglibtestdata/trigonal/POSCAR-164-2",  164),
    std::make_pair("spglibtestdata/trigonal/POSCAR-165",  165),
    std::make_pair("spglibtestdata/trigonal/POSCAR-165-2",  165),
    std::make_pair("spglibtestdata/trigonal/POSCAR-166",  166),
    std::make_pair("spglibtestdata/trigonal/POSCAR-166-2",  166),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167",  167),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167-2",  167),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167-3",  167)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Hexagonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/hexagonal/POSCAR-168",  168),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-169",  169),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-169-2",  169),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-170",  170),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-170-2",  170),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-171",  171),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-171-2",  171),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-172",  172),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-173",  173),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-173-2",  173),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-174",  174),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-174-2",  174),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-175",  175),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-175-2",  175),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-176",  176),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-176-2",  176),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-177",  177),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-179",  179),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-179-2",  179),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-180",  180),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-180-2",  180),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-181",  181),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-181-2",  181),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-182",  182),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-182-2",  182),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-183",  183),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-183-2",  183),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-184",  184),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-184-2",  184),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-185",  185),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-185-2",  185),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-186",  186),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-186-2",  186),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-187",  194),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-187-2",  187),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-188",  188),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-188-2",  188),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-189",  189),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-189-2",  189),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-190",  190),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-190-2",  190),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-191",  191),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-191-2",  191),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-192",  192),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-192-2",  192),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-193",  193),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-193-2",  193),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-194",  194),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-194-2",  194)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Cubic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/cubic/POSCAR-195",  195),
    std::make_pair("spglibtestdata/cubic/POSCAR-195-2",  195),
    std::make_pair("spglibtestdata/cubic/POSCAR-196",  196),
    std::make_pair("spglibtestdata/cubic/POSCAR-196-2",  196),
    std::make_pair("spglibtestdata/cubic/POSCAR-197",  197),
    std::make_pair("spglibtestdata/cubic/POSCAR-197-2",  197),
    std::make_pair("spglibtestdata/cubic/POSCAR-198",  198),
    std::make_pair("spglibtestdata/cubic/POSCAR-198-2",  198),
    std::make_pair("spglibtestdata/cubic/POSCAR-199",  199),
    std::make_pair("spglibtestdata/cubic/POSCAR-199-2",  199),
    std::make_pair("spglibtestdata/cubic/POSCAR-200",  200),
    std::make_pair("spglibtestdata/cubic/POSCAR-200-2",  200),
    std::make_pair("spglibtestdata/cubic/POSCAR-205",  205),
    std::make_pair("spglibtestdata/cubic/POSCAR-205-3",  205),
    std::make_pair("spglibtestdata/cubic/POSCAR-206",  206),
    std::make_pair("spglibtestdata/cubic/POSCAR-206-2",  206),
    std::make_pair("spglibtestdata/cubic/POSCAR-207",  207),
    std::make_pair("spglibtestdata/cubic/POSCAR-208",  208),
    std::make_pair("spglibtestdata/cubic/POSCAR-208-2",  221),
    std::make_pair("spglibtestdata/cubic/POSCAR-209",  209),
    std::make_pair("spglibtestdata/cubic/POSCAR-210",  210),
    std::make_pair("spglibtestdata/cubic/POSCAR-210-2",  210),
    std::make_pair("spglibtestdata/cubic/POSCAR-211",  211),
    std::make_pair("spglibtestdata/cubic/POSCAR-212",  212),
    std::make_pair("spglibtestdata/cubic/POSCAR-212-2",  212),
    std::make_pair("spglibtestdata/cubic/POSCAR-213",  213),
    std::make_pair("spglibtestdata/cubic/POSCAR-213-2",  213),
    std::make_pair("spglibtestdata/cubic/POSCAR-214",  214),
    std::make_pair("spglibtestdata/cubic/POSCAR-214-2",  214),
    std::make_pair("spglibtestdata/cubic/POSCAR-215",  215),
    std::make_pair("spglibtestdata/cubic/POSCAR-215-2",  215),
    std::make_pair("spglibtestdata/cubic/POSCAR-216",  227),
    std::make_pair("spglibtestdata/cubic/POSCAR-216-2",  216),
    std::make_pair("spglibtestdata/cubic/POSCAR-217",  217),
    std::make_pair("spglibtestdata/cubic/POSCAR-217-2",  217),
    std::make_pair("spglibtestdata/cubic/POSCAR-218",  218),
    std::make_pair("spglibtestdata/cubic/POSCAR-218-2",  218),
    std::make_pair("spglibtestdata/cubic/POSCAR-219",  219),
    std::make_pair("spglibtestdata/cubic/POSCAR-219-2",  219),
    std::make_pair("spglibtestdata/cubic/POSCAR-220",  220),
    std::make_pair("spglibtestdata/cubic/POSCAR-220-2",  220),
    std::make_pair("spglibtestdata/cubic/POSCAR-221",  221),
    std::make_pair("spglibtestdata/cubic/POSCAR-221-2",  221),
    std::make_pair("spglibtestdata/cubic/POSCAR-222",  222),
    std::make_pair("spglibtestdata/cubic/POSCAR-222-2",  222),
    std::make_pair("spglibtestdata/cubic/POSCAR-223",  223),
    std::make_pair("spglibtestdata/cubic/POSCAR-223-2",  223),
    std::make_pair("spglibtestdata/cubic/POSCAR-224",  224),
    std::make_pair("spglibtestdata/cubic/POSCAR-224-2",  224),
    std::make_pair("spglibtestdata/cubic/POSCAR-225",  225),
    std::make_pair("spglibtestdata/cubic/POSCAR-225-2",  221),
    std::make_pair("spglibtestdata/cubic/POSCAR-226",  226),
    std::make_pair("spglibtestdata/cubic/POSCAR-226-2",  226),
    std::make_pair("spglibtestdata/cubic/POSCAR-227",  227),
    std::make_pair("spglibtestdata/cubic/POSCAR-227-2",  227),
    std::make_pair("spglibtestdata/cubic/POSCAR-228",  228),
    std::make_pair("spglibtestdata/cubic/POSCAR-228-2",  228),
    std::make_pair("spglibtestdata/cubic/POSCAR-229",  229),
    std::make_pair("spglibtestdata/cubic/POSCAR-229-2",  229),
    std::make_pair("spglibtestdata/cubic/POSCAR-230",  230),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-2",  230),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-3",  230),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-4",  230)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSpacegroup, Virtual)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData =
  {
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-221-33",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-222-33",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-223-33",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-224-33",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-73",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-93",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-99",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-conv-56",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-prim-33",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-bcc-33",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-221-18",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-223-18",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-227-50",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-102-224-13",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-104-222-13",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-105-223-13",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-109-227-13",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-11-227-48",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-conv-15",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-prim-13",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-221-11",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-224-11",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-227-66",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-222-11",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-223-11",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-113-227-68",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-221-14",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-223-14",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-227-33",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-116-230-conv-34",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-117-230-conv-33",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-222-14",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-224-14",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-221-19",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-224-19",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-21",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-83",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-conv-16",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-prim-14",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-conv-13",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-prim-11",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-123-221-05",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-126-222-05",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-222-18",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-224-18",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-227-49",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-230-conv-44",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-131-223-05",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-134-224-05",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-47",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-51",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-230-conv-45",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-conv-05",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-prim-05",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-221-27",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-222-27",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-223-27",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-224-27",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-227-92",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-36",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-55",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-prim-27",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-bcc-27",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-221-15",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-222-15",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-223-15",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-224-15",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-227-70",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-17",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-37",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-prim-15",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-bcc-15",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-222-19",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-223-19",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-21",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-22",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-18",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-19",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-18",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-19",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-221-17",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-222-17",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-223-17",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-224-17",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-227-72",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-19",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-38",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-prim-17",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-bcc-17",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-221-20",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-222-20",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-223-20",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-224-20",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-227-84",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-221-16",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-224-16",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-16",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-71",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-fcc",  166),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-222-16",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-223-16",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-conv-18",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-prim-16",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-bcc-16",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-221-06",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-224-06",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-06",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-38",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-222-06",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-223-06",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-conv-06",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-prim-06",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-bcc-6",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-60",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-85",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-230-conv-46",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-18-227-86",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-59",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-89",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-230-conv-51",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-221-07",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-222-07",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-223-07",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-224-07",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-227-40",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-230-conv-20",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-conv-07",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-prim-07",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-221-28",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-222-28",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-223-28",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-224-28",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-41",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-74",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-94",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-39",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-57",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-prim-28",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-bcc-28",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-53",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-90",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-230-conv-53",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-221-02",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-223-02",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-222-02",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-224-02",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-205-230-conv-08",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-conv-02",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-prim-02",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-221-04",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-222-04",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-223-04",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-224-04",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-221-23",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-222-23",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-223-23",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-224-23",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-230-conv-49",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-212-227-19",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-213-230-conv-09",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-conv-04",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-prim-04",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-221-03",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-224-03",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-227-18",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-216-227-03",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-222-03",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-223-03",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-conv-26",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-prim-23",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-conv-03",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-prim-03",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-221-221-01",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-222-222-01",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-223-223-01",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-224-224-01",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-227-227-01",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-01",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-62",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-prim-01",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-24-230-conv-23",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-24-230-prim-20",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-221-21",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-223-21",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-227-54",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-26-227-64",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-27-230-conv-48",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-28-227-62",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-29-230-conv-52",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-221-29",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-222-29",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-223-29",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-224-29",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-227-82",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-227-95",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-230-conv-58",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-30-227-65",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-31-227-58",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-32-230-conv-47",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-33-227-63",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-34-222-21",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-34-224-21",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-221-22",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-224-22",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-227-87",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-37-222-22",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-37-223-22",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-38-221-26",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-39-224-26",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-77",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-81",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-96",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-230-conv-59",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-40-223-26",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-41-222-26",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-conv-25",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-conv-29",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-prim-22",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-prim-26",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-bcc-22",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-bcc-26",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-44-227-24",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-45-230-conv-24",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-45-230-prim-21",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-46-227-28",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-47-221-08",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-47-223-08",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-48-222-08",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-48-224-08",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-221-32",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-222-32",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-223-32",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-224-32",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-45",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-75",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-98",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-40",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-43",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-61",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-prim-29",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-prim-32",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-bcc-29",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-bcc-32",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-51-227-29",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-53-227-32",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-54-230-conv-30",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-221-30",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-223-30",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-227-79",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-61-230-conv-31",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-62-227-31",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-65-221-09",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-66-223-09",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-67-224-09",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-68-222-09",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-222-30",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-224-30",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-227-78",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-227-80",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-230-conv-60",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-230-conv-11",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-230-prim-09",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-bcc-9",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-73-230-conv-10",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-73-230-prim-08",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-74-227-09",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-75-221-25",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-75-222-25",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-76-227-61",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-77-223-25",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-77-224-25",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-78-227-91",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-78-230-conv-54",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-221-31",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-224-31",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-227-44",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-227-97",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-80-230-conv-28",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-80-230-prim-25",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-221-24",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-222-24",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-223-24",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-224-24",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-227-88",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-230-conv-50",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-82-230-conv-27",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-82-230-prim-24",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-83-221-10",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-84-223-10",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-85-222-10",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-86-224-10",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-88-230-conv-12",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-88-230-prim-10",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-89-221-12",  221),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-89-222-12",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-222-31",  222),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-223-31",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-227-43",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-conv-41",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-conv-42",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-prim-30",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-prim-31",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-bcc-30",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-bcc-31",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-91-227-67",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-92-227-35",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-92-230-conv-35",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-93-223-12",  223),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-93-224-12",  224),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-95-227-36",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-95-230-conv-32",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-96-227-69",  227),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-98-230-conv-14",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-98-230-prim-12",  230),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-99-221-13",  221)
  };

  for (auto& [fileName, spaceGroupTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = true;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      size_t foundHallNumber = spaceGroup->HallNumber;
      SKSpaceGroup foundSpaceGroup = SKSpaceGroup(foundHallNumber);
      size_t foundSpaceGroupNumber = foundSpaceGroup.spaceGroupSetting().number();
      EXPECT_EQ(foundSpaceGroupNumber, spaceGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}
