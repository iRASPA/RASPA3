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
import skpointgroup;

TEST(FindPointgroupNoPartialOccupancies, Triclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {std::make_pair("spglibtestdata/triclinic/POSCAR-001", 1),
                                                  std::make_pair("spglibtestdata/triclinic/POSCAR-002", 2)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Monoclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {std::make_pair("spglibtestdata/monoclinic/POSCAR-003", 3),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-004", 3),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-004-2", 3),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-005", 3),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-005-2", 3),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-006", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-006-2", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-007", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-007-2", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-008", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-008-2", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-009", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-009-2", 4),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-010", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-010-2", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-011", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-011-2", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-012", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-012-2", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-012-3", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-013", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-013-2", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-013-3", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-014", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-014-2", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-015", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-015-2", 5),
                                                  std::make_pair("spglibtestdata/monoclinic/POSCAR-015-3", 5)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Orthogonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {std::make_pair("spglibtestdata/orthorhombic/POSCAR-016", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-016-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-017-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-018", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-018-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-019", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-019-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-020", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-021", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-021-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-022", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-023", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-023-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-024", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-024-2", 6),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-025", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-025-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-026", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-026-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-027", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-027-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-028", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-028-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-029", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-029-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-030", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-030-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-031", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-031-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-032", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-032-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-033", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-3", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-034", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-034-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-035", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-035-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-036", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-036-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-037", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-037-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-038", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-038-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-039", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-039-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-040", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-040-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-041", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-041-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-042", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-043", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-043-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-044", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-044-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-045", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-045-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-046", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-046-2", 7),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-047", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-047-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-048", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-048-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-049", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-049-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-050", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-050-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-051", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-052", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-052-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-053", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-053-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-054", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-054-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-055", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-055-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-056", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-056-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-057", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-057-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-058", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-059", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-059-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-060", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-061", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-061-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-062", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-062-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-063", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-064", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-065", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-066", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-066-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-067", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-3", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-068", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-068-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-069", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-069-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-070", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-070-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-071", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-071-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-072", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-072-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-073", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-073-2", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-074", 8),
                                                  std::make_pair("spglibtestdata/orthorhombic/POSCAR-074-2", 8)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Tetragonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {std::make_pair("spglibtestdata/tetragonal/POSCAR-075", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-075-2", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-076", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-076-2", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-077", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-077-2", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-077-3", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-078", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-078-2", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-079", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-079-2", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-080", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-080-2", 9),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-081", 10),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-081-2", 10),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-082", 10),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-082-2", 10),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-083", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-083-2", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-083-3", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-084", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-084-2", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-085", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-085-2", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-086", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-086-2", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-087", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-087-2", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-088", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-088-2", 11),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-090", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-090-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-091", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-091-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-092", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-092-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-092-3", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-094", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-094-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-094-3", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-095", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-095-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-096", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-096-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-097", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-097-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-098", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-098-2", 12),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-099", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-099-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-100", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-100-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-102", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-102-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-103", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-103-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-104", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-104-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-105", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-105-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-106", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-107", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-107-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-107-3", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-108", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-108-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-109", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-109-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-110", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-110-2", 13),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-111", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-111-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-112", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-112-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-113", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-113-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-114", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-114-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-115", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-115-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-115-3", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-115-4", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-115-5", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-116", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-116-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-117", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-117-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-118", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-118-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-119", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-119-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-120", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-120-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-121", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-121-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-122", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-122-2", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-122-3", 14),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-123", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-123-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-123-3", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-124", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-124-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-125", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-125-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-126", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-126-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-127", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-127-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-128", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-128-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-129", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-129-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-129-3", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-130", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-130-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-131", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-131-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-132", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-132-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-133", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-133-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-134", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-134-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-135", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-135-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-136", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-136-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-136-3", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-136-4", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-136-5", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-137", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-137-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-137-3", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-138", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-138-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-139", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-139-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-140", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-140-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-141", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-141-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-142", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-142-2", 15),
                                                  std::make_pair("spglibtestdata/tetragonal/POSCAR-142-3", 15)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Trigonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {std::make_pair("spglibtestdata/trigonal/POSCAR-143", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-143-2", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-144", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-144-2", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-145", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-145-2", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-146", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-146-2", 16),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-147", 17),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-147-2", 17),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-148", 17),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-148-2", 17),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-149", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-149-2", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-150", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-150-2", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-151", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-151-2", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-152", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-152-2", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-153", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-154", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-154-2", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-154-3", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-155", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-155-2", 18),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-156", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-156-2", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-157", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-157-2", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-158", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-158-2", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-159", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-159-2", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-160", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-160-2", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-161", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-161-2", 19),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-162", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-162-2", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-163", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-163-2", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-164", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-164-2", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-165", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-165-2", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-166", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-166-2", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-167", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-167-2", 20),
                                                  std::make_pair("spglibtestdata/trigonal/POSCAR-167-3", 20)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Hexagonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {std::make_pair("spglibtestdata/hexagonal/POSCAR-168", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-169", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-169-2", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-170", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-170-2", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-171", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-171-2", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-172", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-173", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-173-2", 21),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-174", 22),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-174-2", 22),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-175", 23),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-175-2", 23),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-176", 23),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-176-2", 23),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-177", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-179", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-179-2", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-180", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-180-2", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-181", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-181-2", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-182", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-182-2", 24),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-183", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-183-2", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-184", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-184-2", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-185", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-185-2", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-186", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-186-2", 25),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-187", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-187-2", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-188", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-188-2", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-189", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-189-2", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-190", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-190-2", 26),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-191", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-191-2", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-192", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-192-2", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-193", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-193-2", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-194", 27),
                                                  std::make_pair("spglibtestdata/hexagonal/POSCAR-194-2", 27)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Cubic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {
      std::make_pair("spglibtestdata/cubic/POSCAR-195", 28),   std::make_pair("spglibtestdata/cubic/POSCAR-195-2", 28),
      std::make_pair("spglibtestdata/cubic/POSCAR-196", 28),   std::make_pair("spglibtestdata/cubic/POSCAR-196-2", 28),
      std::make_pair("spglibtestdata/cubic/POSCAR-197", 28),   std::make_pair("spglibtestdata/cubic/POSCAR-197-2", 28),
      std::make_pair("spglibtestdata/cubic/POSCAR-198", 28),   std::make_pair("spglibtestdata/cubic/POSCAR-198-2", 28),
      std::make_pair("spglibtestdata/cubic/POSCAR-199", 28),   std::make_pair("spglibtestdata/cubic/POSCAR-199-2", 28),
      std::make_pair("spglibtestdata/cubic/POSCAR-200", 29),   std::make_pair("spglibtestdata/cubic/POSCAR-200-2", 29),
      std::make_pair("spglibtestdata/cubic/POSCAR-205", 29),   std::make_pair("spglibtestdata/cubic/POSCAR-205-3", 29),
      std::make_pair("spglibtestdata/cubic/POSCAR-206", 29),   std::make_pair("spglibtestdata/cubic/POSCAR-206-2", 29),
      std::make_pair("spglibtestdata/cubic/POSCAR-207", 30),   std::make_pair("spglibtestdata/cubic/POSCAR-208", 30),
      std::make_pair("spglibtestdata/cubic/POSCAR-208-2", 30), std::make_pair("spglibtestdata/cubic/POSCAR-209", 30),
      std::make_pair("spglibtestdata/cubic/POSCAR-210", 30),   std::make_pair("spglibtestdata/cubic/POSCAR-210-2", 30),
      std::make_pair("spglibtestdata/cubic/POSCAR-211", 30),   std::make_pair("spglibtestdata/cubic/POSCAR-212", 30),
      std::make_pair("spglibtestdata/cubic/POSCAR-212-2", 30), std::make_pair("spglibtestdata/cubic/POSCAR-213", 30),
      std::make_pair("spglibtestdata/cubic/POSCAR-213-2", 30), std::make_pair("spglibtestdata/cubic/POSCAR-214", 30),
      std::make_pair("spglibtestdata/cubic/POSCAR-214-2", 30), std::make_pair("spglibtestdata/cubic/POSCAR-215", 31),
      std::make_pair("spglibtestdata/cubic/POSCAR-215-2", 31), std::make_pair("spglibtestdata/cubic/POSCAR-216", 31),
      std::make_pair("spglibtestdata/cubic/POSCAR-216-2", 31), std::make_pair("spglibtestdata/cubic/POSCAR-217", 31),
      std::make_pair("spglibtestdata/cubic/POSCAR-217-2", 31), std::make_pair("spglibtestdata/cubic/POSCAR-218", 31),
      std::make_pair("spglibtestdata/cubic/POSCAR-218-2", 31), std::make_pair("spglibtestdata/cubic/POSCAR-219", 31),
      std::make_pair("spglibtestdata/cubic/POSCAR-219-2", 31), std::make_pair("spglibtestdata/cubic/POSCAR-220", 31),
      std::make_pair("spglibtestdata/cubic/POSCAR-220-2", 31), std::make_pair("spglibtestdata/cubic/POSCAR-221", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-221-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-222", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-222-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-223", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-223-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-224", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-224-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-225", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-225-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-226", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-226-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-227", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-227-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-228", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-228-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-229", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-229-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-230", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-230-2", 32), std::make_pair("spglibtestdata/cubic/POSCAR-230-3", 32),
      std::make_pair("spglibtestdata/cubic/POSCAR-230-4", 32)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindPointgroupNoPartialOccupancies, Virtual)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, size_t> testData = {
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-221-33", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-222-33", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-223-33", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-224-33", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-73", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-93", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-99", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-conv-56", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-prim-33", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-bcc-33", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-221-18", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-223-18", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-227-50", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-102-224-13", 13),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-104-222-13", 13),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-105-223-13", 13),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-109-227-13", 13),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-11-227-48", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-conv-15", 13),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-prim-13", 13),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-221-11", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-224-11", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-227-66", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-222-11", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-223-11", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-113-227-68", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-221-14", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-223-14", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-227-33", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-116-230-conv-34", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-117-230-conv-33", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-222-14", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-224-14", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-221-19", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-224-19", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-21", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-83", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-conv-16", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-prim-14", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-conv-13", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-prim-11", 14),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-123-221-05", 15),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-126-222-05", 15),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-222-18", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-224-18", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-227-49", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-230-conv-44", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-131-223-05", 15),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-134-224-05", 15),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-47", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-51", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-230-conv-45", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-conv-05", 15),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-prim-05", 15),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-221-27", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-222-27", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-223-27", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-224-27", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-227-92", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-36", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-55", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-prim-27", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-bcc-27", 16),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-221-15", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-222-15", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-223-15", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-224-15", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-227-70", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-17", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-37", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-prim-15", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-bcc-15", 17),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-222-19", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-223-19", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-21", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-22", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-18", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-19", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-18", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-19", 5),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-221-17", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-222-17", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-223-17", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-224-17", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-227-72", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-19", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-38", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-prim-17", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-bcc-17", 18),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-221-20", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-222-20", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-223-20", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-224-20", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-227-84", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-221-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-224-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-71", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-fcc", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-222-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-223-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-conv-18", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-prim-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-bcc-16", 19),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-221-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-224-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-38", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-222-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-223-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-conv-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-prim-06", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-bcc-6", 20),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-60", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-85", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-230-conv-46", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-18-227-86", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-59", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-89", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-230-conv-51", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-221-07", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-222-07", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-223-07", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-224-07", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-227-40", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-230-conv-20", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-conv-07", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-prim-07", 28),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-221-28", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-222-28", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-223-28", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-224-28", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-41", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-74", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-94", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-39", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-57", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-prim-28", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-bcc-28", 2),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-53", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-90", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-230-conv-53", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-221-02", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-223-02", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-222-02", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-224-02", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-205-230-conv-08", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-conv-02", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-prim-02", 29),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-221-04", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-222-04", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-223-04", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-224-04", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-221-23", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-222-23", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-223-23", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-224-23", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-230-conv-49", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-212-227-19", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-213-230-conv-09", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-conv-04", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-prim-04", 30),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-221-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-224-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-227-18", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-216-227-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-222-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-223-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-conv-26", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-prim-23", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-conv-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-prim-03", 31),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-221-221-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-222-222-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-223-223-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-224-224-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-227-227-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-62", 1),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-prim-01", 32),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-24-230-conv-23", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-24-230-prim-20", 6),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-221-21", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-223-21", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-227-54", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-26-227-64", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-27-230-conv-48", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-28-227-62", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-29-230-conv-52", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-221-29", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-222-29", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-223-29", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-224-29", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-227-82", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-227-95", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-230-conv-58", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-30-227-65", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-31-227-58", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-32-230-conv-47", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-33-227-63", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-34-222-21", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-34-224-21", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-221-22", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-224-22", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-227-87", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-37-222-22", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-37-223-22", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-38-221-26", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-39-224-26", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-77", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-81", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-96", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-230-conv-59", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-40-223-26", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-41-222-26", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-conv-25", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-conv-29", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-prim-22", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-prim-26", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-bcc-22", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-bcc-26", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-44-227-24", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-45-230-conv-24", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-45-230-prim-21", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-46-227-28", 7),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-47-221-08", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-47-223-08", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-48-222-08", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-48-224-08", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-221-32", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-222-32", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-223-32", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-224-32", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-45", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-75", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-98", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-40", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-43", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-61", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-prim-29", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-prim-32", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-bcc-29", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-bcc-32", 3),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-51-227-29", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-53-227-32", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-54-230-conv-30", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-221-30", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-223-30", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-227-79", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-61-230-conv-31", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-62-227-31", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-65-221-09", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-66-223-09", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-67-224-09", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-68-222-09", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-222-30", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-224-30", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-227-78", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-227-80", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-230-conv-60", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-230-conv-11", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-230-prim-09", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-bcc-9", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-73-230-conv-10", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-73-230-prim-08", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-74-227-09", 8),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-75-221-25", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-75-222-25", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-76-227-61", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-77-223-25", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-77-224-25", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-78-227-91", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-78-230-conv-54", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-221-31", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-224-31", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-227-44", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-227-97", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-80-230-conv-28", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-80-230-prim-25", 9),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-221-24", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-222-24", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-223-24", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-224-24", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-227-88", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-230-conv-50", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-82-230-conv-27", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-82-230-prim-24", 10),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-83-221-10", 11),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-84-223-10", 11),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-85-222-10", 11),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-86-224-10", 11),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-88-230-conv-12", 11),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-88-230-prim-10", 11),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-89-221-12", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-89-222-12", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-222-31", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-223-31", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-227-43", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-conv-41", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-conv-42", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-prim-30", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-prim-31", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-bcc-30", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-bcc-31", 4),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-91-227-67", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-92-227-35", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-92-230-conv-35", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-93-223-12", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-93-224-12", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-95-227-36", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-95-230-conv-32", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-96-227-69", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-98-230-conv-14", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-98-230-prim-12", 12),
      std::make_pair("spglibtestdata/virtual_structure/POSCAR-99-221-13", 13)};

  for (auto& [fileName, pointGroupTargetValue] : testData)
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
                   [randomShift](const std::tuple<double3, size_t, double>& atom)
                   { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKPointGroup> pointGroup =
        SKSpaceGroup::findPointGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (pointGroup)
    {
      size_t foundPointGroup = pointGroup->number();
      EXPECT_EQ(foundPointGroup, pointGroupTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}
