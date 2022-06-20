#include "CppUnitTest.h"

import <map>;
import <tuple>;
import <vector>;
import <iostream>;
import <fstream>;
import <streambuf>;
import <algorithm>;
import <format>;
import <string>;
import <random>;
import symmetrykit;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;



namespace symmetrykittests
{
    TEST_CLASS(symmetrykit_test_transformationmatrix)
    {
    public:

        TEST_METHOD(TestTriclinic)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/triclinic/POSCAR-001",
                L"spglibtestdata/triclinic/POSCAR-002"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;
                    
                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestMonoclinic)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/monoclinic/POSCAR-003",
                L"spglibtestdata/monoclinic/POSCAR-004",
                L"spglibtestdata/monoclinic/POSCAR-004-2",
                L"spglibtestdata/monoclinic/POSCAR-005",
                L"spglibtestdata/monoclinic/POSCAR-005-2",
                L"spglibtestdata/monoclinic/POSCAR-006",
                L"spglibtestdata/monoclinic/POSCAR-006-2",
                L"spglibtestdata/monoclinic/POSCAR-007",
                L"spglibtestdata/monoclinic/POSCAR-007-2",
                L"spglibtestdata/monoclinic/POSCAR-008",
                L"spglibtestdata/monoclinic/POSCAR-008-2",
                L"spglibtestdata/monoclinic/POSCAR-009",
                L"spglibtestdata/monoclinic/POSCAR-009-2",
                L"spglibtestdata/monoclinic/POSCAR-010",
                L"spglibtestdata/monoclinic/POSCAR-010-2",
                L"spglibtestdata/monoclinic/POSCAR-011",
                L"spglibtestdata/monoclinic/POSCAR-011-2",
                L"spglibtestdata/monoclinic/POSCAR-012",
                L"spglibtestdata/monoclinic/POSCAR-012-2",
                L"spglibtestdata/monoclinic/POSCAR-012-3",
                L"spglibtestdata/monoclinic/POSCAR-013",
                L"spglibtestdata/monoclinic/POSCAR-013-2",
                L"spglibtestdata/monoclinic/POSCAR-013-3",
                L"spglibtestdata/monoclinic/POSCAR-014",
                L"spglibtestdata/monoclinic/POSCAR-014-2",
                L"spglibtestdata/monoclinic/POSCAR-015",
                L"spglibtestdata/monoclinic/POSCAR-015-2",
                L"spglibtestdata/monoclinic/POSCAR-015-3"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestOrthorhombic)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/orthorhombic/POSCAR-016",
                L"spglibtestdata/orthorhombic/POSCAR-016-2",
                L"spglibtestdata/orthorhombic/POSCAR-017-2",
                L"spglibtestdata/orthorhombic/POSCAR-018",
                L"spglibtestdata/orthorhombic/POSCAR-018-2",
                L"spglibtestdata/orthorhombic/POSCAR-019",
                L"spglibtestdata/orthorhombic/POSCAR-019-2",
                L"spglibtestdata/orthorhombic/POSCAR-020",
                L"spglibtestdata/orthorhombic/POSCAR-021",
                L"spglibtestdata/orthorhombic/POSCAR-021-2",
                L"spglibtestdata/orthorhombic/POSCAR-022",
                L"spglibtestdata/orthorhombic/POSCAR-023",
                L"spglibtestdata/orthorhombic/POSCAR-023-2",
                L"spglibtestdata/orthorhombic/POSCAR-024",
                L"spglibtestdata/orthorhombic/POSCAR-024-2",
                L"spglibtestdata/orthorhombic/POSCAR-025",
                L"spglibtestdata/orthorhombic/POSCAR-025-2",
                L"spglibtestdata/orthorhombic/POSCAR-026",
                L"spglibtestdata/orthorhombic/POSCAR-026-2",
                L"spglibtestdata/orthorhombic/POSCAR-027",
                L"spglibtestdata/orthorhombic/POSCAR-027-2",
                L"spglibtestdata/orthorhombic/POSCAR-028",
                L"spglibtestdata/orthorhombic/POSCAR-028-2",
                L"spglibtestdata/orthorhombic/POSCAR-029",
                L"spglibtestdata/orthorhombic/POSCAR-029-2",
                L"spglibtestdata/orthorhombic/POSCAR-030",
                L"spglibtestdata/orthorhombic/POSCAR-030-2",
                L"spglibtestdata/orthorhombic/POSCAR-031",
                L"spglibtestdata/orthorhombic/POSCAR-031-2",
                L"spglibtestdata/orthorhombic/POSCAR-032",
                L"spglibtestdata/orthorhombic/POSCAR-032-2",
                L"spglibtestdata/orthorhombic/POSCAR-033",
                L"spglibtestdata/orthorhombic/POSCAR-033-2",
                L"spglibtestdata/orthorhombic/POSCAR-033-3",
                L"spglibtestdata/orthorhombic/POSCAR-034",
                L"spglibtestdata/orthorhombic/POSCAR-034-2",
                L"spglibtestdata/orthorhombic/POSCAR-035",
                L"spglibtestdata/orthorhombic/POSCAR-035-2",
                L"spglibtestdata/orthorhombic/POSCAR-036",
                L"spglibtestdata/orthorhombic/POSCAR-036-2",
                L"spglibtestdata/orthorhombic/POSCAR-037",
                L"spglibtestdata/orthorhombic/POSCAR-037-2",
                L"spglibtestdata/orthorhombic/POSCAR-038",
                L"spglibtestdata/orthorhombic/POSCAR-038-2",
                L"spglibtestdata/orthorhombic/POSCAR-039",
                L"spglibtestdata/orthorhombic/POSCAR-039-2",
                L"spglibtestdata/orthorhombic/POSCAR-040",
                L"spglibtestdata/orthorhombic/POSCAR-040-2",
                L"spglibtestdata/orthorhombic/POSCAR-041",
                L"spglibtestdata/orthorhombic/POSCAR-041-2",
                L"spglibtestdata/orthorhombic/POSCAR-042",
                L"spglibtestdata/orthorhombic/POSCAR-043",
                L"spglibtestdata/orthorhombic/POSCAR-043-2",
                L"spglibtestdata/orthorhombic/POSCAR-044",
                L"spglibtestdata/orthorhombic/POSCAR-044-2",
                L"spglibtestdata/orthorhombic/POSCAR-045",
                L"spglibtestdata/orthorhombic/POSCAR-045-2",
                L"spglibtestdata/orthorhombic/POSCAR-046",
                L"spglibtestdata/orthorhombic/POSCAR-046-2",
                L"spglibtestdata/orthorhombic/POSCAR-047",
                L"spglibtestdata/orthorhombic/POSCAR-047-2",
                L"spglibtestdata/orthorhombic/POSCAR-048",
                L"spglibtestdata/orthorhombic/POSCAR-048-2",
                L"spglibtestdata/orthorhombic/POSCAR-049",
                L"spglibtestdata/orthorhombic/POSCAR-049-2",
                L"spglibtestdata/orthorhombic/POSCAR-050",
                L"spglibtestdata/orthorhombic/POSCAR-050-2",
                L"spglibtestdata/orthorhombic/POSCAR-051",
                L"spglibtestdata/orthorhombic/POSCAR-051-2",
                L"spglibtestdata/orthorhombic/POSCAR-051-3",
                L"spglibtestdata/orthorhombic/POSCAR-052",
                L"spglibtestdata/orthorhombic/POSCAR-052-2",
                L"spglibtestdata/orthorhombic/POSCAR-053",
                L"spglibtestdata/orthorhombic/POSCAR-053-2",
                L"spglibtestdata/orthorhombic/POSCAR-054",
                L"spglibtestdata/orthorhombic/POSCAR-054-2",
                L"spglibtestdata/orthorhombic/POSCAR-055",
                L"spglibtestdata/orthorhombic/POSCAR-055-2",
                L"spglibtestdata/orthorhombic/POSCAR-056",
                L"spglibtestdata/orthorhombic/POSCAR-056-2",
                L"spglibtestdata/orthorhombic/POSCAR-057",
                L"spglibtestdata/orthorhombic/POSCAR-057-2",
                L"spglibtestdata/orthorhombic/POSCAR-058",
                L"spglibtestdata/orthorhombic/POSCAR-058-2",
                L"spglibtestdata/orthorhombic/POSCAR-058-3",
                L"spglibtestdata/orthorhombic/POSCAR-059",
                L"spglibtestdata/orthorhombic/POSCAR-059-2",
                L"spglibtestdata/orthorhombic/POSCAR-060",
                L"spglibtestdata/orthorhombic/POSCAR-060-2",
                L"spglibtestdata/orthorhombic/POSCAR-060-3",
                L"spglibtestdata/orthorhombic/POSCAR-061",
                L"spglibtestdata/orthorhombic/POSCAR-061-2",
                L"spglibtestdata/orthorhombic/POSCAR-062",
                L"spglibtestdata/orthorhombic/POSCAR-062-2",
                L"spglibtestdata/orthorhombic/POSCAR-063",
                L"spglibtestdata/orthorhombic/POSCAR-063-2",
                L"spglibtestdata/orthorhombic/POSCAR-063-3",
                L"spglibtestdata/orthorhombic/POSCAR-064",
                L"spglibtestdata/orthorhombic/POSCAR-064-2",
                L"spglibtestdata/orthorhombic/POSCAR-064-3",
                L"spglibtestdata/orthorhombic/POSCAR-065",
                L"spglibtestdata/orthorhombic/POSCAR-065-2",
                L"spglibtestdata/orthorhombic/POSCAR-065-3",
                L"spglibtestdata/orthorhombic/POSCAR-066",
                L"spglibtestdata/orthorhombic/POSCAR-066-2",
                L"spglibtestdata/orthorhombic/POSCAR-067",
                L"spglibtestdata/orthorhombic/POSCAR-067-2",
                L"spglibtestdata/orthorhombic/POSCAR-067-3",
                L"spglibtestdata/orthorhombic/POSCAR-068",
                L"spglibtestdata/orthorhombic/POSCAR-068-2",
                L"spglibtestdata/orthorhombic/POSCAR-069",
                L"spglibtestdata/orthorhombic/POSCAR-069-2",
                L"spglibtestdata/orthorhombic/POSCAR-070",
                L"spglibtestdata/orthorhombic/POSCAR-070-2",
                L"spglibtestdata/orthorhombic/POSCAR-071",
                L"spglibtestdata/orthorhombic/POSCAR-071-2",
                L"spglibtestdata/orthorhombic/POSCAR-072",
                L"spglibtestdata/orthorhombic/POSCAR-072-2",
                L"spglibtestdata/orthorhombic/POSCAR-073",
                L"spglibtestdata/orthorhombic/POSCAR-073-2",
                L"spglibtestdata/orthorhombic/POSCAR-074",
                L"spglibtestdata/orthorhombic/POSCAR-074-2"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestTetragonal)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/tetragonal/POSCAR-075",
                L"spglibtestdata/tetragonal/POSCAR-075-2",
                L"spglibtestdata/tetragonal/POSCAR-076",
                L"spglibtestdata/tetragonal/POSCAR-076-2",
                L"spglibtestdata/tetragonal/POSCAR-077",
                L"spglibtestdata/tetragonal/POSCAR-077-2",
                L"spglibtestdata/tetragonal/POSCAR-077-3",
                L"spglibtestdata/tetragonal/POSCAR-078",
                L"spglibtestdata/tetragonal/POSCAR-078-2",
                L"spglibtestdata/tetragonal/POSCAR-079",
                L"spglibtestdata/tetragonal/POSCAR-079-2",
                L"spglibtestdata/tetragonal/POSCAR-080",
                L"spglibtestdata/tetragonal/POSCAR-080-2",
                L"spglibtestdata/tetragonal/POSCAR-081",
                L"spglibtestdata/tetragonal/POSCAR-081-2",
                L"spglibtestdata/tetragonal/POSCAR-082",
                L"spglibtestdata/tetragonal/POSCAR-082-2",
                L"spglibtestdata/tetragonal/POSCAR-083",
                L"spglibtestdata/tetragonal/POSCAR-083-2",
                L"spglibtestdata/tetragonal/POSCAR-083-3",
                L"spglibtestdata/tetragonal/POSCAR-084",
                L"spglibtestdata/tetragonal/POSCAR-084-2",
                L"spglibtestdata/tetragonal/POSCAR-085",
                L"spglibtestdata/tetragonal/POSCAR-085-2",
                L"spglibtestdata/tetragonal/POSCAR-086",
                L"spglibtestdata/tetragonal/POSCAR-086-2",
                L"spglibtestdata/tetragonal/POSCAR-087",
                L"spglibtestdata/tetragonal/POSCAR-087-2",
                L"spglibtestdata/tetragonal/POSCAR-088",
                L"spglibtestdata/tetragonal/POSCAR-088-2",
                L"spglibtestdata/tetragonal/POSCAR-090",
                L"spglibtestdata/tetragonal/POSCAR-090-2",
                L"spglibtestdata/tetragonal/POSCAR-091",
                L"spglibtestdata/tetragonal/POSCAR-091-2",
                L"spglibtestdata/tetragonal/POSCAR-092",
                L"spglibtestdata/tetragonal/POSCAR-092-2",
                L"spglibtestdata/tetragonal/POSCAR-092-3",
                L"spglibtestdata/tetragonal/POSCAR-094",
                L"spglibtestdata/tetragonal/POSCAR-094-2",
                L"spglibtestdata/tetragonal/POSCAR-094-3",
                L"spglibtestdata/tetragonal/POSCAR-095",
                L"spglibtestdata/tetragonal/POSCAR-095-2",
                L"spglibtestdata/tetragonal/POSCAR-096",
                L"spglibtestdata/tetragonal/POSCAR-096-2",
                L"spglibtestdata/tetragonal/POSCAR-097",
                L"spglibtestdata/tetragonal/POSCAR-097-2",
                L"spglibtestdata/tetragonal/POSCAR-098",
                L"spglibtestdata/tetragonal/POSCAR-098-2",
                L"spglibtestdata/tetragonal/POSCAR-099",
                L"spglibtestdata/tetragonal/POSCAR-099-2",
                L"spglibtestdata/tetragonal/POSCAR-100",
                L"spglibtestdata/tetragonal/POSCAR-100-2",
                L"spglibtestdata/tetragonal/POSCAR-102",
                L"spglibtestdata/tetragonal/POSCAR-102-2",
                L"spglibtestdata/tetragonal/POSCAR-103",
                L"spglibtestdata/tetragonal/POSCAR-103-2",
                L"spglibtestdata/tetragonal/POSCAR-104",
                L"spglibtestdata/tetragonal/POSCAR-104-2",
                L"spglibtestdata/tetragonal/POSCAR-105",
                L"spglibtestdata/tetragonal/POSCAR-105-2",
                L"spglibtestdata/tetragonal/POSCAR-106",
                L"spglibtestdata/tetragonal/POSCAR-107",
                L"spglibtestdata/tetragonal/POSCAR-107-2",
                L"spglibtestdata/tetragonal/POSCAR-107-3",
                L"spglibtestdata/tetragonal/POSCAR-108",
                L"spglibtestdata/tetragonal/POSCAR-108-2",
                L"spglibtestdata/tetragonal/POSCAR-109",
                L"spglibtestdata/tetragonal/POSCAR-109-2",
                L"spglibtestdata/tetragonal/POSCAR-110",
                L"spglibtestdata/tetragonal/POSCAR-110-2",
                L"spglibtestdata/tetragonal/POSCAR-111",
                L"spglibtestdata/tetragonal/POSCAR-111-2",
                L"spglibtestdata/tetragonal/POSCAR-112",
                L"spglibtestdata/tetragonal/POSCAR-112-2",
                L"spglibtestdata/tetragonal/POSCAR-113",
                L"spglibtestdata/tetragonal/POSCAR-113-2",
                L"spglibtestdata/tetragonal/POSCAR-114",
                L"spglibtestdata/tetragonal/POSCAR-114-2",
                L"spglibtestdata/tetragonal/POSCAR-115",
                L"spglibtestdata/tetragonal/POSCAR-115-2",
                L"spglibtestdata/tetragonal/POSCAR-115-3",
                L"spglibtestdata/tetragonal/POSCAR-115-4",
                L"spglibtestdata/tetragonal/POSCAR-115-5",
                L"spglibtestdata/tetragonal/POSCAR-116",
                L"spglibtestdata/tetragonal/POSCAR-116-2",
                L"spglibtestdata/tetragonal/POSCAR-117",
                L"spglibtestdata/tetragonal/POSCAR-117-2",
                L"spglibtestdata/tetragonal/POSCAR-118",
                L"spglibtestdata/tetragonal/POSCAR-118-2",
                L"spglibtestdata/tetragonal/POSCAR-119",
                L"spglibtestdata/tetragonal/POSCAR-119-2",
                L"spglibtestdata/tetragonal/POSCAR-120",
                L"spglibtestdata/tetragonal/POSCAR-120-2",
                L"spglibtestdata/tetragonal/POSCAR-121",
                L"spglibtestdata/tetragonal/POSCAR-121-2",
                L"spglibtestdata/tetragonal/POSCAR-122",
                L"spglibtestdata/tetragonal/POSCAR-122-2",
                L"spglibtestdata/tetragonal/POSCAR-122-3",
                L"spglibtestdata/tetragonal/POSCAR-123",
                L"spglibtestdata/tetragonal/POSCAR-123-2",
                L"spglibtestdata/tetragonal/POSCAR-123-3",
                L"spglibtestdata/tetragonal/POSCAR-124",
                L"spglibtestdata/tetragonal/POSCAR-124-2",
                L"spglibtestdata/tetragonal/POSCAR-125",
                L"spglibtestdata/tetragonal/POSCAR-125-2",
                L"spglibtestdata/tetragonal/POSCAR-126",
                L"spglibtestdata/tetragonal/POSCAR-126-2",
                L"spglibtestdata/tetragonal/POSCAR-127",
                L"spglibtestdata/tetragonal/POSCAR-127-2",
                L"spglibtestdata/tetragonal/POSCAR-128",
                L"spglibtestdata/tetragonal/POSCAR-128-2",
                L"spglibtestdata/tetragonal/POSCAR-129",
                L"spglibtestdata/tetragonal/POSCAR-129-2",
                L"spglibtestdata/tetragonal/POSCAR-129-3",
                L"spglibtestdata/tetragonal/POSCAR-130",
                L"spglibtestdata/tetragonal/POSCAR-130-2",
                L"spglibtestdata/tetragonal/POSCAR-131",
                L"spglibtestdata/tetragonal/POSCAR-131-2",
                L"spglibtestdata/tetragonal/POSCAR-132",
                L"spglibtestdata/tetragonal/POSCAR-132-2",
                L"spglibtestdata/tetragonal/POSCAR-133",
                L"spglibtestdata/tetragonal/POSCAR-133-2",
                L"spglibtestdata/tetragonal/POSCAR-134",
                L"spglibtestdata/tetragonal/POSCAR-134-2",
                L"spglibtestdata/tetragonal/POSCAR-135",
                L"spglibtestdata/tetragonal/POSCAR-135-2",
                L"spglibtestdata/tetragonal/POSCAR-136",
                L"spglibtestdata/tetragonal/POSCAR-136-2",
                L"spglibtestdata/tetragonal/POSCAR-136-3",
                L"spglibtestdata/tetragonal/POSCAR-136-4",
                L"spglibtestdata/tetragonal/POSCAR-136-5",
                L"spglibtestdata/tetragonal/POSCAR-137",
                L"spglibtestdata/tetragonal/POSCAR-137-2",
                L"spglibtestdata/tetragonal/POSCAR-137-3",
                L"spglibtestdata/tetragonal/POSCAR-138",
                L"spglibtestdata/tetragonal/POSCAR-138-2",
                L"spglibtestdata/tetragonal/POSCAR-139",
                L"spglibtestdata/tetragonal/POSCAR-139-2",
                L"spglibtestdata/tetragonal/POSCAR-140",
                L"spglibtestdata/tetragonal/POSCAR-140-2",
                L"spglibtestdata/tetragonal/POSCAR-141",
                L"spglibtestdata/tetragonal/POSCAR-141-2",
                L"spglibtestdata/tetragonal/POSCAR-142",
                L"spglibtestdata/tetragonal/POSCAR-142-2",
                L"spglibtestdata/tetragonal/POSCAR-142-3"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestTrigonal)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/trigonal/POSCAR-143",
                L"spglibtestdata/trigonal/POSCAR-143-2",
                L"spglibtestdata/trigonal/POSCAR-144",
                L"spglibtestdata/trigonal/POSCAR-144-2",
                L"spglibtestdata/trigonal/POSCAR-145",
                L"spglibtestdata/trigonal/POSCAR-145-2",
                L"spglibtestdata/trigonal/POSCAR-146",
                L"spglibtestdata/trigonal/POSCAR-146-2",
                L"spglibtestdata/trigonal/POSCAR-147",
                L"spglibtestdata/trigonal/POSCAR-147-2",
                L"spglibtestdata/trigonal/POSCAR-148",
                L"spglibtestdata/trigonal/POSCAR-148-2",
                L"spglibtestdata/trigonal/POSCAR-149",
                L"spglibtestdata/trigonal/POSCAR-149-2",
                L"spglibtestdata/trigonal/POSCAR-150",
                L"spglibtestdata/trigonal/POSCAR-150-2",
                L"spglibtestdata/trigonal/POSCAR-151",
                L"spglibtestdata/trigonal/POSCAR-151-2",
                L"spglibtestdata/trigonal/POSCAR-152",
                L"spglibtestdata/trigonal/POSCAR-152-2",
                L"spglibtestdata/trigonal/POSCAR-153",
                L"spglibtestdata/trigonal/POSCAR-154",
                L"spglibtestdata/trigonal/POSCAR-154-2",
                L"spglibtestdata/trigonal/POSCAR-154-3",
                L"spglibtestdata/trigonal/POSCAR-155",
                L"spglibtestdata/trigonal/POSCAR-155-2",
                L"spglibtestdata/trigonal/POSCAR-156",
                L"spglibtestdata/trigonal/POSCAR-156-2",
                L"spglibtestdata/trigonal/POSCAR-157",
                L"spglibtestdata/trigonal/POSCAR-157-2",
                L"spglibtestdata/trigonal/POSCAR-158",
                L"spglibtestdata/trigonal/POSCAR-158-2",
                L"spglibtestdata/trigonal/POSCAR-159",
                L"spglibtestdata/trigonal/POSCAR-159-2",
                L"spglibtestdata/trigonal/POSCAR-160",
                L"spglibtestdata/trigonal/POSCAR-160-2",
                L"spglibtestdata/trigonal/POSCAR-161",
                L"spglibtestdata/trigonal/POSCAR-161-2",
                L"spglibtestdata/trigonal/POSCAR-162",
                L"spglibtestdata/trigonal/POSCAR-162-2",
                L"spglibtestdata/trigonal/POSCAR-163",
                L"spglibtestdata/trigonal/POSCAR-163-2",
                L"spglibtestdata/trigonal/POSCAR-164",
                L"spglibtestdata/trigonal/POSCAR-164-2",
                L"spglibtestdata/trigonal/POSCAR-165",
                L"spglibtestdata/trigonal/POSCAR-165-2",
                L"spglibtestdata/trigonal/POSCAR-166",
                L"spglibtestdata/trigonal/POSCAR-166-2",
                L"spglibtestdata/trigonal/POSCAR-167",
                L"spglibtestdata/trigonal/POSCAR-167-2",
                L"spglibtestdata/trigonal/POSCAR-167-3"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestHexagonal)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/hexagonal/POSCAR-168",
                L"spglibtestdata/hexagonal/POSCAR-169",
                L"spglibtestdata/hexagonal/POSCAR-169-2",
                L"spglibtestdata/hexagonal/POSCAR-170",
                L"spglibtestdata/hexagonal/POSCAR-170-2",
                L"spglibtestdata/hexagonal/POSCAR-171",
                L"spglibtestdata/hexagonal/POSCAR-171-2",
                L"spglibtestdata/hexagonal/POSCAR-172",
                L"spglibtestdata/hexagonal/POSCAR-173",
                L"spglibtestdata/hexagonal/POSCAR-173-2",
                L"spglibtestdata/hexagonal/POSCAR-174",
                L"spglibtestdata/hexagonal/POSCAR-174-2",
                L"spglibtestdata/hexagonal/POSCAR-175",
                L"spglibtestdata/hexagonal/POSCAR-175-2",
                L"spglibtestdata/hexagonal/POSCAR-176",
                L"spglibtestdata/hexagonal/POSCAR-176-2",
                L"spglibtestdata/hexagonal/POSCAR-177",
                L"spglibtestdata/hexagonal/POSCAR-179",
                L"spglibtestdata/hexagonal/POSCAR-179-2",
                L"spglibtestdata/hexagonal/POSCAR-180",
                L"spglibtestdata/hexagonal/POSCAR-180-2",
                L"spglibtestdata/hexagonal/POSCAR-181",
                L"spglibtestdata/hexagonal/POSCAR-181-2",
                L"spglibtestdata/hexagonal/POSCAR-182",
                L"spglibtestdata/hexagonal/POSCAR-182-2",
                L"spglibtestdata/hexagonal/POSCAR-183",
                L"spglibtestdata/hexagonal/POSCAR-183-2",
                L"spglibtestdata/hexagonal/POSCAR-184",
                L"spglibtestdata/hexagonal/POSCAR-184-2",
                L"spglibtestdata/hexagonal/POSCAR-185",
                L"spglibtestdata/hexagonal/POSCAR-185-2",
                L"spglibtestdata/hexagonal/POSCAR-186",
                L"spglibtestdata/hexagonal/POSCAR-186-2",
                L"spglibtestdata/hexagonal/POSCAR-187",
                L"spglibtestdata/hexagonal/POSCAR-187-2",
                L"spglibtestdata/hexagonal/POSCAR-188",
                L"spglibtestdata/hexagonal/POSCAR-188-2",
                L"spglibtestdata/hexagonal/POSCAR-189",
                L"spglibtestdata/hexagonal/POSCAR-189-2",
                L"spglibtestdata/hexagonal/POSCAR-190",
                L"spglibtestdata/hexagonal/POSCAR-190-2",
                L"spglibtestdata/hexagonal/POSCAR-191",
                L"spglibtestdata/hexagonal/POSCAR-191-2",
                L"spglibtestdata/hexagonal/POSCAR-192",
                L"spglibtestdata/hexagonal/POSCAR-192-2",
                L"spglibtestdata/hexagonal/POSCAR-193",
                L"spglibtestdata/hexagonal/POSCAR-193-2",
                L"spglibtestdata/hexagonal/POSCAR-194",
                L"spglibtestdata/hexagonal/POSCAR-194-2"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestCubic)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/cubic/POSCAR-195",
                L"spglibtestdata/cubic/POSCAR-195-2",
                L"spglibtestdata/cubic/POSCAR-196",
                L"spglibtestdata/cubic/POSCAR-196-2",
                L"spglibtestdata/cubic/POSCAR-197",
                L"spglibtestdata/cubic/POSCAR-197-2",
                L"spglibtestdata/cubic/POSCAR-198",
                L"spglibtestdata/cubic/POSCAR-198-2",
                L"spglibtestdata/cubic/POSCAR-199",
                L"spglibtestdata/cubic/POSCAR-199-2",
                L"spglibtestdata/cubic/POSCAR-200",
                L"spglibtestdata/cubic/POSCAR-200-2",
                L"spglibtestdata/cubic/POSCAR-205",
                L"spglibtestdata/cubic/POSCAR-205-3",
                L"spglibtestdata/cubic/POSCAR-206",
                L"spglibtestdata/cubic/POSCAR-206-2",
                L"spglibtestdata/cubic/POSCAR-207",
                L"spglibtestdata/cubic/POSCAR-208",
                L"spglibtestdata/cubic/POSCAR-208-2",
                L"spglibtestdata/cubic/POSCAR-209",
                L"spglibtestdata/cubic/POSCAR-210",
                L"spglibtestdata/cubic/POSCAR-210-2",
                L"spglibtestdata/cubic/POSCAR-211",
                L"spglibtestdata/cubic/POSCAR-212",
                L"spglibtestdata/cubic/POSCAR-212-2",
                L"spglibtestdata/cubic/POSCAR-213",
                L"spglibtestdata/cubic/POSCAR-213-2",
                L"spglibtestdata/cubic/POSCAR-214",
                L"spglibtestdata/cubic/POSCAR-214-2",
                L"spglibtestdata/cubic/POSCAR-215",
                L"spglibtestdata/cubic/POSCAR-215-2",
                L"spglibtestdata/cubic/POSCAR-216",
                L"spglibtestdata/cubic/POSCAR-216-2",
                L"spglibtestdata/cubic/POSCAR-217",
                L"spglibtestdata/cubic/POSCAR-217-2",
                L"spglibtestdata/cubic/POSCAR-218",
                L"spglibtestdata/cubic/POSCAR-218-2",
                L"spglibtestdata/cubic/POSCAR-219",
                L"spglibtestdata/cubic/POSCAR-219-2",
                L"spglibtestdata/cubic/POSCAR-220",
                L"spglibtestdata/cubic/POSCAR-220-2",
                L"spglibtestdata/cubic/POSCAR-221",
                L"spglibtestdata/cubic/POSCAR-221-2",
                L"spglibtestdata/cubic/POSCAR-222",
                L"spglibtestdata/cubic/POSCAR-222-2",
                L"spglibtestdata/cubic/POSCAR-223",
                L"spglibtestdata/cubic/POSCAR-223-2",
                L"spglibtestdata/cubic/POSCAR-224",
                L"spglibtestdata/cubic/POSCAR-224-2",
                L"spglibtestdata/cubic/POSCAR-225",
                L"spglibtestdata/cubic/POSCAR-225-2",
                L"spglibtestdata/cubic/POSCAR-226",
                L"spglibtestdata/cubic/POSCAR-226-2",
                L"spglibtestdata/cubic/POSCAR-227",
                L"spglibtestdata/cubic/POSCAR-227-2",
                L"spglibtestdata/cubic/POSCAR-228",
                L"spglibtestdata/cubic/POSCAR-228-2",
                L"spglibtestdata/cubic/POSCAR-229",
                L"spglibtestdata/cubic/POSCAR-229-2",
                L"spglibtestdata/cubic/POSCAR-230",
                L"spglibtestdata/cubic/POSCAR-230-2",
                L"spglibtestdata/cubic/POSCAR-230-3",
                L"spglibtestdata/cubic/POSCAR-230-4"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

        TEST_METHOD(TestVirtual)
        {
            std::random_device rd;
            std::mt19937_64 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5, 0.5);

            const std::vector<std::wstring> testData =
            {
                L"spglibtestdata/virtual_structure/POSCAR-1-221-33",
                L"spglibtestdata/virtual_structure/POSCAR-1-222-33",
                L"spglibtestdata/virtual_structure/POSCAR-1-223-33",
                L"spglibtestdata/virtual_structure/POSCAR-1-224-33",
                L"spglibtestdata/virtual_structure/POSCAR-1-227-73",
                L"spglibtestdata/virtual_structure/POSCAR-1-227-93",
                L"spglibtestdata/virtual_structure/POSCAR-1-227-99",
                L"spglibtestdata/virtual_structure/POSCAR-1-230-conv-56",
                L"spglibtestdata/virtual_structure/POSCAR-1-230-prim-33",
                L"spglibtestdata/virtual_structure/POSCAR-1-bcc-33",
                L"spglibtestdata/virtual_structure/POSCAR-10-221-18",
                L"spglibtestdata/virtual_structure/POSCAR-10-223-18",
                L"spglibtestdata/virtual_structure/POSCAR-10-227-50",
                L"spglibtestdata/virtual_structure/POSCAR-102-224-13",
                L"spglibtestdata/virtual_structure/POSCAR-104-222-13",
                L"spglibtestdata/virtual_structure/POSCAR-105-223-13",
                L"spglibtestdata/virtual_structure/POSCAR-109-227-13",
                L"spglibtestdata/virtual_structure/POSCAR-11-227-48",
                L"spglibtestdata/virtual_structure/POSCAR-110-230-conv-15",
                L"spglibtestdata/virtual_structure/POSCAR-110-230-prim-13",
                L"spglibtestdata/virtual_structure/POSCAR-111-221-11",
                L"spglibtestdata/virtual_structure/POSCAR-111-224-11",
                L"spglibtestdata/virtual_structure/POSCAR-111-227-66",
                L"spglibtestdata/virtual_structure/POSCAR-112-222-11",
                L"spglibtestdata/virtual_structure/POSCAR-112-223-11",
                L"spglibtestdata/virtual_structure/POSCAR-113-227-68",
                L"spglibtestdata/virtual_structure/POSCAR-115-221-14",
                L"spglibtestdata/virtual_structure/POSCAR-115-223-14",
                L"spglibtestdata/virtual_structure/POSCAR-115-227-33",
                L"spglibtestdata/virtual_structure/POSCAR-116-230-conv-34",
                L"spglibtestdata/virtual_structure/POSCAR-117-230-conv-33",
                L"spglibtestdata/virtual_structure/POSCAR-118-222-14",
                L"spglibtestdata/virtual_structure/POSCAR-118-224-14",
                L"spglibtestdata/virtual_structure/POSCAR-12-221-19",
                L"spglibtestdata/virtual_structure/POSCAR-12-224-19",
                L"spglibtestdata/virtual_structure/POSCAR-12-227-21",
                L"spglibtestdata/virtual_structure/POSCAR-12-227-83",
                L"spglibtestdata/virtual_structure/POSCAR-120-230-conv-16",
                L"spglibtestdata/virtual_structure/POSCAR-120-230-prim-14",
                L"spglibtestdata/virtual_structure/POSCAR-122-230-conv-13",
                L"spglibtestdata/virtual_structure/POSCAR-122-230-prim-11",
                L"spglibtestdata/virtual_structure/POSCAR-123-221-05",
                L"spglibtestdata/virtual_structure/POSCAR-126-222-05",
                L"spglibtestdata/virtual_structure/POSCAR-13-222-18",
                L"spglibtestdata/virtual_structure/POSCAR-13-224-18",
                L"spglibtestdata/virtual_structure/POSCAR-13-227-49",
                L"spglibtestdata/virtual_structure/POSCAR-13-230-conv-44",
                L"spglibtestdata/virtual_structure/POSCAR-131-223-05",
                L"spglibtestdata/virtual_structure/POSCAR-134-224-05",
                L"spglibtestdata/virtual_structure/POSCAR-14-227-47",
                L"spglibtestdata/virtual_structure/POSCAR-14-227-51",
                L"spglibtestdata/virtual_structure/POSCAR-14-230-conv-45",
                L"spglibtestdata/virtual_structure/POSCAR-142-230-conv-05",
                L"spglibtestdata/virtual_structure/POSCAR-142-230-prim-05",
                L"spglibtestdata/virtual_structure/POSCAR-146-221-27",
                L"spglibtestdata/virtual_structure/POSCAR-146-222-27",
                L"spglibtestdata/virtual_structure/POSCAR-146-223-27",
                L"spglibtestdata/virtual_structure/POSCAR-146-224-27",
                L"spglibtestdata/virtual_structure/POSCAR-146-227-92",
                L"spglibtestdata/virtual_structure/POSCAR-146-230-conv-36",
                L"spglibtestdata/virtual_structure/POSCAR-146-230-conv-55",
                L"spglibtestdata/virtual_structure/POSCAR-146-230-prim-27",
                L"spglibtestdata/virtual_structure/POSCAR-146-bcc-27",
                L"spglibtestdata/virtual_structure/POSCAR-148-221-15",
                L"spglibtestdata/virtual_structure/POSCAR-148-222-15",
                L"spglibtestdata/virtual_structure/POSCAR-148-223-15",
                L"spglibtestdata/virtual_structure/POSCAR-148-224-15",
                L"spglibtestdata/virtual_structure/POSCAR-148-227-70",
                L"spglibtestdata/virtual_structure/POSCAR-148-230-conv-17",
                L"spglibtestdata/virtual_structure/POSCAR-148-230-conv-37",
                L"spglibtestdata/virtual_structure/POSCAR-148-230-prim-15",
                L"spglibtestdata/virtual_structure/POSCAR-148-bcc-15",
                L"spglibtestdata/virtual_structure/POSCAR-15-222-19",
                L"spglibtestdata/virtual_structure/POSCAR-15-223-19",
                L"spglibtestdata/virtual_structure/POSCAR-15-230-conv-21",
                L"spglibtestdata/virtual_structure/POSCAR-15-230-conv-22",
                L"spglibtestdata/virtual_structure/POSCAR-15-230-prim-18",
                L"spglibtestdata/virtual_structure/POSCAR-15-230-prim-19",
                L"spglibtestdata/virtual_structure/POSCAR-15-bcc-18",
                L"spglibtestdata/virtual_structure/POSCAR-15-bcc-19",
                L"spglibtestdata/virtual_structure/POSCAR-155-221-17",
                L"spglibtestdata/virtual_structure/POSCAR-155-222-17",
                L"spglibtestdata/virtual_structure/POSCAR-155-223-17",
                L"spglibtestdata/virtual_structure/POSCAR-155-224-17",
                L"spglibtestdata/virtual_structure/POSCAR-155-227-72",
                L"spglibtestdata/virtual_structure/POSCAR-155-230-conv-19",
                L"spglibtestdata/virtual_structure/POSCAR-155-230-conv-38",
                L"spglibtestdata/virtual_structure/POSCAR-155-230-prim-17",
                L"spglibtestdata/virtual_structure/POSCAR-155-bcc-17",
                L"spglibtestdata/virtual_structure/POSCAR-16-221-20",
                L"spglibtestdata/virtual_structure/POSCAR-16-222-20",
                L"spglibtestdata/virtual_structure/POSCAR-16-223-20",
                L"spglibtestdata/virtual_structure/POSCAR-16-224-20",
                L"spglibtestdata/virtual_structure/POSCAR-16-227-84",
                L"spglibtestdata/virtual_structure/POSCAR-160-221-16",
                L"spglibtestdata/virtual_structure/POSCAR-160-224-16",
                L"spglibtestdata/virtual_structure/POSCAR-160-227-16",
                L"spglibtestdata/virtual_structure/POSCAR-160-227-71",
                L"spglibtestdata/virtual_structure/POSCAR-160-fcc",
                L"spglibtestdata/virtual_structure/POSCAR-161-222-16",
                L"spglibtestdata/virtual_structure/POSCAR-161-223-16",
                L"spglibtestdata/virtual_structure/POSCAR-161-230-conv-18",
                L"spglibtestdata/virtual_structure/POSCAR-161-230-prim-16",
                L"spglibtestdata/virtual_structure/POSCAR-161-bcc-16",
                L"spglibtestdata/virtual_structure/POSCAR-166-221-06",
                L"spglibtestdata/virtual_structure/POSCAR-166-224-06",
                L"spglibtestdata/virtual_structure/POSCAR-166-227-06",
                L"spglibtestdata/virtual_structure/POSCAR-166-227-38",
                L"spglibtestdata/virtual_structure/POSCAR-167-222-06",
                L"spglibtestdata/virtual_structure/POSCAR-167-223-06",
                L"spglibtestdata/virtual_structure/POSCAR-167-230-conv-06",
                L"spglibtestdata/virtual_structure/POSCAR-167-230-prim-06",
                L"spglibtestdata/virtual_structure/POSCAR-167-bcc-6",
                L"spglibtestdata/virtual_structure/POSCAR-17-227-60",
                L"spglibtestdata/virtual_structure/POSCAR-17-227-85",
                L"spglibtestdata/virtual_structure/POSCAR-17-230-conv-46",
                L"spglibtestdata/virtual_structure/POSCAR-18-227-86",
                L"spglibtestdata/virtual_structure/POSCAR-19-227-59",
                L"spglibtestdata/virtual_structure/POSCAR-19-227-89",
                L"spglibtestdata/virtual_structure/POSCAR-19-230-conv-51",
                L"spglibtestdata/virtual_structure/POSCAR-195-221-07",
                L"spglibtestdata/virtual_structure/POSCAR-195-222-07",
                L"spglibtestdata/virtual_structure/POSCAR-195-223-07",
                L"spglibtestdata/virtual_structure/POSCAR-195-224-07",
                L"spglibtestdata/virtual_structure/POSCAR-198-227-40",
                L"spglibtestdata/virtual_structure/POSCAR-198-230-conv-20",
                L"spglibtestdata/virtual_structure/POSCAR-199-230-conv-07",
                L"spglibtestdata/virtual_structure/POSCAR-199-230-prim-07",
                L"spglibtestdata/virtual_structure/POSCAR-2-221-28",
                L"spglibtestdata/virtual_structure/POSCAR-2-222-28",
                L"spglibtestdata/virtual_structure/POSCAR-2-223-28",
                L"spglibtestdata/virtual_structure/POSCAR-2-224-28",
                L"spglibtestdata/virtual_structure/POSCAR-2-227-41",
                L"spglibtestdata/virtual_structure/POSCAR-2-227-74",
                L"spglibtestdata/virtual_structure/POSCAR-2-227-94",
                L"spglibtestdata/virtual_structure/POSCAR-2-230-conv-39",
                L"spglibtestdata/virtual_structure/POSCAR-2-230-conv-57",
                L"spglibtestdata/virtual_structure/POSCAR-2-230-prim-28",
                L"spglibtestdata/virtual_structure/POSCAR-2-bcc-28",
                L"spglibtestdata/virtual_structure/POSCAR-20-227-53",
                L"spglibtestdata/virtual_structure/POSCAR-20-227-90",
                L"spglibtestdata/virtual_structure/POSCAR-20-230-conv-53",
                L"spglibtestdata/virtual_structure/POSCAR-200-221-02",
                L"spglibtestdata/virtual_structure/POSCAR-200-223-02",
                L"spglibtestdata/virtual_structure/POSCAR-201-222-02",
                L"spglibtestdata/virtual_structure/POSCAR-201-224-02",
                L"spglibtestdata/virtual_structure/POSCAR-205-230-conv-08",
                L"spglibtestdata/virtual_structure/POSCAR-206-230-conv-02",
                L"spglibtestdata/virtual_structure/POSCAR-206-230-prim-02",
                L"spglibtestdata/virtual_structure/POSCAR-207-221-04",
                L"spglibtestdata/virtual_structure/POSCAR-207-222-04",
                L"spglibtestdata/virtual_structure/POSCAR-208-223-04",
                L"spglibtestdata/virtual_structure/POSCAR-208-224-04",
                L"spglibtestdata/virtual_structure/POSCAR-21-221-23",
                L"spglibtestdata/virtual_structure/POSCAR-21-222-23",
                L"spglibtestdata/virtual_structure/POSCAR-21-223-23",
                L"spglibtestdata/virtual_structure/POSCAR-21-224-23",
                L"spglibtestdata/virtual_structure/POSCAR-21-230-conv-49",
                L"spglibtestdata/virtual_structure/POSCAR-212-227-19",
                L"spglibtestdata/virtual_structure/POSCAR-213-230-conv-09",
                L"spglibtestdata/virtual_structure/POSCAR-214-230-conv-04",
                L"spglibtestdata/virtual_structure/POSCAR-214-230-prim-04",
                L"spglibtestdata/virtual_structure/POSCAR-215-221-03",
                L"spglibtestdata/virtual_structure/POSCAR-215-224-03",
                L"spglibtestdata/virtual_structure/POSCAR-215-227-18",
                L"spglibtestdata/virtual_structure/POSCAR-216-227-03",
                L"spglibtestdata/virtual_structure/POSCAR-218-222-03",
                L"spglibtestdata/virtual_structure/POSCAR-218-223-03",
                L"spglibtestdata/virtual_structure/POSCAR-22-230-conv-26",
                L"spglibtestdata/virtual_structure/POSCAR-22-230-prim-23",
                L"spglibtestdata/virtual_structure/POSCAR-220-230-conv-03",
                L"spglibtestdata/virtual_structure/POSCAR-220-230-prim-03",
                L"spglibtestdata/virtual_structure/POSCAR-221-221-01",
                L"spglibtestdata/virtual_structure/POSCAR-222-222-01",
                L"spglibtestdata/virtual_structure/POSCAR-223-223-01",
                L"spglibtestdata/virtual_structure/POSCAR-224-224-01",
                L"spglibtestdata/virtual_structure/POSCAR-227-227-01",
                L"spglibtestdata/virtual_structure/POSCAR-230-230-conv-01",
                L"spglibtestdata/virtual_structure/POSCAR-230-230-conv-62",
                L"spglibtestdata/virtual_structure/POSCAR-230-230-prim-01",
                L"spglibtestdata/virtual_structure/POSCAR-24-230-conv-23",
                L"spglibtestdata/virtual_structure/POSCAR-24-230-prim-20",
                L"spglibtestdata/virtual_structure/POSCAR-25-221-21",
                L"spglibtestdata/virtual_structure/POSCAR-25-223-21",
                L"spglibtestdata/virtual_structure/POSCAR-25-227-54",
                L"spglibtestdata/virtual_structure/POSCAR-26-227-64",
                L"spglibtestdata/virtual_structure/POSCAR-27-230-conv-48",
                L"spglibtestdata/virtual_structure/POSCAR-28-227-62",
                L"spglibtestdata/virtual_structure/POSCAR-29-230-conv-52",
                L"spglibtestdata/virtual_structure/POSCAR-3-221-29",
                L"spglibtestdata/virtual_structure/POSCAR-3-222-29",
                L"spglibtestdata/virtual_structure/POSCAR-3-223-29",
                L"spglibtestdata/virtual_structure/POSCAR-3-224-29",
                L"spglibtestdata/virtual_structure/POSCAR-3-227-82",
                L"spglibtestdata/virtual_structure/POSCAR-3-227-95",
                L"spglibtestdata/virtual_structure/POSCAR-3-230-conv-58",
                L"spglibtestdata/virtual_structure/POSCAR-30-227-65",
                L"spglibtestdata/virtual_structure/POSCAR-31-227-58",
                L"spglibtestdata/virtual_structure/POSCAR-32-230-conv-47",
                L"spglibtestdata/virtual_structure/POSCAR-33-227-63",
                L"spglibtestdata/virtual_structure/POSCAR-34-222-21",
                L"spglibtestdata/virtual_structure/POSCAR-34-224-21",
                L"spglibtestdata/virtual_structure/POSCAR-35-221-22",
                L"spglibtestdata/virtual_structure/POSCAR-35-224-22",
                L"spglibtestdata/virtual_structure/POSCAR-35-227-87",
                L"spglibtestdata/virtual_structure/POSCAR-37-222-22",
                L"spglibtestdata/virtual_structure/POSCAR-37-223-22",
                L"spglibtestdata/virtual_structure/POSCAR-38-221-26",
                L"spglibtestdata/virtual_structure/POSCAR-39-224-26",
                L"spglibtestdata/virtual_structure/POSCAR-4-227-77",
                L"spglibtestdata/virtual_structure/POSCAR-4-227-81",
                L"spglibtestdata/virtual_structure/POSCAR-4-227-96",
                L"spglibtestdata/virtual_structure/POSCAR-4-230-conv-59",
                L"spglibtestdata/virtual_structure/POSCAR-40-223-26",
                L"spglibtestdata/virtual_structure/POSCAR-41-222-26",
                L"spglibtestdata/virtual_structure/POSCAR-43-230-conv-25",
                L"spglibtestdata/virtual_structure/POSCAR-43-230-conv-29",
                L"spglibtestdata/virtual_structure/POSCAR-43-230-prim-22",
                L"spglibtestdata/virtual_structure/POSCAR-43-230-prim-26",
                L"spglibtestdata/virtual_structure/POSCAR-43-bcc-22",
                L"spglibtestdata/virtual_structure/POSCAR-43-bcc-26",
                L"spglibtestdata/virtual_structure/POSCAR-44-227-24",
                L"spglibtestdata/virtual_structure/POSCAR-45-230-conv-24",
                L"spglibtestdata/virtual_structure/POSCAR-45-230-prim-21",
                L"spglibtestdata/virtual_structure/POSCAR-46-227-28",
                L"spglibtestdata/virtual_structure/POSCAR-47-221-08",
                L"spglibtestdata/virtual_structure/POSCAR-47-223-08",
                L"spglibtestdata/virtual_structure/POSCAR-48-222-08",
                L"spglibtestdata/virtual_structure/POSCAR-48-224-08",
                L"spglibtestdata/virtual_structure/POSCAR-5-221-32",
                L"spglibtestdata/virtual_structure/POSCAR-5-222-32",
                L"spglibtestdata/virtual_structure/POSCAR-5-223-32",
                L"spglibtestdata/virtual_structure/POSCAR-5-224-32",
                L"spglibtestdata/virtual_structure/POSCAR-5-227-45",
                L"spglibtestdata/virtual_structure/POSCAR-5-227-75",
                L"spglibtestdata/virtual_structure/POSCAR-5-227-98",
                L"spglibtestdata/virtual_structure/POSCAR-5-230-conv-40",
                L"spglibtestdata/virtual_structure/POSCAR-5-230-conv-43",
                L"spglibtestdata/virtual_structure/POSCAR-5-230-conv-61",
                L"spglibtestdata/virtual_structure/POSCAR-5-230-prim-29",
                L"spglibtestdata/virtual_structure/POSCAR-5-230-prim-32",
                L"spglibtestdata/virtual_structure/POSCAR-5-bcc-29",
                L"spglibtestdata/virtual_structure/POSCAR-5-bcc-32",
                L"spglibtestdata/virtual_structure/POSCAR-51-227-29",
                L"spglibtestdata/virtual_structure/POSCAR-53-227-32",
                L"spglibtestdata/virtual_structure/POSCAR-54-230-conv-30",
                L"spglibtestdata/virtual_structure/POSCAR-6-221-30",
                L"spglibtestdata/virtual_structure/POSCAR-6-223-30",
                L"spglibtestdata/virtual_structure/POSCAR-6-227-79",
                L"spglibtestdata/virtual_structure/POSCAR-61-230-conv-31",
                L"spglibtestdata/virtual_structure/POSCAR-62-227-31",
                L"spglibtestdata/virtual_structure/POSCAR-65-221-09",
                L"spglibtestdata/virtual_structure/POSCAR-66-223-09",
                L"spglibtestdata/virtual_structure/POSCAR-67-224-09",
                L"spglibtestdata/virtual_structure/POSCAR-68-222-09",
                L"spglibtestdata/virtual_structure/POSCAR-7-222-30",
                L"spglibtestdata/virtual_structure/POSCAR-7-224-30",
                L"spglibtestdata/virtual_structure/POSCAR-7-227-78",
                L"spglibtestdata/virtual_structure/POSCAR-7-227-80",
                L"spglibtestdata/virtual_structure/POSCAR-7-230-conv-60",
                L"spglibtestdata/virtual_structure/POSCAR-70-230-conv-11",
                L"spglibtestdata/virtual_structure/POSCAR-70-230-prim-09",
                L"spglibtestdata/virtual_structure/POSCAR-70-bcc-9",
                L"spglibtestdata/virtual_structure/POSCAR-73-230-conv-10",
                L"spglibtestdata/virtual_structure/POSCAR-73-230-prim-08",
                L"spglibtestdata/virtual_structure/POSCAR-74-227-09",
                L"spglibtestdata/virtual_structure/POSCAR-75-221-25",
                L"spglibtestdata/virtual_structure/POSCAR-75-222-25",
                L"spglibtestdata/virtual_structure/POSCAR-76-227-61",
                L"spglibtestdata/virtual_structure/POSCAR-77-223-25",
                L"spglibtestdata/virtual_structure/POSCAR-77-224-25",
                L"spglibtestdata/virtual_structure/POSCAR-78-227-91",
                L"spglibtestdata/virtual_structure/POSCAR-78-230-conv-54",
                L"spglibtestdata/virtual_structure/POSCAR-8-221-31",
                L"spglibtestdata/virtual_structure/POSCAR-8-224-31",
                L"spglibtestdata/virtual_structure/POSCAR-8-227-44",
                L"spglibtestdata/virtual_structure/POSCAR-8-227-97",
                L"spglibtestdata/virtual_structure/POSCAR-80-230-conv-28",
                L"spglibtestdata/virtual_structure/POSCAR-80-230-prim-25",
                L"spglibtestdata/virtual_structure/POSCAR-81-221-24",
                L"spglibtestdata/virtual_structure/POSCAR-81-222-24",
                L"spglibtestdata/virtual_structure/POSCAR-81-223-24",
                L"spglibtestdata/virtual_structure/POSCAR-81-224-24",
                L"spglibtestdata/virtual_structure/POSCAR-81-227-88",
                L"spglibtestdata/virtual_structure/POSCAR-81-230-conv-50",
                L"spglibtestdata/virtual_structure/POSCAR-82-230-conv-27",
                L"spglibtestdata/virtual_structure/POSCAR-82-230-prim-24",
                L"spglibtestdata/virtual_structure/POSCAR-83-221-10",
                L"spglibtestdata/virtual_structure/POSCAR-84-223-10",
                L"spglibtestdata/virtual_structure/POSCAR-85-222-10",
                L"spglibtestdata/virtual_structure/POSCAR-86-224-10",
                L"spglibtestdata/virtual_structure/POSCAR-88-230-conv-12",
                L"spglibtestdata/virtual_structure/POSCAR-88-230-prim-10",
                L"spglibtestdata/virtual_structure/POSCAR-89-221-12",
                L"spglibtestdata/virtual_structure/POSCAR-89-222-12",
                L"spglibtestdata/virtual_structure/POSCAR-9-222-31",
                L"spglibtestdata/virtual_structure/POSCAR-9-223-31",
                L"spglibtestdata/virtual_structure/POSCAR-9-227-43",
                L"spglibtestdata/virtual_structure/POSCAR-9-230-conv-41",
                L"spglibtestdata/virtual_structure/POSCAR-9-230-conv-42",
                L"spglibtestdata/virtual_structure/POSCAR-9-230-prim-30",
                L"spglibtestdata/virtual_structure/POSCAR-9-230-prim-31",
                L"spglibtestdata/virtual_structure/POSCAR-9-bcc-30",
                L"spglibtestdata/virtual_structure/POSCAR-9-bcc-31",
                L"spglibtestdata/virtual_structure/POSCAR-91-227-67",
                L"spglibtestdata/virtual_structure/POSCAR-92-227-35",
                L"spglibtestdata/virtual_structure/POSCAR-92-230-conv-35",
                L"spglibtestdata/virtual_structure/POSCAR-93-223-12",
                L"spglibtestdata/virtual_structure/POSCAR-93-224-12",
                L"spglibtestdata/virtual_structure/POSCAR-95-227-36",
                L"spglibtestdata/virtual_structure/POSCAR-95-230-conv-32",
                L"spglibtestdata/virtual_structure/POSCAR-96-227-69",
                L"spglibtestdata/virtual_structure/POSCAR-98-230-conv-14",
                L"spglibtestdata/virtual_structure/POSCAR-98-230-prim-12",
                L"spglibtestdata/virtual_structure/POSCAR-99-221-13"
            };

            for (const std::wstring& fileName : testData)
            {
                std::ifstream t(fileName.c_str());
                std::string fileContent((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

                SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
                parser.startParsing();

                std::vector<std::tuple<double3, int, double> > atoms = parser.firstTestFrame();
                double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
                bool allowPartialOccupancies = false;
                double symmetryPrecision = 1e-5;

                double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
                std::vector<std::tuple<double3, int, double>> randomlyShiftedAtoms{};
                std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
                    [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

                std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

                if (spaceGroup)
                {
                    SKSymmetryCell cell = spaceGroup->cell;
                    double3x3 transformationMatrix = spaceGroup->transformationMatrix;
                    double3x3 rotationMatrix = spaceGroup->rotationMatrix;

                    double3x3 originalUnitCell = rotationMatrix * cell.unitCell() * transformationMatrix;

                    std::wstring failMessage = std::format(L"file: {}, result should be: {}, got {}",
                        fileName, originalUnitCell, unitCell);
                    Assert::IsTrue(originalUnitCell == unitCell, failMessage.c_str());
                }
                else
                {
                    std::wstring failMessage = std::format(L"failed: {}", fileName);
                    Assert::Fail(failMessage.c_str());
                }
            }
        }

    };
}