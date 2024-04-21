#include <gtest/gtest.h>
#include <filesystem>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <tuple>
#include <random>
#include <algorithm>
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

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Triclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
      std::make_pair("spglibtestdata/triclinic/POSCAR-001" , double3x3(double3(4.9159976868, 0.0000000000, 0.0000000000), double3(2.4574988436, 4.2582449067, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.4069974558))),
      std::make_pair("spglibtestdata/triclinic/POSCAR-002" , double3x3(double3(-5.5089974078, -0.0000000000, -0.0000000000), double3(1.7074361851, -2.5382443466, -6.0543179877), double3(2.3104709548, -6.6161727417, -0.0000000000)))
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Monoclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/monoclinic/POSCAR-003", double3x3(double3(-0.0000000000, -4.1293980569, -0.0000000000), double3(-4.1604980423, -0.0000000000, -0.0000000000), double3(1.4636598779, -0.0000000000, -7.2753263256))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-004", double3x3(double3(0.0192362784, -5.0120607273, 0.0000000000), double3(4.3671298014, 2.5060303637, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2140961349))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-004-2", double3x3(double3(0.0000000000, 7.3439965443, 0.0000000000), double3(-4.3146051629, 0.0000000000, 10.9420608705), double3(11.8809944095, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-005", double3x3(double3(0.0000000000, -3.8299981978, 0.0000000000), double3(6.2599970544, 1.9149990989, 0.0000000000), double3(-2.0057067389, 0.0000000000, 6.3612890682))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-005-2", double3x3(double3(-2.4344899442, 0.0000000000, 7.7684668591), double3(-6.4309969739, -5.6024973638, -0.0000000000), double3(6.4309969739, -5.6024973638, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-006", double3x3(double3(-0.2213964213, 6.9674800963, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.6699954499), double3(10.9429948509, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-006-2", double3x3(double3(0.0000000000, 3.2087984901, 0.0000000000), double3(-2.1919524756, 0.0000000000, 6.1584385798), double3(9.3991955773, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-007", double3x3(double3(-3.4242907681, -0.0000000000, -5.8698437366), double3(3.4507059969, -0.0000000000, -5.8698437366), double3(-0.0000000000, -22.5499893893, -0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-007-2", double3x3(double3(0.0000000000, 5.4049974567, 0.0000000000), double3(-2.0448445974, 0.0000000000, 12.9252406329), double3(16.4529922582, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-008", double3x3(double3(8.3249960827, 0.0000000000, 7.0409966869), double3(-8.3249960827, -0.0000000000, 7.0409966869), double3(2.3753536987, -10.6441730009, -0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-008-2", double3x3(double3(-7.0439966855, -4.0688980854, -0.0000000000), double3(7.0439966855, -4.0688980854, 0.0000000000), double3(-2.2889694552, 4.0688980854, 26.6955687405))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-009", double3x3(double3(0.0000000000, -5.6320973499, 0.0000000000), double3(8.1389961703, 2.8160486749, 0.0000000000), double3(-1.2035882678, 2.8160486749, 9.3759078039))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-009-2", double3x3(double3(1.3424000453, -0.0000000000, 9.1237692915), double3(-10.4229950955, -0.0000000000, -0.0000000000), double3(4.5402975251, -9.3434956035, -4.5618846457))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-010", double3x3(double3(-0.0000000000, -3.7769982228, -0.0000000000), double3(-12.3929941686, -0.0000000000, -0.0000000000), double3(5.9123807571, -0.0000000000, -14.2035825069))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-010-2", double3x3(double3(-0.0000000000, -3.7769982228, -0.0000000000), double3(-12.3929941686, -0.0000000000, -0.0000000000), double3(5.9123807571, -0.0000000000, -14.2035825069))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-011", double3x3(double3(0.0000000000, 4.1669980393, 0.0000000000), double3(-4.7272549382, 0.0000000000, 10.0459281057), double3(11.4066946327, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-011-2", double3x3(double3(-0.2256254105, 0.0000000000, 4.8747790476), double3(7.0129967001, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.5389955115, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012", double3x3(double3(2.5087734371, -4.3370210905, -0.0001864450), double3(2.5087734371, 4.3370210905, -0.0001864450), double3(-1.7018015993, 0.0000000000, 4.8033156256))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012-2", double3x3(double3(-2.5086743348, -4.3368272641, 0.0000428228), double3(-2.5086743348, 4.3368272641, 0.0000428228), double3(1.7014607644, -0.0000000000, -4.8030263272))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012-3", double3x3(double3(-6.6449968732, -4.2114980183, -0.0000000000), double3(-6.6449968732, 4.2114980183, -0.0000000000), double3(2.0520515829, -0.0000000000, -10.2230773735))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013", double3x3(double3(-4.8589977136, -0.0000000000, -0.0000000000), double3(0.5498746160, -0.0000000000, -5.8170658220), double3(-0.0000000000, -6.7559968210, -0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013-2", double3x3(double3(-0.0000000000, -7.6279964107, -0.0000000000), double3(-11.5259945765, -0.0000000000, -0.0000000000), double3(4.3627821524, -0.0000000000, -11.2946738742))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013-3", double3x3(double3(0.0000000000, 6.5669969100, 0.0000000000), double3(-0.5056791541, 0.0000000000, 7.9930162785), double3(9.7019954348, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-014", double3x3(double3(-5.0699976144, -0.0000000000, -0.0000000000), double3(-2.2121897783, -0.0000000000, -5.7823347551), double3(-0.0000000000, -13.8299934924, -0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-014-2", double3x3(double3(0.3494222389, -0.0000000000, -7.1444569386), double3(-0.0000000000, -9.9939952974, -0.0000000000), double3(-11.1929947332, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015", double3x3(double3(-5.1896717218, -0.0000000000, 0.0189686342), double3(-2.5948358609, 4.5638432340, 0.0094843171), double3(0.2840727811, -0.0000000000, -10.3538966582))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015-2", double3x3(double3(5.1896717218, 0.0000000000, -0.0189686342), double3(2.5948358609, 4.5638432340, -0.0094843171), double3(-0.2840727811, 0.0000000000, 10.3538966582))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015-3", double3x3(double3(0.0925408600, -0.0000000000, -5.0491496501), double3(-4.7064977854, -5.7609972892, -0.0000000000), double3(-4.7064977854, 5.7609972892, -0.0000000000))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Orthorhombic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-016", double3x3(double3(10.7049949629, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.7339949492, 0.0000000000), double3(0.0000000000, 0.0000000000, 31.6299851168))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-016-2", double3x3(double3(5.6099973603, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.6699973320, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.0499957416))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-017-2", double3x3(double3(0.0000000000, 0.0000000000, 4.3299979626), double3(7.0499966827, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.8499963062, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-018", double3x3(double3(-0.0000000000, -0.0000000000, -5.0639976172), double3(-0.0000000000, -8.3339960785, -0.0000000000), double3(-13.9949934148, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-018-2", double3x3(double3(7.3489965420, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.5149964639, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.8939962855))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-019", double3x3(double3(3.5183598275, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.6304070169, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.3802740222))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-019-2", double3x3(double3(4.8089977372, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.9569967264, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.4659960164))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-020", double3x3(double3(-4.3699979437, -2.5249988119, -0.0000000000), double3(-4.3699979437, 2.5249988119, -0.0000000000), double3(-0.0000000000, -0.0000000000, -8.2399961227))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-021", double3x3(double3(-0.0000000000, -0.0000000000, -3.7999982119), double3(-3.1929984976, -5.2149975461, -0.0000000000), double3(-3.1929984976, 5.2149975461, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-021-2", double3x3(double3(-6.5079969377, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -6.5179969330), double3(-3.2539984689, -7.5819964324, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-022", double3x3(double3(-0.0000000000, -0.0000000000, 5.8307972564), double3(0.0000000000, -6.4444969676, -2.9153986282), double3(6.6689968620, -0.0000000000, 2.9153986282))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-023", double3x3(double3(5.0869976064, 5.0874976061, 5.0869976064), double3(5.0869976064, -5.0874976061, -5.0869976064), double3(-5.0869976064, 5.0874976061, -5.0869976064))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-023-2", double3x3(double3(-0.0000000000, -0.0000000000, -6.0439971560), double3(-8.3459960729, -0.0000000000, -0.0000000000), double3(4.1729980364, 8.8229958484, 3.0219985780))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-024", double3x3(double3(-7.0509966822, -0.0000000000, -0.0000000000), double3(-3.5254983411, -4.9839976548, -3.6424982861), double3(3.5254983411, -4.9839976548, 3.6424982861))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-024-2", double3x3(double3(0.0000000000, 0.0000000000, -12.8849939371), double3(7.9359962658, 7.9609962540, 6.4424969685), double3(7.9359962658, -7.9609962540, -6.4424969685))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-025", double3x3(double3(-2.9189986265, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -3.0659985573), double3(-0.0000000000, -5.6179973565, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-025-2", double3x3(double3(0.0000000000, 5.7675972861, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.8488972478), double3(8.2032961400, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-026", double3x3(double3(0.0000000000, 0.0000000000, 4.0099981131), double3(11.1499947535, 0.0000000000, 0.0000000000), double3(0.0000000000, 11.5099945841, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-026-2", double3x3(double3(-0.0000000000, -0.0000000000, -7.4452964967), double3(-0.0000000000, -8.1759961529, -0.0000000000), double3(-8.6489959303, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-027", double3x3(double3(0.0000000000, 0.0000000000, 9.1609956894), double3(13.0279938698, 0.0000000000, 0.0000000000), double3(0.0000000000, 13.0369938655, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-027-2", double3x3(double3(0.0000000000, 0.0000000000, 8.4167960395), double3(13.7939935093, 0.0000000000, 0.0000000000), double3(0.0000000000, 23.8999887541, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-028", double3x3(double3(0.0000000000, 6.2579970553, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.2029966107), double3(7.9549962568, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-028-2", double3x3(double3(-5.1847975603, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -5.1915975571), double3(-0.0000000000, -6.0979971306, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-029", double3x3(double3(-0.0000000000, -6.7899968050, -0.0000000000), double3(-11.2269947172, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -21.1859900311))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-029-2", double3x3(double3(0.0000000000, 0.0000000000, 5.2589975254), double3(6.4719969547, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.7239949539, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-030", double3x3(double3(0.0000000000, 0.0000000000, 7.6379964060), double3(8.8619958301, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.1869952066, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-030-2", double3x3(double3(4.4809978915, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.6719963900, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.3289932576))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-031", double3x3(double3(4.6539978101, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.9359967363, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.8789958221))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-031-2", double3x3(double3(-0.0000000000, -4.9149976873, -0.0000000000), double3(-5.7429972977, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -9.3609955953))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-032", double3x3(double3(10.3880951120, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.4193950972, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.7006949649))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-032-2", double3x3(double3(-5.8839972313, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -8.2199961321), double3(-0.0000000000, -11.7679944627, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033", double3x3(double3(0.0000000000, 0.0000000000, 4.0556376958), double3(4.1085543821, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.5910714829, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-2", double3x3(double3(-0.0000000000, -4.8139977348, -0.0000000000), double3(-5.4559974327, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -11.7869944537))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-3", double3x3(double3(-6.9989967067, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -9.0939957209), double3(-0.0000000000, -13.8479934839, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-034", double3x3(double3(0.0000000000, 0.0000000000, 5.9199972144), double3(10.8899948758, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.0299943394, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-034-2", double3x3(double3(0.0000000000, 0.0000000000, 7.9459962611), double3(10.3479951308, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.5219950490, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-035", double3x3(double3(-0.0000000000, 3.6199982966, -0.0000000000), double3(-0.0000000000, -0.0000000000, -4.1299980567), double3(-9.6999954357, -1.8099991483, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-035-2", double3x3(double3(-0.0000000000, -0.0000000000, -5.3209974962), double3(-7.1789966220, -8.4124960416, -0.0000000000), double3(-7.1789966220, 8.4124960416, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-036", double3x3(double3(0.0000000000, 0.0000000000, -7.7459963552), double3(0.0000000000, 17.9429915571, 0.0000000000), double3(17.6464916966, 0.0000000000, 3.8729981776))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-036-2", double3x3(double3(-8.6399959345, -4.9899976520, -0.0000000000), double3(-8.6399959345, 4.9899976520, -0.0000000000), double3(-0.0000000000, -0.0000000000, -13.5499936242))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-037", double3x3(double3(-0.0000000000, -0.0000000000, -5.8759972351), double3(-6.0364971596, -9.5114955244, -0.0000000000), double3(-6.0364971596, 9.5114955244, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-037-2", double3x3(double3(0.0000000000, 0.0000000000, 4.7729977541), double3(5.8069972676, 0.0000000000, 0.0000000000), double3(2.9034986338, 7.2909965693, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-038", double3x3(double3(0.0000000000, 0.0000000000, -4.4759978939), double3(0.0000000000, 6.9469967311, 0.0000000000), double3(9.4249955651, 0.0000000000, 2.2379989469))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-038-2", double3x3(double3(4.1479980482, 0.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, -6.7169968394), double3(0.0000000000, 5.9839971843, 3.3584984197))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-039", double3x3(double3(5.4199974497, 0.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, -5.5269973993), double3(0.0000000000, 19.2899909232, 2.7634986997))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-039-2", double3x3(double3(0.0000000000, 0.0000000000, -5.5359973951), double3(0.0000000000, 8.2208461317, 2.7679986975), double3(11.5685945565, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-040", double3x3(double3(0.0000000000, 0.0000000000, -3.5799983155), double3(0.0000000000, 4.9249976826, 1.7899991577), double3(9.5399955110, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-040-2", double3x3(double3(5.0859976068, 0.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, -5.8989972243), double3(0.0000000000, 5.1189975913, 2.9494986121))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-041", double3x3(double3(-0.0000000000, -5.5874973708, -4.5204978729), double3(-0.0000000000, -5.5874973708, 4.5204978729), double3(-11.0619947949, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-041-2", double3x3(double3(-0.0000000000, -3.7899982166, -5.0199976379), double3(-0.0000000000, -3.7899982166, 5.0199976379), double3(-10.2299951864, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-042", double3x3(double3(-2.6559987502, 2.6814987382, -0.0000000000), double3(-2.6559987502, -2.6814987382, -0.0000000000), double3(2.6559987502, 0.0000000000, 5.9344972076))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-043", double3x3(double3(0.0000000000, -4.0784980809, -5.7899972756), double3(0.0000000000, -4.0784980809, 5.7899972756), double3(-19.6469907553, 4.0784980809, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-043-2", double3x3(double3(5.5909973692, 0.0000000000, 5.2864975125), double3(-5.5909973692, 0.0000000000, 5.2864975125), double3(0.0000000000, -11.4364946186, -5.2864975125))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-044", double3x3(double3(-3.6519982816, -0.0000000000, -0.0000000000), double3(-1.8259991408, -2.8259986702, -2.6809987385), double3(1.8259991408, -2.8259986702, 2.6809987385))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-044-2", double3x3(double3(4.3609979480, 0.0000000000, 0.0000000000), double3(-2.1804989740, -7.8899962874, -4.8599977132), double3(-2.1804989740, 7.8899962874, -4.8599977132))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-045", double3x3(double3(0.0000000000, -0.0000000000, -5.5719973781), double3(-11.1029947756, -0.0000000000, -0.0000000000), double3(5.5514973878, 9.4619955477, 2.7859986891))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-045-2", double3x3(double3(-5.9199972144, -0.0000000000, -0.0000000000), double3(2.9599986072, 5.7349973014, 7.0799966686), double3(2.9599986072, 5.7349973014, -7.0799966686))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-046", double3x3(double3(0.0000000000, 5.0899976049, 0.0000000000), double3(11.4199946264, 0.0000000000, 0.0000000000), double3(-5.7099973132, -2.5449988025, -10.9749948358))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-046-2", double3x3(double3(0.0000000000, -0.0000000000, -6.2099970779), double3(3.9899981225, 5.0949976026, 3.1049985390), double3(3.9899981225, -5.0949976026, 3.1049985390))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-047", double3x3(double3(0.0000000000, 0.0000000000, 3.0049985860), double3(3.5849983131, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.8489972478, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-047-2", double3x3(double3(0.0000000000, 0.0000000000, 3.8339981959), double3(7.2899965698, 0.0000000000, 0.0000000000), double3(0.0000000000, 25.2599881141, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-048", double3x3(double3(6.3299970215, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.3299970215, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.5399955110))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-048-2", double3x3(double3(0.0000000000, 4.4792978923, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.0665962043), double3(9.3330956084, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-049", double3x3(double3(0.0000000000, 5.1399975814, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2599961133), double3(9.5599955016, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-049-2", double3x3(double3(-0.0000000000, -3.6769982698, -0.0000000000), double3(-6.2169970746, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -7.7939963326))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-050", double3x3(double3(0.0000000000, 0.0000000000, 4.1149980637), double3(7.3269965523, 0.0000000000, 0.0000000000), double3(0.0000000000, 20.0799905515, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-050-2", double3x3(double3(5.4768974229, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.4768974229, 0.0000000000), double3(0.0000000000, 0.0000000000, 20.7962902145))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051", double3x3(double3(0.0000000000, 3.8359981950, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.9279957990), double3(16.7049921396, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-2", double3x3(double3(-0.0000000000, -0.0000000000, -4.5409978633), double3(-0.0000000000, -8.4099960427, -0.0000000000), double3(-9.3279956108, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-3", double3x3(double3(-2.9659986044, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -4.3699979437), double3(-0.0000000000, -4.5219978722, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-052", double3x3(double3(8.5879959590, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.7659958752, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.3429956037))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-052-2", double3x3(double3(-0.0000000000, -4.8929976976, -0.0000000000), double3(-5.1829975612, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -8.4909960046))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-053", double3x3(double3(8.7939958621, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.0249957534, 0.0000000000), double3(0.0000000000, 0.0000000000, 13.0659938519))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-053-2", double3x3(double3(0.0000000000, 0.0000000000, 3.7299982449), double3(7.3949965203, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.0149962286, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-054", double3x3(double3(-0.0000000000, -5.8099972662, -0.0000000000), double3(-10.1199952381, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -10.9499948476))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-054-2", double3x3(double3(-0.0000000000, -5.7809972798, -0.0000000000), double3(-9.9809953035, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -10.5219950490))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-055", double3x3(double3(0.0000000000, 0.0000000000, 3.9739981301), double3(11.5419945690, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.6899940288, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-055-2", double3x3(double3(-0.0000000000, -0.0000000000, -3.9599981367), double3(-0.0000000000, -7.9149962757, -0.0000000000), double3(-11.2199947205, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-056", double3x3(double3(-4.9109976892, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -5.4119974534), double3(-0.0000000000, -12.4639941352, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-056-2", double3x3(double3(0.0000000000, 7.8029963284, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.5659959693), double3(10.2729951661, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-057", double3x3(double3(-5.2609975245, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -5.7149973109), double3(-0.0000000000, -11.4249946241, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-057-2", double3x3(double3(-6.3589970078, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -6.6669968629), double3(-0.0000000000, -9.7649954052, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058", double3x3(double3(0.0000000000, 0.0000000000, 6.9019967523), double3(10.8399948993, 0.0000000000, 0.0000000000), double3(0.0000000000, 23.6929888515, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-2", double3x3(double3(-0.0000000000, -0.0000000000, -2.9629986058), double3(-0.0000000000, -4.3319979616, -0.0000000000), double3(-4.8729977070, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-3", double3x3(double3(-0.0000000000, -0.0000000000, -2.5724917260), double3(-0.0000000000, -3.8912820645, -0.0000000000), double3(-4.0215465742, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-059", double3x3(double3(-0.0000000000, -0.0000000000, -3.5629983235), double3(-0.0000000000, -4.3689979442, -0.0000000000), double3(-11.5099945841, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-059-2", double3x3(double3(0.0000000000, 2.8649986519, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.0449980967), double3(4.6449978143, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060", double3x3(double3(0.0000000000, 6.2489970596, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.3069970323), double3(16.1519923998, 0.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-2", double3x3(double3(0.0000000000, 0.0000000000, 8.0309962211), double3(9.5179955214, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.7489954127, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-3", double3x3(double3(0.0000000000, 0.0000000000, 3.8726073154), double3(5.3672739848, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.5883789815, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-061", double3x3(double3(0.0000000000, 0.0000000000, 10.5369950419), double3(12.1989942599, 0.0000000000, 0.0000000000), double3(0.0000000000, 13.0469938608, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-061-2", double3x3(double3(5.9929971800, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.8189963208, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.0109962305))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-062", double3x3(double3(-0.0000000000, -6.8979967542, -0.0000000000), double3(-7.4899964756, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -10.9419948513))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-062-2", double3x3(double3(0.0000000000, 0.0000000000, 9.3749955887), double3(9.4949955322, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.1499952240, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063", double3x3(double3(-4.6004978353, -3.5794983157, -0.0000000000), double3(-4.6004978353, 3.5794983157, -0.0000000000), double3(-0.0000000000, -0.0000000000, -9.7709954023))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-2", double3x3(double3(5.5689973796, 0.0000000000, 0.0000000000), double3(2.7844986898, 6.3979969895, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.3199965556))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-3", double3x3(double3(-4.4499979061, -2.8799986448, -0.0000000000), double3(-4.4499979061, 2.8799986448, -0.0000000000), double3(-0.0000000000, -0.0000000000, -13.3299937277))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064", double3x3(double3(-0.0000000000, -5.3699974732, -0.0000000000), double3(-5.4059974563, -0.0000000000, -0.0000000000), double3(-0.0000000000, -2.6849987366, -6.5749969062))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-2", double3x3(double3(5.4649974285, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.4709974257, 0.0000000000), double3(0.0000000000, 2.7354987128, 6.1054971271))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-3", double3x3(double3(-2.8785911754, -2.5705256182, -0.0000000000), double3(-2.8785911754, 2.5705256182, -0.0000000000), double3(-0.0000000000, -0.0000000000, -6.8086137766))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065", double3x3(double3(-2.4189988618, -4.0734980832, -0.0000000000), double3(-2.4189988618, 4.0734980832, -0.0000000000), double3(-0.0000000000, -0.0000000000, -6.1069971264))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-2", double3x3(double3(-0.0000000000, -0.0000000000, -4.0599980896), double3(-2.8799986448, -3.1899984990, -0.0000000000), double3(-2.8799986448, 3.1899984990, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-3", double3x3(double3(-0.0000000000, -3.8712981784, -0.0000000000), double3(-3.9091981606, -0.0000000000, -0.0000000000), double3(0.0511698229, -0.0000000000, -3.9088632501))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-066", double3x3(double3(-3.1649985107, -5.2399975344, -0.0000000000), double3(-3.1649985107, 5.2399975344, -0.0000000000), double3(-0.0000000000, -0.0000000000, -10.5299950452))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-066-2", double3x3(double3(0.0000000000, 0.0000000000, 7.0268966936), double3(7.0696966734, 0.0000000000, 0.0000000000), double3(3.5348483367, 12.7461440024, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067", double3x3(double3(-0.0000000000, -0.0000000000, -4.3739979419), double3(-3.9849981249, -5.8609972422, -0.0000000000), double3(-3.9849981249, 5.8609972422, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-2", double3x3(double3(-5.6829973259, -4.9579976671, -0.0000000000), double3(-5.6829973259, 4.9579976671, -0.0000000000), double3(-0.0000000000, -0.0000000000, -8.2709961081))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-3", double3x3(double3(3.8899981696, 0.0000000000, 3.8199982025), double3(3.8899981696, 0.0000000000, -3.8199982025), double3(0.0000000000, 5.4899974167, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-068", double3x3(double3(-7.4099965133, -0.0000000000, -0.0000000000), double3(-0.0000000000, -0.0000000000, -7.4399964992), double3(-3.7049982566, -11.1299947629, -0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-068-2", double3x3(double3(0.0000000000, 0.0000000000, 6.3839969961), double3(6.4179969801, 0.0000000000, 0.0000000000), double3(3.2089984900, 5.6829973259, 0.0000000000))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-069", double3x3(double3(5.4299974450, 0.0000000000, 3.1949984966), double3(-5.4299974450, 0.0000000000, 3.1949984966), double3(0.0000000000, -6.7999968003, -3.1949984966))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-069-2", double3x3(double3(-0.0000000000, -0.0000000000, -2.7382087116), double3(0.0000000000, 5.6303973507, 1.3691043558), double3(6.2133470764, -0.0000000000, -1.3691043558))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-070", double3x3(double3(4.1779980341, 0.0000000000, 3.5194983439), double3(-4.1779980341, 0.0000000000, 3.5194983439), double3(0.0000000000, -5.0929976035, -3.5194983439))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-070-2", double3x3(double3(0.0000000000, -4.8014977407, -3.7309982444), double3(-0.0000000000, 4.8014977407, -3.7309982444), double3(4.8494977181, -0.0000000000, 3.7309982444))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-071", double3x3(double3(-0.0000000000, 0.0000000000, -2.8749986472), double3(-4.7149977814, -0.0000000000, -0.0000000000), double3(2.3574988907, 7.8534963046, 1.4374993236))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-071-2", double3x3(double3(3.5419983333, 0.0000000000, 0.0000000000), double3(0.0000000000, -0.0000000000, 3.8269981992), double3(-1.7709991667, -6.3479970130, -1.9134990996))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-072", double3x3(double3(0.0000000000, 0.0000000000, 4.8579977141), double3(0.0000000000, 7.5009964705, 0.0000000000), double3(-7.9829962437, -3.7504982352, -2.4289988571))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-072-2", double3x3(double3(-0.0000000000, 0.0000000000, -5.4019974581), double3(-5.9669971923, -0.0000000000, -0.0000000000), double3(2.9834985961, 5.2399975344, 2.7009987291))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-073", double3x3(double3(-8.2701961085, -0.0000000000, -0.0000000000), double3(-0.0000000000, -8.3114960891, -0.0000000000), double3(4.1350980543, 4.1557480445, 10.3034951518))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-073-2", double3x3(double3(-0.0000000000, -5.7499972944, -0.0000000000), double3(-0.0000000000, 0.0000000000, -5.9499972003), double3(10.1149952405, 2.8749986472, 2.9749986001))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-074", double3x3(double3(-0.0000000000, -0.0000000000, -5.6959973198), double3(4.1239980595, 5.7219973076, 2.8479986599), double3(4.1239980595, -5.7219973076, 2.8479986599))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-074-2", double3x3(double3(-5.9119972182, -0.0000000000, -0.0000000000), double3(2.9559986091, 2.9724986013, 4.1939980265), double3(2.9559986091, 2.9724986013, -4.1939980265))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Tetragonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/tetragonal/POSCAR-075", double3x3(double3(0.0000000000, 0.0000000000, 3.9439981442), double3(17.4899917702, 0.0000000000, 0.0000000000), double3(0.0000000000, 17.4899917702, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-075-2", double3x3(double3(9.1969956724, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.1969956724, 0.0000000000), double3(0.0000000000, 0.0000000000, 20.5049903515))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-076", double3x3(double3(3.9809981268, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.9809981268, 0.0000000000), double3(0.0000000000, 0.0000000000, 15.3499927772))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-076-2", double3x3(double3(8.4479960249, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.4479960249, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.9119929833))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077", double3x3(double3(0.0000000000, 0.0000000000, 10.6379949944), double3(11.1639947469, 0.0000000000, 0.0000000000), double3(0.0000000000, 11.1639947469, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077-2", double3x3(double3(7.9799962451, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.9799962451, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.7799953981))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077-3", double3x3(double3(0.0000000000, 0.0000000000, 10.6379949944), double3(11.1639947469, 0.0000000000, 0.0000000000), double3(0.0000000000, 11.1639947469, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-078", double3x3(double3(10.8649948876, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.8649948876, 0.0000000000), double3(0.0000000000, 0.0000000000, 28.3569866568))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-078-2", double3x3(double3(7.6289964102, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.6289964102, 0.0000000000), double3(0.0000000000, 0.0000000000, 29.4969861204))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-079", double3x3(double3(0.0000000000, 0.0000000000, -5.8129972647), double3(4.2419980040, 4.2419980040, 2.9064986324), double3(4.2419980040, -4.2419980040, -2.9064986324))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-079-2", double3x3(double3(0.0000000000, -0.0000000000, -7.6089964196), double3(7.4594964900, 7.4594964900, 3.8044982098), double3(7.4594964900, -7.4594964900, -3.8044982098))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-080", double3x3(double3(0.0000000000, 0.0000000000, -14.6679930981), double3(10.1689952151, 10.1689952151, 7.3339965490), double3(10.1689952151, -10.1689952151, -7.3339965490))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-080-2", double3x3(double3(0.0000000000, 0.0000000000, -5.9849971838), double3(4.8464977195, 4.8464977195, 2.9924985919), double3(4.8464977195, -4.8464977195, -2.9924985919))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-081", double3x3(double3(0.0000000000, 0.0000000000, 6.3199970262), double3(7.6208964140, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.6208964140, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-081-2", double3x3(double3(0.0000000000, 0.0000000000, 5.2949975085), double3(10.1814952092, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.1814952092, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-082", double3x3(double3(-5.5479973894, -0.0000000000, -0.0000000000), double3(-0.0000000000, -5.5479973894, -0.0000000000), double3(2.7739986947, 2.7739986947, 5.0849976073))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-082-2", double3x3(double3(-6.3219970252, -0.0000000000, -0.0000000000), double3(-0.0000000000, -6.3219970252, -0.0000000000), double3(3.1609985126, 3.1609985126, 6.3024970344))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083", double3x3(double3(0.0000000000, 0.0000000000, 3.1344985251), double3(8.3279960813, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.3279960813, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083-2", double3x3(double3(0.0000000000, 0.0000000000, 3.9999981178), double3(12.5999940712, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.5999940712, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083-3", double3x3(double3(-0.000000000, -0.0000000000, -3.6799982684), double3(-3.9099981602, -3.9099981602, -0.0000000000), double3(-3.9099981602, 3.9099981602, -0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-084", double3x3(double3(6.4299969744, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.4299969744, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.6299968803))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-084-2", double3x3(double3(0.0000000000, 0.0000000000, 6.5979968954), double3(7.1669966276, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.1669966276, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-085", double3x3(double3(0.0000000000, 0.0000000000, 4.1009980703), double3(6.2609970539, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.2609970539, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-085-2", double3x3(double3(0.0000000000, 0.0000000000, 7.4899964756), double3(8.3799960569, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.3799960569, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-086", double3x3(double3(0.0000000000, 0.0000000000, 6.1889970878), double3(11.1869947360, 0.0000000000, 0.0000000000), double3(0.0000000000, 11.1869947360, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-086-2", double3x3(double3(0.0000000000, 0.0000000000, 3.8219982016), double3(7.0699966733, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.0699966733, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-087", double3x3(double3(5.7039973160, 5.7039973160, 5.1279975871), double3(5.7039973160, -5.7039973160, -5.1279975871), double3(-5.7039973160, 5.7039973160, -5.1279975871))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-087-2", double3x3(double3(0.0000000000, 0.0000000000, -3.1269985286), double3(4.9424976743, 4.9424976743, 1.5634992643), double3(4.9424976743, -4.9424976743, -1.5634992643))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-088", double3x3(double3(0.0000000000, 0.0000000000, -5.9809971857), double3(6.8479967777, 6.8479967777, 2.9904985928), double3(6.8479967777, -6.8479967777, -2.9904985928))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-088-2", double3x3(double3(-5.7409972986, -0.0000000000, -0.0000000000), double3(-0.0000000000, -5.7409972986, -0.0000000000), double3(2.8704986493, 2.8704986493, 6.5604969130))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-090", double3x3(double3(0.0000000000, 0.0000000000, 7.1599966309), double3(9.5597955017, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.5597955017, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-090-2", double3x3(double3(-3.7149982519, -3.7149982519, -0.0000000000), double3(-3.7149982519, 3.7149982519, -0.0000000000), double3(-0.0000000000, -0.0000000000, -10.6099950076))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-091", double3x3(double3(7.0081967023, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.0081967023, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.6517959290))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-091-2", double3x3(double3(12.4639941352, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.4639941352, 0.0000000000), double3(0.0000000000, 0.0000000000, 26.2229876610))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092", double3x3(double3(7.1039966573, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.1039966573, 0.0000000000), double3(0.0000000000, 0.0000000000, 36.5969827796))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092-2", double3x3(double3(6.5899968991, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.5899968991, 0.0000000000), double3(0.0000000000, 0.0000000000, 17.0399919820))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092-3", double3x3(double3(4.0643712984, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.0643712984, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.6346473255))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094", double3x3(double3(4.6862977949, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.6862977949, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.1909956753))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094-2", double3x3(double3(9.9619953125, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.9619953125, 0.0000000000), double3(0.0000000000, 0.0000000000, 13.4139936882))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094-3", double3x3(double3(7.3449965439, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.3449965439, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.3999951064))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-095", double3x3(double3(6.1699970968, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.1699970968, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.5639959703))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-095-2", double3x3(double3(6.0803971389, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.0803971389, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.3987960480))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-096", double3x3(double3(7.4899964756, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.4899964756, 0.0000000000), double3(0.0000000000, 0.0000000000, 13.2399937700))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-096-2", double3x3(double3(5.9969971782, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.9969971782, 0.0000000000), double3(0.0000000000, 0.0000000000, 39.1909815590))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-097", double3x3(double3(-7.4829964789, -0.0000000000, -0.0000000000), double3(-0.0000000000, -7.4829964789, -0.0000000000), double3(3.7414982395, 3.7414982395, 7.4464964961))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-097-2", double3x3(double3(4.7654977576, 4.7654977576, 6.4119969829), double3(4.7654977576, -4.7654977576, -6.4119969829), double3(-4.7654977576, 4.7654977576, -6.4119969829))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-098", double3x3(double3(0.0000000000, -0.0000000000, -4.6779977988), double3(3.9769981287, 3.9769981287, 2.3389988994), double3(3.9769981287, -3.9769981287, -2.3389988994))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-098-2", double3x3(double3(-9.3829955849, -0.0000000000, -0.0000000000), double3(-0.0000000000, -9.3829955849, -0.0000000000), double3(4.6914977925, 4.6914977925, 27.2999871542))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-099", double3x3(double3(3.8989981654, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.8989981654, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.1669980393))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-099-2", double3x3(double3(3.8071982086, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.8071982086, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.6981977893))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-100", double3x3(double3(0.0000000000, 0.0000000000, 5.2149975461), double3(8.8699958263, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.8699958263, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-100-2", double3x3(double3(8.3119960889, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.3119960889, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0699952616))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-102", double3x3(double3(0.0000000000, 0.0000000000, 8.0199962263), double3(8.8639958291, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.8639958291, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-102-2", double3x3(double3(0.0000000000, 0.0000000000, 5.6559973386), double3(10.7589949374, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.7589949374, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-103", double3x3(double3(6.5509969175, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.5509969175, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.8469967782))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-103-2", double3x3(double3(6.5139969349, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.5139969349, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.8089967961))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-104", double3x3(double3(0.0000000000, 0.0000000000, 9.2369956536), double3(9.4159955694, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.4159955694, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-104-2", double3x3(double3(0.0000000000, 0.0000000000, 9.0089957609), double3(10.6199950028, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.6199950028, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-105", double3x3(double3(7.6179964154, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.6179964154, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.4999960004))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-105-2", double3x3(double3(5.6399973461, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.6399973461, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.0499957416))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-106", double3x3(double3(0.0000000000, 0.0000000000, 5.3079975024), double3(10.8389948998, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.8389948998, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107", double3x3(double3(-7.7619963477, -0.0000000000, -0.0000000000), double3(-0.0000000000, -7.7619963477, -0.0000000000), double3(3.8809981738, 3.8809981738, 5.8729972365))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107-2", double3x3(double3(-4.1199980614, -0.0000000000, -0.0000000000), double3(-0.0000000000, -4.1199980614, -0.0000000000), double3(2.0599990307, 2.0599990307, 5.2374975355))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107-3", double3x3(double3(-3.5518983287, -0.0000000000, -0.0000000000), double3(-0.0000000000, -3.5518983287, -0.0000000000), double3(1.7759491643, 1.7759491643, 12.8226939664))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-108", double3x3(double3(-6.1399971109, -0.0000000000, -0.0000000000), double3(-0.0000000000, -6.1399971109, -0.0000000000), double3(3.0699985554, 3.0699985554, 5.2859975127))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-108-2", double3x3(double3(-8.0549962098, -0.0000000000, -0.0000000000), double3(-0.0000000000, -8.0549962098, -0.0000000000), double3(4.0274981049, 4.0274981049, 7.8439963091))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-109", double3x3(double3(-3.4516983758, -0.0000000000, -0.0000000000), double3(-0.0000000000, -3.4516983758, -0.0000000000), double3(1.7258491879, 1.7258491879, 5.8399972520))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-109-2", double3x3(double3(0.0000000000, 0.0000000000, -8.7499958828), double3(5.6899973226, 5.6899973226, 4.3749979414), double3(5.6899973226, -5.6899973226, -4.3749979414))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-110", double3x3(double3(0.0000000000, 0.0000000000, -9.0999957181), double3(6.8099967956, 6.8099967956, 4.5499978590), double3(6.8099967956, -6.8099967956, -4.5499978590))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-110-2", double3x3(double3(-11.7822944559, -0.0000000000, -0.0000000000), double3(-0.0000000000, -11.7822944559, -0.0000000000), double3(5.8911472280, 5.8911472280, 11.8192444385))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-111", double3x3(double3(5.1599975720, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.1599975720, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0699952616))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-111-2", double3x3(double3(0.0000000000, 0.0000000000, 5.8013972702), double3(5.8150972638, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.8150972638, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-112", double3x3(double3(5.4299974450, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.4299974450, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0959952494))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-112-2", double3x3(double3(5.4149974520, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.4149974520, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.1969952019))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-113", double3x3(double3(0.0000000000, 0.0000000000, 4.7159977809), double3(5.6619973358, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.6619973358, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-113-2", double3x3(double3(0.0000000000, 0.0000000000, 4.2785979867), double3(6.4023969874, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.4023969874, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-114", double3x3(double3(0.0000000000, 0.0000000000, 6.3489970125), double3(7.4839964785, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.4839964785, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-114-2", double3x3(double3(0.0000000000, 0.0000000000, 6.8099967956), double3(10.8099949134, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.8099949134, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115", double3x3(double3(3.3269984345, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.3269984345, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.1509971057))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-2", double3x3(double3(3.8589981842, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.8589981842, 0.0000000000), double3(0.0000000000, 0.0000000000, 11.7099944900))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-3", double3x3(double3(-1.6494930782, -1.6494930782, -0.0000000000), double3(-1.6494930782, 1.6494930782, -0.0000000000), double3(-0.0000000000, -0.0000000000, -4.7886334410))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-4", double3x3(double3(5.0189976384, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.0189976384, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.4279960343))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-5", double3x3(double3(-2.6399987578, -2.6399987578, -0.0000000000), double3(-2.6399987578, 2.6399987578, -0.0000000000), double3(-0.0000000000, -0.0000000000, -5.2049975508))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-116", double3x3(double3(5.5249974003, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.5249974003, 0.0000000000), double3(0.0000000000, 0.0000000000, 17.4629917829))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-116-2", double3x3(double3(10.5659950283, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.5659950283, 0.0000000000), double3(0.0000000000, 0.0000000000, 25.2199881329))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-117", double3x3(double3(-0.0000000000, -0.0000000000, -5.6199973556), double3(-5.4649974285, -5.4649974285, -0.0000000000), double3(-5.4649974285, 5.4649974285, -0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-117-2", double3x3(double3(5.7099973132, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.7099973132, 0.0000000000), double3(0.0000000000, 0.0000000000, 16.4579922558))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-118", double3x3(double3(0.0000000000, 0.0000000000, 5.1646975698), double3(6.5959968963, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.5959968963, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-118-2", double3x3(double3(7.8229963190, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.8229963190, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.0529957402))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-119", double3x3(double3(-6.2679970506, -0.0000000000, -0.0000000000), double3(-0.0000000000, -6.2679970506, -0.0000000000), double3(3.1339985253, 3.1339985253, 7.3909965222))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-119-2", double3x3(double3(-3.9199981555, -0.0000000000, -0.0000000000), double3(-0.0000000000, -3.9199981555, -0.0000000000), double3(1.9599990777, 1.9599990777, 7.6099964192))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-120", double3x3(double3(-8.4669960159, -0.0000000000, -0.0000000000), double3(-0.0000000000, -8.4669960159, -0.0000000000), double3(4.2334980080, 4.2334980080, 6.3724970015))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-120-2", double3x3(double3(6.2237970714, 6.2237970714, 6.4638469585), double3(6.2237970714, -6.2237970714, -6.4638469585), double3(-6.2237970714, 6.2237970714, -6.4638469585))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-121", double3x3(double3(2.4879988293, 2.4879988293, 3.3729984129), double3(2.4879988293, -2.4879988293, -3.3729984129), double3(-2.4879988293, 2.4879988293, -3.3729984129))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-121-2", double3x3(double3(2.5234988126, 2.5234988126, 3.2429984740), double3(2.5234988126, -2.5234988126, -3.2429984740), double3(-2.5234988126, 2.5234988126, -3.2429984740))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122", double3x3(double3(0.0000000000, 0.0000000000, -5.6199973556), double3(5.7419972981, 5.7419972981, 2.8099986778), double3(5.7419972981, -5.7419972981, -2.8099986778))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122-2", double3x3(double3(7.4664964867, 7.4664964867, 7.7734963422), double3(7.4664964867, -7.4664964867, -7.7734963422), double3(-7.4664964867, 7.4664964867, -7.7734963422))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122-3", double3x3(double3(-3.9576830213, -0.0000000000, -0.0000000000), double3(-0.0000000000, -3.9576830213, -0.0000000000), double3(1.9788415106, 1.9788415106, 2.9882372600))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123", double3x3(double3(0.0000000000, 0.0000000000, 3.2789984571), double3(4.0189981089, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.0189981089, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123-2", double3x3(double3(8.9809957741, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.9809957741, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.6379940533))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123-3", double3x3(double3(4.1999980237, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.1999980237, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.9599962545))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-124", double3x3(double3(0.0000000000, 0.0000000000, 4.9959976492), double3(6.1179971212, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.1179971212, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-124-2", double3x3(double3(6.2099970779, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.2099970779, 0.0000000000), double3(0.0000000000, 0.0000000000, 11.0019948231))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-125", double3x3(double3(0.0000000000, 0.0000000000, 6.7129968413), double3(8.5159959929, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.5159959929, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-125-2", double3x3(double3(6.3759969998, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.3759969998, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.3289960809))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-126", double3x3(double3(0.0000000000, 0.0000000000, 6.1899970873), double3(11.3499946594, 0.0000000000, 0.0000000000), double3(0.0000000000, 11.3499946594, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-126-2", double3x3(double3(6.3123970298, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.3123970298, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.5494955066))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-127", double3x3(double3(5.8939972266, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.8939972266, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.3479960719))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-127-2", double3x3(double3(0.0000000000, 0.0000000000, 4.1439980501), double3(7.1049966568, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.1049966568, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-128", double3x3(double3(7.0576966791, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.0576966791, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.9783953047))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-128-2", double3x3(double3(7.7435963563, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.7435963563, 0.0000000000), double3(0.0000000000, 0.0000000000, 11.6402945228))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129", double3x3(double3(4.2819979851, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.2819979851, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.1819970911))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129-2", double3x3(double3(0.0000000000, 0.0000000000, 4.8099977367), double3(5.0049976449, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.0049976449, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129-3", double3x3(double3(4.2819979851, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.2819979851, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.1819970911))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-130", double3x3(double3(5.9239972125, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.9239972125, 0.0000000000), double3(0.0000000000, 0.0000000000, 18.1299914691))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-130-2", double3x3(double3(7.3771965287, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.3771965287, 0.0000000000), double3(0.0000000000, 0.0000000000, 15.1230928839))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-131", double3x3(double3(3.0199985790, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.0199985790, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.3099975014))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-131-2", double3x3(double3(4.9262976820, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.9262976820, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2849961016))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-132", double3x3(double3(0.0000000000, 0.0000000000, 6.0519971523), double3(6.1659970986, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.1659970986, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-132-2", double3x3(double3(5.8519972464, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.8519972464, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.2339933023))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-133", double3x3(double3(0.0000000000, 0.0000000000, 4.6629978059), double3(9.3809955858, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.3809955858, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-133-2", double3x3(double3(11.7889944528, 0.0000000000, 0.0000000000), double3(0.0000000000, 11.7889944528, 0.0000000000), double3(0.0000000000, 0.0000000000, 23.6349888787))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-134", double3x3(double3(8.4269960347, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.4269960347, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.4919931809))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-134-2", double3x3(double3(6.1689970972, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.1689970972, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.2129970765))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-135", double3x3(double3(0.0000000000, 0.0000000000, 5.9129972177), double3(8.5899959580, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.5899959580, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-135-2", double3x3(double3(0.0000000000, 0.0000000000, 5.9419972040), double3(8.5269959877, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.5269959877, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136", double3x3(double3(0.0000000000, 0.0000000000, 2.8729986481), double3(4.3982979304, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.3982979304, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-2", double3x3(double3(0.0000000000, 0.0000000000, 2.9532986103), double3(4.5844978428, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.5844978428, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-3", double3x3(double3(0.0000000000, 0.0000000000, 2.6888359272), double3(4.2266540200, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.2266540200, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-4", double3x3(double3(4.6685978032, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.6685978032, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.2149975461))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-5", double3x3(double3(0.0000000000, 0.0000000000, 4.1099980661), double3(6.7539968220, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.7539968220, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137", double3x3(double3(0.0000000000, 0.0000000000, 5.4499974355), double3(8.0899961933, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.0899961933, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137-2", double3x3(double3(3.6399982872, 0.0000000000, 0.0000000000), double3(0.0000000000, 3.6399982872, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.2699975202))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137-3", double3x3(double3(2.2423146613, 0.0000000000, 0.0000000000), double3(0.0000000000, 2.2423146613, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.8505089202))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-138", double3x3(double3(0.0000000000, 0.0000000000, 7.6781963871), double3(8.4335960316, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.4335960316, 0.0000000000))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-138-2", double3x3(double3(4.3499979531, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.3499979531, 0.0000000000), double3(0.0000000000, 0.0000000000, 13.7299935395))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-139", double3x3(double3(-11.9399943817, -0.0000000000, -0.0000000000), double3(-0.0000000000, -11.9399943817, -0.0000000000), double3(5.9699971909, 5.9699971909, 8.6999959063))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-139-2", double3x3(double3(-4.1699980378, -0.0000000000, -0.0000000000), double3(-0.0000000000, -4.1699980378, -0.0000000000), double3(2.0849990189, 2.0849990189, 5.4399974403))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-140", double3x3(double3(-11.0759947883, -0.0000000000, -0.0000000000), double3(-0.0000000000, -11.0759947883, -0.0000000000), double3(5.5379973941, 5.5379973941, 18.4664913107))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-140-2", double3x3(double3(0.0000000000, -0.0000000000, -5.7269973052), double3(5.5059974092, 5.5059974092, 2.8634986526), double3(5.5059974092, -5.5059974092, -2.8634986526))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-141", double3x3(double3(3.5885983114, 3.5885983114, 3.1644485110), double3(3.5885983114, -3.5885983114, -3.1644485110), double3(-3.5885983114, 3.5885983114, -3.1644485110))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-141-2", double3x3(double3(-6.9012967527, -0.0000000000, -0.0000000000), double3(-0.0000000000, -6.9012967527, -0.0000000000), double3(3.4506483763, 3.4506483763, 9.9876953004))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142", double3x3(double3(-10.3299951393, -0.0000000000, -0.0000000000), double3(-0.0000000000, -10.3299951393, -0.0000000000), double3(5.1649975697, 5.1649975697, 10.1899952052))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142-2", double3x3(double3(-12.2839942199, -0.0000000000, -0.0000000000), double3(-0.0000000000, -12.2839942199, -0.0000000000), double3(6.1419971099, 6.1419971099, 11.7909944518))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142-3", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Trigonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/trigonal/POSCAR-143", double3x3(double3(0.0000000000, 0.0000000000, 6.7835968080), double3(7.2487965891, 0.0000000000, 0.0000000000), double3(-3.6243982946, 6.2776419931, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-143-2", double3x3(double3(0.0000000000, 0.0000000000, 7.3612965362), double3(7.9541862572, 0.0000000000, 0.0000000000), double3(-3.9770931286, 6.8885273652, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-144", double3x3(double3(6.8672967686, 0.0000000000, 0.0000000000), double3(-3.4336483843, 5.9472534570, 0.0000000000), double3(0.0000000000, 0.0000000000, 17.0619919716))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-144-2", double3x3(double3(4.3367979594, 0.0000000000, 0.0000000000), double3(-2.1683989797, 3.7557772039, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.3396960758))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-145", double3x3(double3(12.6919940279, 0.0000000000, 0.0000000000), double3(-6.3459970139, 10.9915892528, 0.0000000000), double3(0.0000000000, 0.0000000000, 19.1859909722))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-145-2", double3x3(double3(0.0000000000, 0.0000000000, 7.4729964836), double3(10.5019950584, 0.0000000000, 0.0000000000), double3(-5.2509975292, 9.0949945110, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-146", double3x3(double3(5.4184974504, 3.1283709616, 2.7189987206), double3(-5.4184974504, 3.1283709616, 2.7189987206), double3(-0.0000000000, -6.2567419231, 2.7189987206))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-146-2", double3x3(double3(-2.9999985884, 1.7320499926, 4.7766644190), double3(-2.9999985884, -1.7320499926, -4.7766644190), double3(0.0000000000, -3.4640999851, 4.7766644190))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-147", double3x3(double3(0.0000000000, 0.0000000000, 8.8887958174), double3(16.9905920052, 0.0000000000, 0.0000000000), double3(-8.4952960026, 14.7142843019, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-147-2", double3x3(double3(0.0000000000, 0.0000000000, 7.2224966015), double3(9.3961955787, 0.0000000000, 0.0000000000), double3(-4.6980977893, 8.1373440701, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-148", double3x3(double3(-3.5533483280, 2.0515266137, 5.6955639867), double3(-3.5533483280, -2.0515266137, -5.6955639867), double3(-0.0000000000, -4.1030532274, 5.6955639867))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-148-2", double3x3(double3(-0.0000000000, 0.0000000000, -12.1339942904), double3(-10.4869950654, 6.0546694240, 4.0446647635), double3(10.4869950654, 6.0546694240, 4.0446647635))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-149", double3x3(double3(5.0219976369, 0.0000000000, 0.0000000000), double3(-2.5109988185, 4.3491775313, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.3759969998))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-149-2", double3x3(double3(7.1509966352, 0.0000000000, 0.0000000000), double3(-3.5754983176, 6.1929447484, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.1757961529))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-150", double3x3(double3(0.0000000000, 0.0000000000, 4.9839976548), double3(9.0699957322, 0.0000000000, 0.0000000000), double3(-4.5349978661, 7.8548467163, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-150-2", double3x3(double3(0.0000000000, 0.0000000000, 4.7379977706), double3(8.6379959355, 0.0000000000, 0.0000000000), double3(-4.3189979677, 7.4807239179, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-151", double3x3(double3(5.9599971956, 0.0000000000, 0.0000000000), double3(-2.9799985978, 5.1615089778, 0.0000000000), double3(0.0000000000, 0.0000000000, 17.1999919067))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-151-2", double3x3(double3(5.0339976313, 0.0000000000, 0.0000000000), double3(-2.5169988156, 4.3595698313, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.1409933461))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-152", double3x3(double3(9.2039956691, 0.0000000000, 0.0000000000), double3(-4.6019978346, 7.9708940658, 0.0000000000), double3(0.0000000000, 0.0000000000, 24.8179883221))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-152-2", double3x3(double3(5.0359976304, 0.0000000000, 0.0000000000), double3(-2.5179988152, 4.3613018813, 0.0000000000), double3(0.0000000000, 0.0000000000, 11.2549947041))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-153", double3x3(double3(6.0199971673, 0.0000000000, 0.0000000000), double3(-3.0099985837, 5.2134704776, 0.0000000000), double3(0.0000000000, 0.0000000000, 17.2999918596))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154", double3x3(double3(4.9133976880, 0.0000000000, 0.0000000000), double3(-2.4566988440, 4.2551272167, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.4051974566))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154-2", double3x3(double3(0.0000000000, 0.0000000000, 8.3999960474), double3(13.0399938641, 0.0000000000, 0.0000000000), double3(-6.5199969321, 11.2929659515, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154-3", double3x3(double3(3.5185386123, -2.0314292150, 0.0000000000), double3(0.0000000000, 4.0628584300, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.6155767045))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-155", double3x3(double3(4.6707978022, 2.6966863684, 2.4351655208), double3(-4.6707978022, 2.6966863684, 2.4351655208), double3(-0.0000000000, -5.3933727369, 2.4351655208))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-155-2", double3x3(double3(-4.5614978536, 2.6335820137, 5.6609973363), double3(-4.5614978536, -2.6335820137, -5.6609973363), double3(-0.0000000000, -5.2671640274, 5.6609973363))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-156", double3x3(double3(3.7332982433, 0.0000000000, 0.0000000000), double3(-1.8666491217, 3.2331311186, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.0979971306))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-156-2", double3x3(double3(3.9149981578, 0.0000000000, 0.0000000000), double3(-1.9574990789, 3.3904878604, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.7249940124))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-157", double3x3(double3(12.1969942608, 0.0000000000, 0.0000000000), double3(-6.0984971304, 10.5629068797, 0.0000000000), double3(0.0000000000, 0.0000000000, 19.3589908908))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-157-2", double3x3(double3(0.0000000000, 0.0000000000, 3.9659981338), double3(8.7529958813, 0.0000000000, 0.0000000000), double3(-4.3764979407, 7.5803167925, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-158", double3x3(double3(0.0000000000, 0.0000000000, 5.6579973377), double3(6.1199971203, 0.0000000000, 0.0000000000), double3(-3.0599985601, 5.3000729773, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-158-2", double3x3(double3(0.0000000000, 0.0000000000, 9.2389956527), double3(12.8739939422, 0.0000000000, 0.0000000000), double3(-6.4369969711, 11.1492058022, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-159", double3x3(double3(10.1999952005, 0.0000000000, 0.0000000000), double3(-5.0999976002, 8.8334549621, 0.0000000000), double3(0.0000000000, 0.0000000000, 30.3509857186))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-159-2", double3x3(double3(0.0000000000, 0.0000000000, 5.3679974741), double3(10.5629950297, 0.0000000000, 0.0000000000), double3(-5.2814975148, 9.1478220357, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-160", double3x3(double3(7.8050963274, 0.0000000000, 0.0000000000), double3(-2.5690227481, 7.3701866190, 0.0000000000), double3(-2.5690227481, -3.6161021795, 6.4221068059))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-160-2", double3x3(double3(-2.7434987091, 1.5839597182, 3.0519985639), double3(-2.7434987091, -1.5839597182, -3.0519985639), double3(-0.0000000000, -3.1679194364, 3.0519985639))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-161", double3x3(double3(5.2189975442, -9.0395689112, -0.0000000000), double3(5.2189975442, 9.0395689112, 0.0000000000), double3(-5.2189975442, 3.0131896371, 12.3833275065))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-161-2", double3x3(double3(2.5799987860, -4.4686889808, -0.0000000000), double3(2.5799987860, 4.4686889808, -0.0000000000), double3(-2.5799987860, 1.4895629936, 5.5266640661))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-162", double3x3(double3(5.4499974355, 0.0000000000, 0.0000000000), double3(-2.7249987178, 4.7198362297, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.1009961881))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-162-2", double3x3(double3(0.0000000000, 0.0000000000, 4.6219978252), double3(4.9899976520, 0.0000000000, 0.0000000000), double3(-2.4949988260, 4.3214647315, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-163", double3x3(double3(5.8899972285, 0.0000000000, 0.0000000000), double3(-2.9449986143, 5.1008872281, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.5909954870))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-163-2", double3x3(double3(5.3099975014, 0.0000000000, 0.0000000000), double3(-2.6549987507, 4.5985927303, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.2499932948))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-164", double3x3(double3(4.0469980957, 0.0000000000, 0.0000000000), double3(-2.0234990479, 3.5048031600, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.3299974920))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-164-2", double3x3(double3(0.0000000000, 0.0000000000, 5.0859976068), double3(6.2489970596, 0.0000000000, 0.0000000000), double3(-3.1244985298, 5.4117902018, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-165", double3x3(double3(7.1849966192, 0.0000000000, 0.0000000000), double3(-3.5924983096, 6.2223895983, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.3509965410))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-165-2", double3x3(double3(0.0000000000, 0.0000000000, 10.1399952287), double3(12.1899942641, 0.0000000000, 0.0000000000), double3(-6.0949971320, 10.5568447047, 0.0000000000))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-166", double3x3(double3(3.1214985312, -5.4065940518, -0.0000000000), double3(3.1214985312, 5.4065940518, -0.0000000000), double3(-3.1214985312, 1.8021980173, 9.9999952946))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-166-2", double3x3(double3(-2.7124987237, 1.5660618683, 3.2786651239), double3(-2.7124987237, -1.5660618683, -3.2786651239), double3(-0.0000000000, -3.1321237366, 3.2786651239))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167", double3x3(double3(5.7549972920, -9.9679477072, -0.0000000000), double3(5.7549972920, 9.9679477072, -0.0000000000), double3(-5.7549972920, 3.3226492357, 20.1899904998))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167-2", double3x3(double3(5.0114976419, -8.6801685377, -0.0000000000), double3(5.0114976419, 8.6801685377, -0.0000000000), double3(-5.0114976419, 2.8933895126, 8.4903293383))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167-3", double3x3(double3(2.4745988356, -4.2861309116, -0.0000000000), double3(2.4745988356, 4.2861309116, -0.0000000000), double3(-2.4745988356, 1.4287103039, 4.6659978045))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Hexagonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/hexagonal/POSCAR-168", double3x3(double3(0.0000000000, 0.0000000000, 3.8919981687), double3(15.9359925014, 0.0000000000, 0.0000000000), double3(-7.9679962507, 13.8009743408, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-169", double3x3(double3(7.1099966544, 0.0000000000, 0.0000000000), double3(-3.5549983272, 6.1574377236, 0.0000000000), double3(0.0000000000, 0.0000000000, 19.3399908997))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-169-2", double3x3(double3(9.7089954315, 0.0000000000, 0.0000000000), double3(-4.8544977158, 8.4082366889, 0.0000000000), double3(0.0000000000, 0.0000000000, 19.3429908983))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-170", double3x3(double3(7.1099966544, 0.0000000000, 0.0000000000), double3(-3.5549983272, 6.1574377236, 0.0000000000), double3(0.0000000000, 0.0000000000, 19.2999909185))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-170-2", double3x3(double3(10.5125950534, 0.0000000000, 0.0000000000), double3(-5.2562975267, 9.1041743759, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.9375929712))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-171", double3x3(double3(0.0000000000, 0.0000000000, 13.0359938660), double3(17.3899918173, 0.0000000000, 0.0000000000), double3(-8.6949959086, 15.0601746854, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-171-2", double3x3(double3(6.3199970262, 0.0000000000, 0.0000000000), double3(-3.1599985131, 5.4732779765, 0.0000000000), double3(0.0000000000, 0.0000000000, 19.2899909232))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-172", double3x3(double3(6.1979970836, 0.0000000000, 0.0000000000), double3(-3.0989985418, 5.3676229270, 0.0000000000), double3(0.0000000000, 0.0000000000, 18.7269911882))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-173", double3x3(double3(7.1329966436, 0.0000000000, 0.0000000000), double3(-3.5664983218, 6.1773562985, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.4139965114))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-173-2", double3x3(double3(0.0000000000, 0.0000000000, 5.2239975419), double3(9.2249956593, 0.0000000000, 0.0000000000), double3(-4.6124978296, 7.9890805907, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-174", double3x3(double3(0.0000000000, 0.0000000000, 3.9874981237), double3(10.2742951655, 0.0000000000, 0.0000000000), double3(-5.1371475828, 8.8978006193, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-174-2", double3x3(double3(0.0000000000, 0.0000000000, 9.8799953510), double3(12.3199942029, 0.0000000000, 0.0000000000), double3(-6.1599971015, 10.6694279542, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-175", double3x3(double3(0.0000000000, 0.0000000000, 9.1380957001), double3(12.6520940467, 0.0000000000, 0.0000000000), double3(-6.3260470233, 10.9570348555, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-175-2", double3x3(double3(5.4589974313, 0.0000000000, 0.0000000000), double3(-2.7294987157, 4.7276304547, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.1789970925))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-176", double3x3(double3(0.0000000000, 0.0000000000, 3.7429982388), double3(6.4179969801, 0.0000000000, 0.0000000000), double3(-3.2089984900, 5.5581484261, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-176-2", double3x3(double3(0.0000000000, 0.0000000000, 9.9499953181), double3(11.6699945088, 0.0000000000, 0.0000000000), double3(-5.8349972544, 10.1065117066, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-177", double3x3(double3(6.3412970162, 0.0000000000, 0.0000000000), double3(-3.1706485081, 5.4917243089, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.4621969593))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-179", double3x3(double3(0.0000000000, 0.0000000000, 6.9071967499), double3(7.2212966021, 0.0000000000, 0.0000000000), double3(-3.6106483010, 6.2538263057, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-179-2", double3x3(double3(10.4119951007, 0.0000000000, 0.0000000000), double3(-5.2059975504, 9.0170522613, 0.0000000000), double3(0.0000000000, 0.0000000000, 15.1839928553))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-180", double3x3(double3(4.8189977325, 0.0000000000, 0.0000000000), double3(-2.4094988662, 4.1733744571, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.5919968982))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-180-2", double3x3(double3(4.8999976943, 0.0000000000, 0.0000000000), double3(-2.4499988472, 4.2435224818, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.3799974685))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-181", double3x3(double3(4.4282979163, 0.0000000000, 0.0000000000), double3(-2.2141489581, 3.8350184910, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.3679970036))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-181-2", double3x3(double3(10.4817950679, 0.0000000000, 0.0000000000), double3(-5.2408975339, 9.0775008060, 0.0000000000), double3(0.0000000000, 0.0000000000, 11.1749947417))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-182", double3x3(double3(5.3099975014, 0.0000000000, 0.0000000000), double3(-2.6549987507, 4.5985927303, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.2499932948))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-182-2", double3x3(double3(5.4579974318, 0.0000000000, 0.0000000000), double3(-2.7289987159, 4.7267644297, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.0159957576))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-183", double3x3(double3(0.0000000000, 0.0000000000, 10.2999951534), double3(19.4999908244, 0.0000000000, 0.0000000000), double3(-9.7499954122, 16.8874874275, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-183-2", double3x3(double3(3.3959984020, 0.0000000000, 0.0000000000), double3(-1.6979992010, 2.9410208874, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.0919976040))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-184", double3x3(double3(0.0000000000, 0.0000000000, 8.4525960227), double3(13.7179935451, 0.0000000000, 0.0000000000), double3(-6.8589967726, 11.8801308990, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-184-2", double3x3(double3(0.0000000000, 0.0000000000, 8.5029959990), double3(13.8019935056, 0.0000000000, 0.0000000000), double3(-6.9009967528, 11.9528769987, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-185", double3x3(double3(9.8849953487, 0.0000000000, 0.0000000000), double3(-4.9424976743, 8.5606570883, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.8049949158))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-185-2", double3x3(double3(6.2599970544, 0.0000000000, 0.0000000000), double3(-3.1299985272, 5.4213164767, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.2489942363))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-186", double3x3(double3(0.0000000000, 0.0000000000, 7.6399964051), double3(9.9799953040, 0.0000000000, 0.0000000000), double3(-4.9899976520, 8.6429294629, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-186-2", double3x3(double3(8.0999961886, 0.0000000000, 0.0000000000), double3(-4.0499980943, 7.0148024699, 0.0000000000), double3(0.0000000000, 0.0000000000, 13.3399937230))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-187", double3x3(double3(0.0000000000, 0.0000000000, 2.8365986653), double3(2.9064986324, 0.0000000000, 0.0000000000), double3(-1.4532493162, 2.5171016517, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-187-2", double3x3(double3(5.4479974365, 0.0000000000, 0.0000000000), double3(-2.7239987182, 4.7181041798, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.0909961928))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-188", double3x3(double3(0.0000000000, 0.0000000000, 5.6579973377), double3(6.1199971203, 0.0000000000, 0.0000000000), double3(-3.0599985601, 5.3000729773, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-188-2", double3x3(double3(9.2179956625, 0.0000000000, 0.0000000000), double3(-4.6089978313, 7.9830184157, 0.0000000000), double3(0.0000000000, 0.0000000000, 18.0419915105))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-189", double3x3(double3(0.0000000000, 0.0000000000, 6.1369971123), double3(8.1539961632, 0.0000000000, 0.0000000000), double3(-4.0769980816, 7.0615678197, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-189-2", double3x3(double3(0.0000000000, 0.0000000000, 3.8569981851), double3(9.6499954593, 0.0000000000, 0.0000000000), double3(-4.8249977296, 8.3571412141, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-190", double3x3(double3(4.5960978373, 0.0000000000, 0.0000000000), double3(-2.2980489187, 3.9803374854, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.9299957981))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-190-2", double3x3(double3(10.5609950306, 0.0000000000, 0.0000000000), double3(-5.2804975153, 9.1460899857, 0.0000000000), double3(0.0000000000, 0.0000000000, 13.5219936373))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-191", double3x3(double3(0.0000000000, 0.0000000000, 3.8439981912), double3(3.9599981367, 0.0000000000, 0.0000000000), double3(-1.9799990683, 3.4294589853, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-191-2", double3x3(double3(0.0000000000, 0.0000000000, 5.8529972459), double3(11.2569947031, 0.0000000000, 0.0000000000), double3(-5.6284973516, 9.7488433832, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-192", double3x3(double3(10.4569950795, 0.0000000000, 0.0000000000), double3(-5.2284975398, 9.0560233861, 0.0000000000), double3(0.0000000000, 0.0000000000, 14.2379933004))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-192-2", double3x3(double3(0.0000000000, 0.0000000000, 9.3407956048), double3(9.7682954036, 0.0000000000, 0.0000000000), double3(-4.8841477018, 8.4595919712, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-193", double3x3(double3(0.0000000000, 0.0000000000, 6.0839971372), double3(8.4899960051, 0.0000000000, 0.0000000000), double3(-4.2449980026, 7.3525522184, 0.0000000000))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-193-2", double3x3(double3(9.7489954127, 0.0000000000, 0.0000000000), double3(-4.8744977063, 8.4428776888, 0.0000000000), double3(0.0000000000, 0.0000000000, 16.4699922502))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-194", double3x3(double3(3.5869983122, 0.0000000000, 0.0000000000), double3(-1.7934991561, 3.1064316617, 0.0000000000), double3(0.0000000000, 0.0000000000, 15.4919927104))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-194-2", double3x3(double3(3.4699983672, 0.0000000000, 0.0000000000), double3(-1.7349991836, 3.0051067371, 0.0000000000), double3(0.0000000000, 0.0000000000, 28.4499866131))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Cubic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/cubic/POSCAR-195", double3x3(double3(10.3499951299, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.3499951299, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.3499951299))),
    std::make_pair("spglibtestdata/cubic/POSCAR-195-2", double3x3(double3(7.2659965810, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.2659965810, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.2659965810))),
    std::make_pair("spglibtestdata/cubic/POSCAR-196", double3x3(double3(6.0769971405, -6.0769971405, 0.0000000000), double3(-6.0769971405, -0.0000000000, -6.0769971405), double3(6.0769971405, 6.0769971405, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-196-2", double3x3(double3(9.3749955887, -9.3749955887, 0.0000000000), double3(-9.3749955887, -0.0000000000, -9.3749955887), double3(9.3749955887, 9.3749955887, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-197", double3x3(double3(5.0726976131, 5.0726976131, 5.0726976131), double3(5.0726976131, -5.0726976131, -5.0726976131), double3(-5.0726976131, 5.0726976131, -5.0726976131))),
    std::make_pair("spglibtestdata/cubic/POSCAR-197-2", double3x3(double3(5.1249975885, 5.1249975885, 5.1249975885), double3(5.1249975885, -5.1249975885, -5.1249975885), double3(-5.1249975885, 5.1249975885, -5.1249975885))),
    std::make_pair("spglibtestdata/cubic/POSCAR-198", double3x3(double3(7.8399963110, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.8399963110, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.8399963110))),
    std::make_pair("spglibtestdata/cubic/POSCAR-198-2", double3x3(double3(12.7529939992, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.7529939992, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.7529939992))),
    std::make_pair("spglibtestdata/cubic/POSCAR-199", double3x3(double3(5.4649974285, 5.4649974285, 5.4649974285), double3(5.4649974285, -5.4649974285, -5.4649974285), double3(-5.4649974285, 5.4649974285, -5.4649974285))),
    std::make_pair("spglibtestdata/cubic/POSCAR-199-2", double3x3(double3(4.2094980193, 4.2094980193, 4.2094980193), double3(4.2094980193, -4.2094980193, -4.2094980193), double3(-4.2094980193, 4.2094980193, -4.2094980193))),
    std::make_pair("spglibtestdata/cubic/POSCAR-200", double3x3(double3(7.4869964771, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.4869964771, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.4869964771))),
    std::make_pair("spglibtestdata/cubic/POSCAR-200-2", double3x3(double3(5.4499974355, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.4499974355, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.4499974355))),
    std::make_pair("spglibtestdata/cubic/POSCAR-205", double3x3(double3(5.6239973537, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.6239973537, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.6239973537))),
    std::make_pair("spglibtestdata/cubic/POSCAR-205-3", double3x3(double3(5.6239973537, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.6239973537, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.6239973537))),
    std::make_pair("spglibtestdata/cubic/POSCAR-206", double3x3(double3(5.4899974167, 5.4899974167, 5.4899974167), double3(5.4899974167, -5.4899974167, -5.4899974167), double3(-5.4899974167, 5.4899974167, -5.4899974167))),
    std::make_pair("spglibtestdata/cubic/POSCAR-206-2", double3x3(double3(5.5149974050, 5.5149974050, 5.5149974050), double3(5.5149974050, -5.5149974050, -5.5149974050), double3(-5.5149974050, 5.5149974050, -5.5149974050))),
    std::make_pair("spglibtestdata/cubic/POSCAR-207", double3x3(double3(4.3999979296, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.3999979296, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.3999979296))),
    std::make_pair("spglibtestdata/cubic/POSCAR-208", double3x3(double3(6.3099970309, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.3099970309, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.3099970309))),
    std::make_pair("spglibtestdata/cubic/POSCAR-208-2", double3x3(double3(9.5429955096, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.5429955096, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.5429955096))),
    std::make_pair("spglibtestdata/cubic/POSCAR-209", double3x3(double3(3.7114982536, -3.7114982536, 0.0000000000), double3(-3.7114982536, -0.0000000000, -3.7114982536), double3(3.7114982536, 3.7114982536, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-210", double3x3(double3(9.9549953158, -9.9549953158, 0.0000000000), double3(-9.9549953158, -0.0000000000, -9.9549953158), double3(9.9549953158, 9.9549953158, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-210-2", double3x3(double3(7.8494963065, -7.8494963065, 0.0000000000), double3(-7.8494963065, -0.0000000000, -7.8494963065), double3(7.8494963065, 7.8494963065, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-211", double3x3(double3(4.8443977205, 4.8443977205, 4.8443977205), double3(4.8443977205, -4.8443977205, -4.8443977205), double3(-4.8443977205, 4.8443977205, -4.8443977205))),
    std::make_pair("spglibtestdata/cubic/POSCAR-212", double3x3(double3(6.7149968403, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.7149968403, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.7149968403))),
    std::make_pair("spglibtestdata/cubic/POSCAR-212-2", double3x3(double3(6.7149968403, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.7149968403, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.7149968403))),
    std::make_pair("spglibtestdata/cubic/POSCAR-213", double3x3(double3(10.2799951628, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.2799951628, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.2799951628))),
    std::make_pair("spglibtestdata/cubic/POSCAR-213-2", double3x3(double3(7.9359962658, 0.0000000000, 0.0000000000), double3(0.0000000000, 7.9359962658, 0.0000000000), double3(0.0000000000, 0.0000000000, 7.9359962658))),
    std::make_pair("spglibtestdata/cubic/POSCAR-214", double3x3(double3(10.8799948805, 10.8799948805, 10.8799948805), double3(10.8799948805, -10.8799948805, -10.8799948805), double3(-10.8799948805, 10.8799948805, -10.8799948805))),
    std::make_pair("spglibtestdata/cubic/POSCAR-214-2", double3x3(double3(6.1574971026, 6.1574971026, 6.1574971026), double3(6.1574971026, -6.1574971026, -6.1574971026), double3(-6.1574971026, 6.1574971026, -6.1574971026))),
    std::make_pair("spglibtestdata/cubic/POSCAR-215", double3x3(double3(5.3929974624, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.3929974624, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.3929974624))),
    std::make_pair("spglibtestdata/cubic/POSCAR-215-2", double3x3(double3(8.3199960851, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.3199960851, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.3199960851))),
    std::make_pair("spglibtestdata/cubic/POSCAR-216", double3x3(double3(0.0000000000, 3.5879983117, -3.5879983117), double3(-3.5879983117, -3.5879983117, -0.0000000000), double3(0.0000000000, 3.5879983117, 3.5879983117))),
    std::make_pair("spglibtestdata/cubic/POSCAR-216-2", double3x3(double3(3.7884982174, -3.7884982174, 0.0000000000), double3(-3.7884982174, -0.0000000000, -3.7884982174), double3(3.7884982174, 3.7884982174, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-217", double3x3(double3(6.3499970121, 6.3499970121, 6.3499970121), double3(6.3499970121, -6.3499970121, -6.3499970121), double3(-6.3499970121, 6.3499970121, -6.3499970121))),
    std::make_pair("spglibtestdata/cubic/POSCAR-217-2", double3x3(double3(5.0839976078, 5.0839976078, 5.0839976078), double3(5.0839976078, -5.0839976078, -5.0839976078), double3(-5.0839976078, 5.0839976078, -5.0839976078))),
    std::make_pair("spglibtestdata/cubic/POSCAR-218", double3x3(double3(8.2939960973, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2939960973, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2939960973))),
    std::make_pair("spglibtestdata/cubic/POSCAR-218-2", double3x3(double3(6.0259971645, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.0259971645, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.0259971645))),
    std::make_pair("spglibtestdata/cubic/POSCAR-219", double3x3(double3(8.6719959195, -8.6719959195, 0.0000000000), double3(-8.6719959195, -0.0000000000, -8.6719959195), double3(8.6719959195, 8.6719959195, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-219-2", double3x3(double3(6.0704971436, -6.0704971436, 0.0000000000), double3(-6.0704971436, -0.0000000000, -6.0704971436), double3(6.0704971436, 6.0704971436, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-220", double3x3(double3(4.9089976901, 4.9089976901, 4.9089976901), double3(4.9089976901, -4.9089976901, -4.9089976901), double3(-4.9089976901, 4.9089976901, -4.9089976901))),
    std::make_pair("spglibtestdata/cubic/POSCAR-220-2", double3x3(double3(4.2669979922, 4.2669979922, 4.2669979922), double3(4.2669979922, -4.2669979922, -4.2669979922), double3(-4.2669979922, 4.2669979922, -4.2669979922))),
    std::make_pair("spglibtestdata/cubic/POSCAR-221", double3x3(double3(9.6379954649, 0.0000000000, 0.0000000000), double3(0.0000000000, 9.6379954649, 0.0000000000), double3(0.0000000000, 0.0000000000, 9.6379954649))),
    std::make_pair("spglibtestdata/cubic/POSCAR-221-2", double3x3(double3(5.7949972732, 0.0000000000, 0.0000000000), double3(0.0000000000, 5.7949972732, 0.0000000000), double3(0.0000000000, 0.0000000000, 5.7949972732))),
    std::make_pair("spglibtestdata/cubic/POSCAR-222", double3x3(double3(10.9899948287, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.9899948287, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.9899948287))),
    std::make_pair("spglibtestdata/cubic/POSCAR-222-2", double3x3(double3(16.2559923509, 0.0000000000, 0.0000000000), double3(0.0000000000, 16.2559923509, 0.0000000000), double3(0.0000000000, 0.0000000000, 16.2559923509))),
    std::make_pair("spglibtestdata/cubic/POSCAR-223", double3x3(double3(6.6699968615, 0.0000000000, 0.0000000000), double3(0.0000000000, 6.6699968615, 0.0000000000), double3(0.0000000000, 0.0000000000, 6.6699968615))),
    std::make_pair("spglibtestdata/cubic/POSCAR-223-2", double3x3(double3(10.2999951534, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.2999951534, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.2999951534))),
    std::make_pair("spglibtestdata/cubic/POSCAR-224", double3x3(double3(4.9039976925, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.9039976925, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.9039976925))),
    std::make_pair("spglibtestdata/cubic/POSCAR-224-2", double3x3(double3(4.9039976925, 0.0000000000, 0.0000000000), double3(0.0000000000, 4.9039976925, 0.0000000000), double3(0.0000000000, 0.0000000000, 4.9039976925))),
    std::make_pair("spglibtestdata/cubic/POSCAR-225", double3x3(double3(4.9949976496, -4.9949976496, 0.0000000000), double3(-4.9949976496, -0.0000000000, -4.9949976496), double3(4.9949976496, 4.9949976496, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-225-2", double3x3(double3(4.0964980724, -4.0964980724, 0.0000000000), double3(-4.0964980724, -0.0000000000, -4.0964980724), double3(4.0964980724, 4.0964980724, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-226", double3x3(double3(12.5299941041, -12.5299941041, 0.0000000000), double3(-12.5299941041, -0.0000000000, -12.5299941041), double3(12.5299941041, 12.5299941041, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-226-2", double3x3(double3(5.0229976365, -5.0229976365, 0.0000000000), double3(-5.0229976365, -0.0000000000, -5.0229976365), double3(5.0229976365, 5.0229976365, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-227", double3x3(double3(-5.0649976167, 0.0000000000, 5.0649976167), double3(-0.0000000000, -5.0649976167, -5.0649976167), double3(5.0649976167, 0.0000000000, 5.0649976167))),
    std::make_pair("spglibtestdata/cubic/POSCAR-227-2", double3x3(double3(11.6274945288, -11.6274945288, 0.0000000000), double3(-11.6274945288, -0.0000000000, -11.6274945288), double3(11.6274945288, 11.6274945288, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-228", double3x3(double3(7.8524963051, -7.8524963051, 0.0000000000), double3(-7.8524963051, -0.0000000000, -7.8524963051), double3(7.8524963051, 7.8524963051, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-228-2", double3x3(double3(10.9049948687, -10.9049948687, 0.0000000000), double3(-10.9049948687, -0.0000000000, -10.9049948687), double3(10.9049948687, 10.9049948687, 0.0000000000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-229", double3x3(double3(9.1349957016, 9.1349957016, 9.1349957016), double3(9.1349957016, -9.1349957016, -9.1349957016), double3(-9.1349957016, 9.1349957016, -9.1349957016))),
    std::make_pair("spglibtestdata/cubic/POSCAR-229-2", double3x3(double3(3.1104985364, 3.1104985364, 3.1104985364), double3(3.1104985364, -3.1104985364, -3.1104985364), double3(-3.1104985364, 3.1104985364, -3.1104985364))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230", double3x3(double3(6.3009970351, 6.3009970351, 6.3009970351), double3(6.3009970351, -6.3009970351, -6.3009970351), double3(-6.3009970351, 6.3009970351, -6.3009970351))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-2", double3x3(double3(6.1879970883, 6.1879970883, 6.1879970883), double3(6.1879970883, -6.1879970883, -6.1879970883), double3(-6.1879970883, 6.1879970883, -6.1879970883))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-3", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-4", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindSmallestPrimitiveCellNoPartialOccupancies, Virtual)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-221-33", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-222-33", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-223-33", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-224-33", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-73", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-93", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-99", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-conv-56", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-prim-33", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-bcc-33", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-221-18", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-223-18", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-227-50", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-102-224-13", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-104-222-13", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-105-223-13", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-109-227-13", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-11-227-48", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-conv-15", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-prim-13", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-221-11", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-224-11", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-227-66", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-222-11", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-223-11", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-113-227-68", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-221-14", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-223-14", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-227-33", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-116-230-conv-34", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-117-230-conv-33", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-222-14", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-224-14", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-221-19", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-224-19", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-21", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-83", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-conv-16", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-prim-14", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-conv-13", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-prim-11", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-123-221-05", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-126-222-05", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-222-18", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-224-18", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-227-49", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-230-conv-44", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-131-223-05", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-134-224-05", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-47", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-51", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-230-conv-45", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-conv-05", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-prim-05", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-221-27", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-222-27", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-223-27", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-224-27", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-227-92", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-36", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-55", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-prim-27", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-bcc-27", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-221-15", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-222-15", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-223-15", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-224-15", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-227-70", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-17", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-37", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-prim-15", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-bcc-15", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-222-19", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-223-19", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-21", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-22", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-18", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-19", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-18", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-19", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-221-17", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-222-17", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-223-17", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-224-17", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-227-72", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-19", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-38", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-prim-17", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-bcc-17", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-221-20", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-222-20", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-223-20", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-224-20", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-227-84", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-221-16", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-224-16", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-16", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-71", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-fcc", double3x3(double3(5.0000000000, -5.0000000000, 0.0000000000), double3(-5.0000000000, -0.0000000000, -5.0000000000), double3(5.0000000000, 5.0000000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-222-16", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-223-16", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-conv-18", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-prim-16", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-bcc-16", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-221-06", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-224-06", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-06", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-38", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-222-06", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-223-06", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-conv-06", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-prim-06", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-bcc-6", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-60", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-85", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-230-conv-46", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-18-227-86", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-59", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-89", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-230-conv-51", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-221-07", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-222-07", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-223-07", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-224-07", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-227-40", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-230-conv-20", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-conv-07", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-prim-07", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-221-28", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-222-28", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-223-28", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-224-28", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-41", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-74", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-94", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-39", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-57", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-prim-28", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-bcc-28", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-53", double3x3(double3(-0.0000000000, -4.1200000000, -4.1200000000), double3(-0.0000000000, -4.1200000000, 4.1200000000), double3(-8.2400000000, -0.0000000000, -0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-90", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-230-conv-53", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-221-02", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-223-02", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-222-02", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-224-02", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-205-230-conv-08", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-conv-02", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-prim-02", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-221-04", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-222-04", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-223-04", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-224-04", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-221-23", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-222-23", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-223-23", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-224-23", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-230-conv-49", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-212-227-19", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-213-230-conv-09", double3x3(double3(12.8065400000, 0.0000000000, 0.0000000000), double3(0.0000000000, 12.8065400000, 0.0000000000), double3(0.0000000000, 0.0000000000, 12.8065400000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-conv-04", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-prim-04", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-221-03", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-224-03", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-227-18", double3x3(double3(8.2400000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 8.2400000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 8.2400000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-216-227-03", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-222-03", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-223-03", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-conv-26", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-prim-23", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-conv-03", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-prim-03", double3x3(double3(-6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, 6.4032700000), double3(6.4032700000, 6.4032700000, -6.4032700000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-221-221-01", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-222-222-01", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-223-223-01", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-224-224-01", double3x3(double3(10.0000000000, 0.0000000000, 0.0000000000), double3(0.0000000000, 10.0000000000, 0.0000000000), double3(0.0000000000, 0.0000000000, 10.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-227-227-01", double3x3(double3(4.1200000000, -4.1200000000, 0.0000000000), double3(-4.1200000000, -0.0000000000, -4.1200000000), double3(4.1200000000, 4.1200000000, 0.0000000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-01", double3x3(double3(6.4032700000, 6.4032700000, 6.4032700000), double3(6.4032700000, -6.4032700000, -6.4032700000), double3(-6.4032700000, 6.4032700000, -6.4032700000))) 
  };

  for (auto& [fileName, unitCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();

    std::map<size_t, size_t> histogram{};
    for (const std::tuple<double3, size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<size_t, size_t>::iterator index = std::min_element(histogram.begin(), histogram.end(),
      [](const auto& l, const auto& r) { return l.second < r.second; });
    size_t leastOccuringAtomType = index->first;

    std::vector<std::tuple<double3, size_t, double> > reducedAtoms{};
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms), [leastOccuringAtomType](std::tuple<double3, size_t, double> a) {return std::get<1>(a) == leastOccuringAtomType; });
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, int, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, randomlyShiftedAtoms, unitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<double3x3> DelaunayCell = SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

    if (DelaunayCell)
    {
      EXPECT_EQ(*DelaunayCell, unitCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}
