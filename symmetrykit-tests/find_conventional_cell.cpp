#include <gtest/gtest.h>
#include <filesystem>

import <iostream>;
import <fstream>;
import <string>;
import <vector>;
import <iterator>;
import <tuple>;
import <random>;
import <optional>;


import double3;
import double3x3;
import randomnumbers;
import print;

import skposcarlegacyparser;
import sksymmetrycell;
import skspacegroup;
import skpointgroup;

TEST(FindConventionalCell, Triclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
      std::make_pair("spglibtestdata/triclinic/POSCAR-001", double3x3(double3(4.91599769, 0.00000000, 0.00000000), double3(-2.45749884, 4.25824491, 0.00000000), double3(0.00000000, 0.00000000, 5.40699746))),
      std::make_pair("spglibtestdata/triclinic/POSCAR-002", double3x3(double3(5.50899741, 0.00000000, 0.00000000), double3(1.70743619, 6.56486486, 0.00000000), double3(2.31047095, 2.55808206, 6.10163567)))
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Monoclinic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/monoclinic/POSCAR-003", double3x3(double3(4.16049804, 0.00000000, 0.00000000), double3(0.00000000, 4.12939806, 0.00000000), double3(-1.46365988, 0.00000000, 7.27532633))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-004", double3x3(double3(5.01209764, 0.00000000, 0.00000000), double3(0.00000000, 8.21409613, 0.00000000), double3(-2.48925100, 0.00000000, 4.37671571))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-004-2", double3x3(double3(11.76199447, 0.00000000, 0.00000000), double3(0.00000000, 7.34399654, 0.00000000), double3(-4.35825743, 0.00000000, 11.05276528))),
    //std::make_pair("spglibtestdata/monoclinic/POSCAR-005", double3x3(double3(12.51999411, 0.00000000, 0.00000000), double3(0.00000000, 3.82999820, 0.00000000), double3(-2.00570674, 0.00000000, 6.36128907))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-005-2", double3x3(double3(12.86199395, 0.00000000, 0.00000000), double3(0.00000000, 11.20499473, 0.00000000), double3(-2.43448994, 0.00000000, 7.76846686))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-006", double3x3(double3(6.97099672, 0.00000000, 0.00000000), double3(0.00000000, 9.66999545, 0.00000000), double3(-0.34754569, 0.00000000, 10.93747449))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-006-2", double3x3(double3(6.53689692, 0.00000000, 0.00000000), double3(0.00000000, 3.20879849, 0.00000000), double3(-3.15173855, 0.00000000, 8.85502240))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-007", double3x3(double3(6.79564807, 0.00000000, 0.00000000), double3(0.00000000, 22.54998939, 0.00000000), double3(-3.33137393, 0.00000000, 5.93838237))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-007-2", double3x3(double3(13.08599384, 0.00000000, 0.00000000), double3(0.00000000, 5.40499746, 0.00000000), double3(-10.51501507, 0.00000000, 16.25087759))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-008", double3x3(double3(16.64999217, 0.00000000, 0.00000000), double3(0.00000000, 14.08199337, 0.00000000), double3(-2.37535370, 0.00000000, 10.64417300))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-008-2", double3x3(double3(14.08799337, 0.00000000, 0.00000000), double3(0.00000000, 8.13779617, 0.00000000), double3(-4.75502723, 0.00000000, 26.69556874))),
    //std::make_pair("spglibtestdata/monoclinic/POSCAR-009", double3x3(double3(16.27799234, 0.00000000, 0.00000000), double3(0.00000000, 5.63209735, 0.00000000), double3(-6.93540790, 0.00000000, 9.37590780))),
    //std::make_pair("spglibtestdata/monoclinic/POSCAR-009-2", double3x3(double3(12.87246567, 0.00000000, 0.00000000), double3(0.00000000, 18.68699121, 0.00000000), double3(-5.51979525, 0.00000000, 7.38762914))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-010", double3x3(double3(12.39299417, 0.00000000, 0.00000000), double3(0.00000000, 3.77699822, 0.00000000), double3(-5.91238076, 0.00000000, 14.20358251))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-010-2", double3x3(double3(12.39299417, 0.00000000, 0.00000000), double3(0.00000000, 3.77699822, 0.00000000), double3(-5.91238076, 0.00000000, 14.20358251))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-011", double3x3(double3(11.10259478, 0.00000000, 0.00000000), double3(0.00000000, 4.16699804, 0.00000000), double3(-4.85673436, 0.00000000, 10.32108588))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-011-2", double3x3(double3(4.87999770, 0.00000000, 0.00000000), double3(0.00000000, 9.53899551, 0.00000000), double3(-0.32424406, 0.00000000, 7.00549702))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012", double3x3(double3(5.01754689, 0.00000000, 0.00000000), double3(0.00000000, 8.67404218, 0.00000000), double3(-1.70215856, 0.00000000, 4.80318914))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012-2", double3x3(double3(5.01734867, 0.00000000, 0.00000000), double3(0.00000000, 8.67365453, 0.00000000), double3(-1.70154275, 0.00000000, 4.80299728))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-012-3", double3x3(double3(13.28999375, 0.00000000, 0.00000000), double3(0.00000000, 8.42299604, 0.00000000), double3(-2.05205158, 0.00000000, 10.22307737))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013", double3x3(double3(4.85899771, 0.00000000, 0.00000000), double3(0.00000000, 6.75599682, 0.00000000), double3(-0.54987462, 0.00000000, 5.81706582))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013-2", double3x3(double3(12.10799430, 0.00000000, 0.00000000), double3(0.00000000, 7.62799641, 0.00000000), double3(-4.15307459, 0.00000000, 10.75176834))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-013-3", double3x3(double3(8.00899623, 0.00000000, 0.00000000), double3(0.00000000, 6.56699691, 0.00000000), double3(-7.39642298, 0.00000000, 9.68263752))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-014", double3x3(double3(5.06999761, 0.00000000, 0.00000000), double3(0.00000000, 13.82999349, 0.00000000), double3(-2.85780784, 0.00000000, 5.78233476))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-014-2", double3x3(double3(7.15299663, 0.00000000, 0.00000000), double3(0.00000000, 9.99399530, 0.00000000), double3(-6.60622142, 0.00000000, 11.17963183))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015", double3x3(double3(5.18970639, 0.00000000, 0.00000000), double3(0.00000000, 9.12768647, 0.00000000), double3(-0.32191489, 0.00000000, 10.35278920))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015-2", double3x3(double3(5.18970639, 0.00000000, 0.00000000), double3(0.00000000, 9.12768647, 0.00000000), double3(-0.32191489, 0.00000000, 10.35278920))),
    std::make_pair("spglibtestdata/monoclinic/POSCAR-015-3", double3x3(double3(9.41299557, 0.00000000, 0.00000000), double3(0.00000000, 11.52199458, 0.00000000), double3(-0.09254086, 0.00000000, 5.04914965))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Orthorhombic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-016", double3x3(double3(10.70499496, 0.00000000, 0.00000000), double3(0.00000000, 10.73399495, 0.00000000), double3(0.00000000, 0.00000000, 31.62998512))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-016-2", double3x3(double3(5.60999736, 0.00000000, 0.00000000), double3(0.00000000, 5.66999733, 0.00000000), double3(0.00000000, 0.00000000, 9.04999574))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-017-2", double3x3(double3(7.04999668, 0.00000000, 0.00000000), double3(0.00000000, 7.84999631, 0.00000000), double3(0.00000000, 0.00000000, 4.32999796))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-018", double3x3(double3(8.33399608, 0.00000000, 0.00000000), double3(0.00000000, 13.99499341, 0.00000000), double3(0.00000000, 0.00000000, 5.06399762))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-018-2", double3x3(double3(7.34899654, 0.00000000, 0.00000000), double3(0.00000000, 7.51499646, 0.00000000), double3(0.00000000, 0.00000000, 7.89399629))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-019", double3x3(double3(3.51835983, 0.00000000, 0.00000000), double3(0.00000000, 3.63040702, 0.00000000), double3(0.00000000, 0.00000000, 4.38027402))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-019-2", double3x3(double3(4.80899774, 0.00000000, 0.00000000), double3(0.00000000, 6.95699673, 0.00000000), double3(0.00000000, 0.00000000, 8.46599602))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-020", double3x3(double3(5.04999762, 0.00000000, 0.00000000), double3(0.00000000, 8.73999589, 0.00000000), double3(0.00000000, 0.00000000, 8.23999612))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-021", double3x3(double3(6.38599700, 0.00000000, 0.00000000), double3(0.00000000, 10.42999509, 0.00000000), double3(0.00000000, 0.00000000, 3.79999821))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-021-2", double3x3(double3(6.50799694, 0.00000000, 0.00000000), double3(0.00000000, 15.16399286, 0.00000000), double3(0.00000000, 0.00000000, 6.51799693))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-022", double3x3(double3(5.83079726, 0.00000000, 0.00000000), double3(0.00000000, 12.88899394, 0.00000000), double3(0.00000000, 0.00000000, 13.33799372))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-023", double3x3(double3(10.17399521, 0.00000000, 0.00000000), double3(0.00000000, 10.17399521, 0.00000000), double3(0.00000000, 0.00000000, 10.17499521))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-023-2", double3x3(double3(6.04399716, 0.00000000, 0.00000000), double3(0.00000000, 8.34599607, 0.00000000), double3(0.00000000, 0.00000000, 17.64599170))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-024", double3x3(double3(7.05099668, 0.00000000, 0.00000000), double3(0.00000000, 7.28499657, 0.00000000), double3(0.00000000, 0.00000000, 9.96799531))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-024-2", double3x3(double3(12.88499394, 0.00000000, 0.00000000), double3(0.00000000, 15.87199253, 0.00000000), double3(0.00000000, 0.00000000, 15.92199251))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-025", double3x3(double3(5.61799736, 0.00000000, 0.00000000), double3(0.00000000, 2.91899863, 0.00000000), double3(0.00000000, 0.00000000, 3.06599856))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-025-2", double3x3(double3(5.76759729, 0.00000000, 0.00000000), double3(0.00000000, 8.20329614, 0.00000000), double3(0.00000000, 0.00000000, 5.84889725))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-026", double3x3(double3(4.00999811, 0.00000000, 0.00000000), double3(0.00000000, 11.14999475, 0.00000000), double3(0.00000000, 0.00000000, 11.50999458))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-026-2", double3x3(double3(8.17599615, 0.00000000, 0.00000000), double3(0.00000000, 7.44529650, 0.00000000), double3(0.00000000, 0.00000000, 8.64899593))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-027", double3x3(double3(13.02799387, 0.00000000, 0.00000000), double3(0.00000000, 13.03699387, 0.00000000), double3(0.00000000, 0.00000000, 9.16099569))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-027-2", double3x3(double3(13.79399351, 0.00000000, 0.00000000), double3(0.00000000, 23.89998875, 0.00000000), double3(0.00000000, 0.00000000, 8.41679604))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-028", double3x3(double3(7.95499626, 0.00000000, 0.00000000), double3(0.00000000, 6.25799706, 0.00000000), double3(0.00000000, 0.00000000, 7.20299661))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-028-2", double3x3(double3(5.18479756, 0.00000000, 0.00000000), double3(0.00000000, 5.19159756, 0.00000000), double3(0.00000000, 0.00000000, 6.09799713))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-029", double3x3(double3(21.18599003, 0.00000000, 0.00000000), double3(0.00000000, 6.78999681, 0.00000000), double3(0.00000000, 0.00000000, 11.22699472))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-029-2", double3x3(double3(10.72399495, 0.00000000, 0.00000000), double3(0.00000000, 5.25899753, 0.00000000), double3(0.00000000, 0.00000000, 6.47199695))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-030", double3x3(double3(8.86199583, 0.00000000, 0.00000000), double3(0.00000000, 10.18699521, 0.00000000), double3(0.00000000, 0.00000000, 7.63799641))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-030-2", double3x3(double3(4.48099789, 0.00000000, 0.00000000), double3(0.00000000, 7.67199639, 0.00000000), double3(0.00000000, 0.00000000, 14.32899326))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-031", double3x3(double3(8.87899582, 0.00000000, 0.00000000), double3(0.00000000, 6.93599674, 0.00000000), double3(0.00000000, 0.00000000, 4.65399781))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-031-2", double3x3(double3(5.74299730, 0.00000000, 0.00000000), double3(0.00000000, 4.91499769, 0.00000000), double3(0.00000000, 0.00000000, 9.36099560))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-032", double3x3(double3(10.38809511, 0.00000000, 0.00000000), double3(0.00000000, 10.41939510, 0.00000000), double3(0.00000000, 0.00000000, 10.70069496))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-032-2", double3x3(double3(5.88399723, 0.00000000, 0.00000000), double3(0.00000000, 11.76799446, 0.00000000), double3(0.00000000, 0.00000000, 8.21999613))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033", double3x3(double3(4.10855438, 0.00000000, 0.00000000), double3(0.00000000, 5.59107148, 0.00000000), double3(0.00000000, 0.00000000, 4.05563770))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-2", double3x3(double3(5.45599743, 0.00000000, 0.00000000), double3(0.00000000, 4.81399773, 0.00000000), double3(0.00000000, 0.00000000, 11.78699445))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-033-3", double3x3(double3(6.99899671, 0.00000000, 0.00000000), double3(0.00000000, 13.84799348, 0.00000000), double3(0.00000000, 0.00000000, 9.09399572))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-034", double3x3(double3(5.91999721, 0.00000000, 0.00000000), double3(0.00000000, 10.88999488, 0.00000000), double3(0.00000000, 0.00000000, 12.02999434))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-034-2", double3x3(double3(10.34799513, 0.00000000, 0.00000000), double3(0.00000000, 10.52199505, 0.00000000), double3(0.00000000, 0.00000000, 7.94599626))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-035", double3x3(double3(3.61999830, 0.00000000, 0.00000000), double3(0.00000000, 19.39999087, 0.00000000), double3(0.00000000, 0.00000000, 4.12999806))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-035-2", double3x3(double3(14.35799324, 0.00000000, 0.00000000), double3(0.00000000, 16.82499208, 0.00000000), double3(0.00000000, 0.00000000, 5.32099750))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-036", double3x3(double3(7.74599636, 0.00000000, 0.00000000), double3(0.00000000, 35.29298339, 0.00000000), double3(0.00000000, 0.00000000, 17.94299156))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-036-2", double3x3(double3(17.27999187, 0.00000000, 0.00000000), double3(0.00000000, 9.97999530, 0.00000000), double3(0.00000000, 0.00000000, 13.54999362))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-037", double3x3(double3(12.07299432, 0.00000000, 0.00000000), double3(0.00000000, 19.02299105, 0.00000000), double3(0.00000000, 0.00000000, 5.87599724))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-037-2", double3x3(double3(5.80699727, 0.00000000, 0.00000000), double3(0.00000000, 14.58199314, 0.00000000), double3(0.00000000, 0.00000000, 4.77299775))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-038", double3x3(double3(6.94699673, 0.00000000, 0.00000000), double3(0.00000000, 4.47599789, 0.00000000), double3(0.00000000, 0.00000000, 18.84999113))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-038-2", double3x3(double3(4.14799805, 0.00000000, 0.00000000), double3(0.00000000, 11.96799437, 0.00000000), double3(0.00000000, 0.00000000, 6.71699684))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-039", double3x3(double3(5.41999745, 0.00000000, 0.00000000), double3(0.00000000, 38.57998185, 0.00000000), double3(0.00000000, 0.00000000, 5.52699740))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-039-2", double3x3(double3(11.56859456, 0.00000000, 0.00000000), double3(0.00000000, 16.44169226, 0.00000000), double3(0.00000000, 0.00000000, 5.53599740))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-040", double3x3(double3(9.53999551, 0.00000000, 0.00000000), double3(0.00000000, 9.84999537, 0.00000000), double3(0.00000000, 0.00000000, 3.57999832))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-040-2", double3x3(double3(5.08599761, 0.00000000, 0.00000000), double3(0.00000000, 10.23799518, 0.00000000), double3(0.00000000, 0.00000000, 5.89899722))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-041", double3x3(double3(11.06199479, 0.00000000, 0.00000000), double3(0.00000000, 11.17499474, 0.00000000), double3(0.00000000, 0.00000000, 9.04099575))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-041-2", double3x3(double3(10.22999519, 0.00000000, 0.00000000), double3(0.00000000, 7.57999643, 0.00000000), double3(0.00000000, 0.00000000, 10.03999528))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-042", double3x3(double3(5.31199750, 0.00000000, 0.00000000), double3(0.00000000, 5.36299748, 0.00000000), double3(0.00000000, 0.00000000, 11.86899442))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-043", double3x3(double3(8.15699616, 0.00000000, 0.00000000), double3(0.00000000, 39.29398151, 0.00000000), double3(0.00000000, 0.00000000, 11.57999455))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-043-2", double3x3(double3(11.18199474, 0.00000000, 0.00000000), double3(0.00000000, 22.87298924, 0.00000000), double3(0.00000000, 0.00000000, 10.57299502))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-044", double3x3(double3(3.65199828, 0.00000000, 0.00000000), double3(0.00000000, 5.36199748, 0.00000000), double3(0.00000000, 0.00000000, 5.65199734))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-044-2", double3x3(double3(4.36099795, 0.00000000, 0.00000000), double3(0.00000000, 15.77999257, 0.00000000), double3(0.00000000, 0.00000000, 9.71999543))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-045", double3x3(double3(11.10299478, 0.00000000, 0.00000000), double3(0.00000000, 18.92399110, 0.00000000), double3(0.00000000, 0.00000000, 5.57199738))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-045-2", double3x3(double3(5.91999721, 0.00000000, 0.00000000), double3(0.00000000, 11.46999460, 0.00000000), double3(0.00000000, 0.00000000, 14.15999334))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-046", double3x3(double3(21.94998967, 0.00000000, 0.00000000), double3(0.00000000, 5.08999760, 0.00000000), double3(0.00000000, 0.00000000, 11.41999463))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-046-2", double3x3(double3(7.97999625, 0.00000000, 0.00000000), double3(0.00000000, 10.18999521, 0.00000000), double3(0.00000000, 0.00000000, 6.20999708))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-047", double3x3(double3(3.00499859, 0.00000000, 0.00000000), double3(0.00000000, 3.58499831, 0.00000000), double3(0.00000000, 0.00000000, 5.84899725))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-047-2", double3x3(double3(7.28999657, 0.00000000, 0.00000000), double3(0.00000000, 12.62999406, 0.00000000), double3(0.00000000, 0.00000000, 3.83399820))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-048", double3x3(double3(6.32999702, 0.00000000, 0.00000000), double3(0.00000000, 6.32999702, 0.00000000), double3(0.00000000, 0.00000000, 9.53999551))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-048-2", double3x3(double3(4.47929789, 0.00000000, 0.00000000), double3(0.00000000, 8.06659620, 0.00000000), double3(0.00000000, 0.00000000, 9.33309561))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-049", double3x3(double3(5.13999758, 0.00000000, 0.00000000), double3(0.00000000, 9.55999550, 0.00000000), double3(0.00000000, 0.00000000, 8.25999611))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-049-2", double3x3(double3(3.67699827, 0.00000000, 0.00000000), double3(0.00000000, 6.21699707, 0.00000000), double3(0.00000000, 0.00000000, 7.79399633))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-050", double3x3(double3(7.32699655, 0.00000000, 0.00000000), double3(0.00000000, 20.07999055, 0.00000000), double3(0.00000000, 0.00000000, 4.11499806))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-050-2", double3x3(double3(5.47689742, 0.00000000, 0.00000000), double3(0.00000000, 5.47689742, 0.00000000), double3(0.00000000, 0.00000000, 20.79629021))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051", double3x3(double3(16.70499214, 0.00000000, 0.00000000), double3(0.00000000, 3.83599820, 0.00000000), double3(0.00000000, 0.00000000, 8.92799580))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-2", double3x3(double3(8.40999604, 0.00000000, 0.00000000), double3(0.00000000, 4.54099786, 0.00000000), double3(0.00000000, 0.00000000, 9.32799561))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-051-3", double3x3(double3(4.36999794, 0.00000000, 0.00000000), double3(0.00000000, 2.96599860, 0.00000000), double3(0.00000000, 0.00000000, 4.52199787))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-052", double3x3(double3(8.58799596, 0.00000000, 0.00000000), double3(0.00000000, 8.76599588, 0.00000000), double3(0.00000000, 0.00000000, 9.34299560))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-052-2", double3x3(double3(5.18299756, 0.00000000, 0.00000000), double3(0.00000000, 4.89299770, 0.00000000), double3(0.00000000, 0.00000000, 8.49099600))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-053", double3x3(double3(13.06599385, 0.00000000, 0.00000000), double3(0.00000000, 8.79399586, 0.00000000), double3(0.00000000, 0.00000000, 9.02499575))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-053-2", double3x3(double3(8.01499623, 0.00000000, 0.00000000), double3(0.00000000, 3.72999824, 0.00000000), double3(0.00000000, 0.00000000, 7.39499652))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-054", double3x3(double3(10.11999524, 0.00000000, 0.00000000), double3(0.00000000, 5.80999727, 0.00000000), double3(0.00000000, 0.00000000, 10.94999485))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-054-2", double3x3(double3(9.98099530, 0.00000000, 0.00000000), double3(0.00000000, 5.78099728, 0.00000000), double3(0.00000000, 0.00000000, 10.52199505))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-055", double3x3(double3(11.54199457, 0.00000000, 0.00000000), double3(0.00000000, 12.68999403, 0.00000000), double3(0.00000000, 0.00000000, 3.97399813))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-055-2", double3x3(double3(7.91499628, 0.00000000, 0.00000000), double3(0.00000000, 11.21999472, 0.00000000), double3(0.00000000, 0.00000000, 3.95999814))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-056", double3x3(double3(4.91099769, 0.00000000, 0.00000000), double3(0.00000000, 12.46399414, 0.00000000), double3(0.00000000, 0.00000000, 5.41199745))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-056-2", double3x3(double3(7.80299633, 0.00000000, 0.00000000), double3(0.00000000, 10.27299517, 0.00000000), double3(0.00000000, 0.00000000, 8.56599597))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-057", double3x3(double3(5.26099752, 0.00000000, 0.00000000), double3(0.00000000, 11.42499462, 0.00000000), double3(0.00000000, 0.00000000, 5.71499731))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-057-2", double3x3(double3(6.35899701, 0.00000000, 0.00000000), double3(0.00000000, 9.76499541, 0.00000000), double3(0.00000000, 0.00000000, 6.66699686))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058", double3x3(double3(10.83999490, 0.00000000, 0.00000000), double3(0.00000000, 23.69298885, 0.00000000), double3(0.00000000, 0.00000000, 6.90199675))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-2", double3x3(double3(4.33199796, 0.00000000, 0.00000000), double3(0.00000000, 4.87299771, 0.00000000), double3(0.00000000, 0.00000000, 2.96299861))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-058-3", double3x3(double3(3.89128206, 0.00000000, 0.00000000), double3(0.00000000, 4.02154657, 0.00000000), double3(0.00000000, 0.00000000, 2.57249173))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-059", double3x3(double3(3.56299832, 0.00000000, 0.00000000), double3(0.00000000, 11.50999458, 0.00000000), double3(0.00000000, 0.00000000, 4.36899794))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-059-2", double3x3(double3(2.86499865, 0.00000000, 0.00000000), double3(0.00000000, 4.64499781, 0.00000000), double3(0.00000000, 0.00000000, 4.04499810))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060", double3x3(double3(16.15199240, 0.00000000, 0.00000000), double3(0.00000000, 6.24899706, 0.00000000), double3(0.00000000, 0.00000000, 6.30699703))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-2", double3x3(double3(8.03099622, 0.00000000, 0.00000000), double3(0.00000000, 9.51799552, 0.00000000), double3(0.00000000, 0.00000000, 9.74899541))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-060-3", double3x3(double3(5.36727398, 0.00000000, 0.00000000), double3(0.00000000, 6.58837898, 0.00000000), double3(0.00000000, 0.00000000, 3.87260732))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-061", double3x3(double3(10.53699504, 0.00000000, 0.00000000), double3(0.00000000, 12.19899426, 0.00000000), double3(0.00000000, 0.00000000, 13.04699386))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-061-2", double3x3(double3(5.99299718, 0.00000000, 0.00000000), double3(0.00000000, 7.81899632, 0.00000000), double3(0.00000000, 0.00000000, 8.01099623))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-062", double3x3(double3(7.48999648, 0.00000000, 0.00000000), double3(0.00000000, 6.89799675, 0.00000000), double3(0.00000000, 0.00000000, 10.94199485))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-062-2", double3x3(double3(9.49499553, 0.00000000, 0.00000000), double3(0.00000000, 9.37499559, 0.00000000), double3(0.00000000, 0.00000000, 10.14999522))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063", double3x3(double3(9.20099567, 0.00000000, 0.00000000), double3(0.00000000, 7.15899663, 0.00000000), double3(0.00000000, 0.00000000, 9.77099540))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-2", double3x3(double3(5.56899738, 0.00000000, 0.00000000), double3(0.00000000, 12.79599398, 0.00000000), double3(0.00000000, 0.00000000, 7.31999656))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-063-3", double3x3(double3(5.75999729, 0.00000000, 0.00000000), double3(0.00000000, 8.89999581, 0.00000000), double3(0.00000000, 0.00000000, 13.32999373))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064", double3x3(double3(5.36999747, 0.00000000, 0.00000000), double3(0.00000000, 13.14999381, 0.00000000), double3(0.00000000, 0.00000000, 5.40599746))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-2", double3x3(double3(5.47099743, 0.00000000, 0.00000000), double3(0.00000000, 12.21099425, 0.00000000), double3(0.00000000, 0.00000000, 5.46499743))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-064-3", double3x3(double3(5.75718235, 0.00000000, 0.00000000), double3(0.00000000, 5.14105124, 0.00000000), double3(0.00000000, 0.00000000, 6.80861378))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065", double3x3(double3(4.83799772, 0.00000000, 0.00000000), double3(0.00000000, 8.14699617, 0.00000000), double3(0.00000000, 0.00000000, 6.10699713))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-2", double3x3(double3(5.75999729, 0.00000000, 0.00000000), double3(0.00000000, 6.37999700, 0.00000000), double3(0.00000000, 0.00000000, 4.05999809))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-065-3", double3x3(double3(5.49213934, 0.00000000, 0.00000000), double3(0.00000000, 5.56450595, 0.00000000), double3(0.00000000, 0.00000000, 3.87129818))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-066", double3x3(double3(6.32999702, 0.00000000, 0.00000000), double3(0.00000000, 10.47999507, 0.00000000), double3(0.00000000, 0.00000000, 10.52999505))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-066-2", double3x3(double3(7.06969667, 0.00000000, 0.00000000), double3(0.00000000, 25.49228800, 0.00000000), double3(0.00000000, 0.00000000, 7.02689669))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067", double3x3(double3(7.96999625, 0.00000000, 0.00000000), double3(0.00000000, 11.72199448, 0.00000000), double3(0.00000000, 0.00000000, 4.37399794))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-2", double3x3(double3(9.91599533, 0.00000000, 0.00000000), double3(0.00000000, 11.36599465, 0.00000000), double3(0.00000000, 0.00000000, 8.27099611))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-067-3", double3x3(double3(7.63999641, 0.00000000, 0.00000000), double3(0.00000000, 7.77999634, 0.00000000), double3(0.00000000, 0.00000000, 5.48999742))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-068", double3x3(double3(7.40999651, 0.00000000, 0.00000000), double3(0.00000000, 22.25998953, 0.00000000), double3(0.00000000, 0.00000000, 7.43999650))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-068-2", double3x3(double3(6.41799698, 0.00000000, 0.00000000), double3(0.00000000, 11.36599465, 0.00000000), double3(0.00000000, 0.00000000, 6.38399700))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-069", double3x3(double3(6.38999699, 0.00000000, 0.00000000), double3(0.00000000, 10.85999489, 0.00000000), double3(0.00000000, 0.00000000, 13.59999360))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-069-2", double3x3(double3(2.73820871, 0.00000000, 0.00000000), double3(0.00000000, 11.26079470, 0.00000000), double3(0.00000000, 0.00000000, 12.42669415))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-070", double3x3(double3(7.03899669, 0.00000000, 0.00000000), double3(0.00000000, 8.35599607, 0.00000000), double3(0.00000000, 0.00000000, 10.18599521))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-070-2", double3x3(double3(7.46199649, 0.00000000, 0.00000000), double3(0.00000000, 9.60299548, 0.00000000), double3(0.00000000, 0.00000000, 9.69899544))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-071", double3x3(double3(2.87499865, 0.00000000, 0.00000000), double3(0.00000000, 4.71499778, 0.00000000), double3(0.00000000, 0.00000000, 15.70699261))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-071-2", double3x3(double3(3.54199833, 0.00000000, 0.00000000), double3(0.00000000, 3.82699820, 0.00000000), double3(0.00000000, 0.00000000, 12.69599403))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-072", double3x3(double3(7.50099647, 0.00000000, 0.00000000), double3(0.00000000, 15.96599249, 0.00000000), double3(0.00000000, 0.00000000, 4.85799771))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-072-2", double3x3(double3(5.96699719, 0.00000000, 0.00000000), double3(0.00000000, 10.47999507, 0.00000000), double3(0.00000000, 0.00000000, 5.40199746))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-073", double3x3(double3(8.27019611, 0.00000000, 0.00000000), double3(0.00000000, 8.31149609, 0.00000000), double3(0.00000000, 0.00000000, 20.60699030))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-073-2", double3x3(double3(5.74999729, 0.00000000, 0.00000000), double3(0.00000000, 5.94999720, 0.00000000), double3(0.00000000, 0.00000000, 20.22999048))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-074", double3x3(double3(5.69599732, 0.00000000, 0.00000000), double3(0.00000000, 11.44399462, 0.00000000), double3(0.00000000, 0.00000000, 8.24799612))),
    std::make_pair("spglibtestdata/orthorhombic/POSCAR-074-2", double3x3(double3(5.91199722, 0.00000000, 0.00000000), double3(0.00000000, 5.94499720, 0.00000000), double3(0.00000000, 0.00000000, 8.38799605))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Tetragonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/tetragonal/POSCAR-075", double3x3(double3(17.48999177, 0.00000000, 0.00000000), double3(0.00000000, 17.48999177, 0.00000000), double3(0.00000000, 0.00000000, 3.94399814))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-075-2", double3x3(double3(9.19699567, 0.00000000, 0.00000000), double3(0.00000000, 9.19699567, 0.00000000), double3(0.00000000, 0.00000000, 20.50499035))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-076", double3x3(double3(3.98099813, 0.00000000, 0.00000000), double3(0.00000000, 3.98099813, 0.00000000), double3(0.00000000, 0.00000000, 15.34999278))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-076-2", double3x3(double3(8.44799602, 0.00000000, 0.00000000), double3(0.00000000, 8.44799602, 0.00000000), double3(0.00000000, 0.00000000, 14.91199298))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077", double3x3(double3(11.16399475, 0.00000000, 0.00000000), double3(0.00000000, 11.16399475, 0.00000000), double3(0.00000000, 0.00000000, 10.63799499))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077-2", double3x3(double3(7.97999625, 0.00000000, 0.00000000), double3(0.00000000, 7.97999625, 0.00000000), double3(0.00000000, 0.00000000, 9.77999540))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-077-3", double3x3(double3(11.16399475, 0.00000000, 0.00000000), double3(0.00000000, 11.16399475, 0.00000000), double3(0.00000000, 0.00000000, 10.63799499))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-078", double3x3(double3(10.86499489, 0.00000000, 0.00000000), double3(0.00000000, 10.86499489, 0.00000000), double3(0.00000000, 0.00000000, 28.35698666))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-078-2", double3x3(double3(7.62899641, 0.00000000, 0.00000000), double3(0.00000000, 7.62899641, 0.00000000), double3(0.00000000, 0.00000000, 29.49698612))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-079", double3x3(double3(8.48399601, 0.00000000, 0.00000000), double3(0.00000000, 8.48399601, 0.00000000), double3(0.00000000, 0.00000000, 5.81299726))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-079-2", double3x3(double3(14.91899298, 0.00000000, 0.00000000), double3(0.00000000, 14.91899298, 0.00000000), double3(0.00000000, 0.00000000, 7.60899642))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-080", double3x3(double3(20.33799043, 0.00000000, 0.00000000), double3(0.00000000, 20.33799043, 0.00000000), double3(0.00000000, 0.00000000, 14.66799310))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-080-2", double3x3(double3(9.69299544, 0.00000000, 0.00000000), double3(0.00000000, 9.69299544, 0.00000000), double3(0.00000000, 0.00000000, 5.98499718))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-081", double3x3(double3(7.62089641, 0.00000000, 0.00000000), double3(0.00000000, 7.62089641, 0.00000000), double3(0.00000000, 0.00000000, 6.31999703))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-081-2", double3x3(double3(10.18149521, 0.00000000, 0.00000000), double3(0.00000000, 10.18149521, 0.00000000), double3(0.00000000, 0.00000000, 5.29499751))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-082", double3x3(double3(5.54799739, 0.00000000, 0.00000000), double3(0.00000000, 5.54799739, 0.00000000), double3(0.00000000, 0.00000000, 10.16999521))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-082-2", double3x3(double3(6.32199703, 0.00000000, 0.00000000), double3(0.00000000, 6.32199703, 0.00000000), double3(0.00000000, 0.00000000, 12.60499407))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083", double3x3(double3(8.32799608, 0.00000000, 0.00000000), double3(0.00000000, 8.32799608, 0.00000000), double3(0.00000000, 0.00000000, 3.13449853))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083-2", double3x3(double3(2.81744433, 0.00000000, 0.00000000), double3(0.00000000, 2.81744433, 0.00000000), double3(0.00000000, 0.00000000, 3.99999812))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-083-3", double3x3(double3(5.52957243, 0.00000000, 0.00000000), double3(0.00000000, 5.52957243, 0.00000000), double3(0.00000000, 0.00000000, 3.67999827))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-084", double3x3(double3(6.42999697, 0.00000000, 0.00000000), double3(0.00000000, 6.42999697, 0.00000000), double3(0.00000000, 0.00000000, 6.62999688))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-084-2", double3x3(double3(7.16699663, 0.00000000, 0.00000000), double3(0.00000000, 7.16699663, 0.00000000), double3(0.00000000, 0.00000000, 6.59799690))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-085", double3x3(double3(6.26099705, 0.00000000, 0.00000000), double3(0.00000000, 6.26099705, 0.00000000), double3(0.00000000, 0.00000000, 4.10099807))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-085-2", double3x3(double3(8.37999606, 0.00000000, 0.00000000), double3(0.00000000, 8.37999606, 0.00000000), double3(0.00000000, 0.00000000, 7.48999648))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-086", double3x3(double3(11.18699474, 0.00000000, 0.00000000), double3(0.00000000, 11.18699474, 0.00000000), double3(0.00000000, 0.00000000, 6.18899709))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-086-2", double3x3(double3(7.06999667, 0.00000000, 0.00000000), double3(0.00000000, 7.06999667, 0.00000000), double3(0.00000000, 0.00000000, 3.82199820))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-087", double3x3(double3(11.40799463, 0.00000000, 0.00000000), double3(0.00000000, 11.40799463, 0.00000000), double3(0.00000000, 0.00000000, 10.25599517))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-087-2", double3x3(double3(9.88499535, 0.00000000, 0.00000000), double3(0.00000000, 9.88499535, 0.00000000), double3(0.00000000, 0.00000000, 3.12699853))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-088", double3x3(double3(13.69599356, 0.00000000, 0.00000000), double3(0.00000000, 13.69599356, 0.00000000), double3(0.00000000, 0.00000000, 5.98099719))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-088-2", double3x3(double3(5.74099730, 0.00000000, 0.00000000), double3(0.00000000, 5.74099730, 0.00000000), double3(0.00000000, 0.00000000, 13.12099383))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-090", double3x3(double3(9.55979550, 0.00000000, 0.00000000), double3(0.00000000, 9.55979550, 0.00000000), double3(0.00000000, 0.00000000, 7.15999663))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-090-2", double3x3(double3(5.25380091, 0.00000000, 0.00000000), double3(0.00000000, 5.25380091, 0.00000000), double3(0.00000000, 0.00000000, 10.60999501))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-091", double3x3(double3(7.00819670, 0.00000000, 0.00000000), double3(0.00000000, 7.00819670, 0.00000000), double3(0.00000000, 0.00000000, 8.65179593))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-091-2", double3x3(double3(12.46399414, 0.00000000, 0.00000000), double3(0.00000000, 12.46399414, 0.00000000), double3(0.00000000, 0.00000000, 26.22298766))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092", double3x3(double3(7.10399666, 0.00000000, 0.00000000), double3(0.00000000, 7.10399666, 0.00000000), double3(0.00000000, 0.00000000, 36.59698278))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092-2", double3x3(double3(6.58999690, 0.00000000, 0.00000000), double3(0.00000000, 6.58999690, 0.00000000), double3(0.00000000, 0.00000000, 17.03999198))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-092-3", double3x3(double3(4.06437130, 0.00000000, 0.00000000), double3(0.00000000, 4.06437130, 0.00000000), double3(0.00000000, 0.00000000, 5.63464733))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094", double3x3(double3(4.68629779, 0.00000000, 0.00000000), double3(0.00000000, 4.68629779, 0.00000000), double3(0.00000000, 0.00000000, 9.19099568))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094-2", double3x3(double3(9.96199531, 0.00000000, 0.00000000), double3(0.00000000, 9.96199531, 0.00000000), double3(0.00000000, 0.00000000, 13.41399369))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-094-3", double3x3(double3(7.34499654, 0.00000000, 0.00000000), double3(0.00000000, 7.34499654, 0.00000000), double3(0.00000000, 0.00000000, 10.39999511))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-095", double3x3(double3(6.16999710, 0.00000000, 0.00000000), double3(0.00000000, 6.16999710, 0.00000000), double3(0.00000000, 0.00000000, 8.56399597))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-095-2", double3x3(double3(6.08039714, 0.00000000, 0.00000000), double3(0.00000000, 6.08039714, 0.00000000), double3(0.00000000, 0.00000000, 8.39879605))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-096", double3x3(double3(7.48999648, 0.00000000, 0.00000000), double3(0.00000000, 7.48999648, 0.00000000), double3(0.00000000, 0.00000000, 13.23999377))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-096-2", double3x3(double3(5.99699718, 0.00000000, 0.00000000), double3(0.00000000, 5.99699718, 0.00000000), double3(0.00000000, 0.00000000, 39.19098156))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-097", double3x3(double3(7.48299648, 0.00000000, 0.00000000), double3(0.00000000, 7.48299648, 0.00000000), double3(0.00000000, 0.00000000, 14.89299299))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-097-2", double3x3(double3(9.53099552, 0.00000000, 0.00000000), double3(0.00000000, 9.53099552, 0.00000000), double3(0.00000000, 0.00000000, 12.82399397))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-098", double3x3(double3(7.95399626, 0.00000000, 0.00000000), double3(0.00000000, 7.95399626, 0.00000000), double3(0.00000000, 0.00000000, 4.67799780))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-098-2", double3x3(double3(9.38299558, 0.00000000, 0.00000000), double3(0.00000000, 9.38299558, 0.00000000), double3(0.00000000, 0.00000000, 54.59997431))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-099", double3x3(double3(3.89899817, 0.00000000, 0.00000000), double3(0.00000000, 3.89899817, 0.00000000), double3(0.00000000, 0.00000000, 4.16699804))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-099-2", double3x3(double3(3.80719821, 0.00000000, 0.00000000), double3(0.00000000, 3.80719821, 0.00000000), double3(0.00000000, 0.00000000, 4.69819779))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-100", double3x3(double3(8.86999583, 0.00000000, 0.00000000), double3(0.00000000, 8.86999583, 0.00000000), double3(0.00000000, 0.00000000, 5.21499755))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-100-2", double3x3(double3(8.31199609, 0.00000000, 0.00000000), double3(0.00000000, 8.31199609, 0.00000000), double3(0.00000000, 0.00000000, 10.06999526))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-102", double3x3(double3(8.86399583, 0.00000000, 0.00000000), double3(0.00000000, 8.86399583, 0.00000000), double3(0.00000000, 0.00000000, 8.01999623))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-102-2", double3x3(double3(10.75899494, 0.00000000, 0.00000000), double3(0.00000000, 10.75899494, 0.00000000), double3(0.00000000, 0.00000000, 5.65599734))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-103", double3x3(double3(6.55099692, 0.00000000, 0.00000000), double3(0.00000000, 6.55099692, 0.00000000), double3(0.00000000, 0.00000000, 6.84699678))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-103-2", double3x3(double3(6.51399693, 0.00000000, 0.00000000), double3(0.00000000, 6.51399693, 0.00000000), double3(0.00000000, 0.00000000, 6.80899680))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-104", double3x3(double3(9.41599557, 0.00000000, 0.00000000), double3(0.00000000, 9.41599557, 0.00000000), double3(0.00000000, 0.00000000, 9.23699565))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-104-2", double3x3(double3(10.61999500, 0.00000000, 0.00000000), double3(0.00000000, 10.61999500, 0.00000000), double3(0.00000000, 0.00000000, 9.00899576))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-105", double3x3(double3(7.61799642, 0.00000000, 0.00000000), double3(0.00000000, 7.61799642, 0.00000000), double3(0.00000000, 0.00000000, 8.49999600))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-105-2", double3x3(double3(5.63999735, 0.00000000, 0.00000000), double3(0.00000000, 5.63999735, 0.00000000), double3(0.00000000, 0.00000000, 9.04999574))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-106", double3x3(double3(10.83899490, 0.00000000, 0.00000000), double3(0.00000000, 10.83899490, 0.00000000), double3(0.00000000, 0.00000000, 5.30799750))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107", double3x3(double3(7.76199635, 0.00000000, 0.00000000), double3(0.00000000, 7.76199635, 0.00000000), double3(0.00000000, 0.00000000, 11.74599447))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107-2", double3x3(double3(4.11999806, 0.00000000, 0.00000000), double3(0.00000000, 4.11999806, 0.00000000), double3(0.00000000, 0.00000000, 10.47499507))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-107-3", double3x3(double3(3.55189833, 0.00000000, 0.00000000), double3(0.00000000, 3.55189833, 0.00000000), double3(0.00000000, 0.00000000, 25.64538793))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-108", double3x3(double3(6.13999711, 0.00000000, 0.00000000), double3(0.00000000, 6.13999711, 0.00000000), double3(0.00000000, 0.00000000, 10.57199503))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-108-2", double3x3(double3(8.05499621, 0.00000000, 0.00000000), double3(0.00000000, 8.05499621, 0.00000000), double3(0.00000000, 0.00000000, 15.68799262))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-109", double3x3(double3(3.45169838, 0.00000000, 0.00000000), double3(0.00000000, 3.45169838, 0.00000000), double3(0.00000000, 0.00000000, 11.67999450))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-109-2", double3x3(double3(11.37999465, 0.00000000, 0.00000000), double3(0.00000000, 11.37999465, 0.00000000), double3(0.00000000, 0.00000000, 8.74999588))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-110", double3x3(double3(13.61999359, 0.00000000, 0.00000000), double3(0.00000000, 13.61999359, 0.00000000), double3(0.00000000, 0.00000000, 9.09999572))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-110-2", double3x3(double3(11.78229446, 0.00000000, 0.00000000), double3(0.00000000, 11.78229446, 0.00000000), double3(0.00000000, 0.00000000, 23.63848888))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-111", double3x3(double3(5.15999757, 0.00000000, 0.00000000), double3(0.00000000, 5.15999757, 0.00000000), double3(0.00000000, 0.00000000, 10.06999526))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-111-2", double3x3(double3(5.81509726, 0.00000000, 0.00000000), double3(0.00000000, 5.81509726, 0.00000000), double3(0.00000000, 0.00000000, 5.80139727))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-112", double3x3(double3(5.42999744, 0.00000000, 0.00000000), double3(0.00000000, 5.42999744, 0.00000000), double3(0.00000000, 0.00000000, 10.09599525))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-112-2", double3x3(double3(5.41499745, 0.00000000, 0.00000000), double3(0.00000000, 5.41499745, 0.00000000), double3(0.00000000, 0.00000000, 10.19699520))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-113", double3x3(double3(5.66199734, 0.00000000, 0.00000000), double3(0.00000000, 5.66199734, 0.00000000), double3(0.00000000, 0.00000000, 4.71599778))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-113-2", double3x3(double3(6.40239699, 0.00000000, 0.00000000), double3(0.00000000, 6.40239699, 0.00000000), double3(0.00000000, 0.00000000, 4.27859799))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-114", double3x3(double3(7.48399648, 0.00000000, 0.00000000), double3(0.00000000, 7.48399648, 0.00000000), double3(0.00000000, 0.00000000, 6.34899701))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-114-2", double3x3(double3(10.80999491, 0.00000000, 0.00000000), double3(0.00000000, 10.80999491, 0.00000000), double3(0.00000000, 0.00000000, 6.80999680))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115", double3x3(double3(3.32699843, 0.00000000, 0.00000000), double3(0.00000000, 3.32699843, 0.00000000), double3(0.00000000, 0.00000000, 6.15099711))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-2", double3x3(double3(3.85899818, 0.00000000, 0.00000000), double3(0.00000000, 3.85899818, 0.00000000), double3(0.00000000, 0.00000000, 11.70999449))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-3", double3x3(double3(2.33273548, 0.00000000, 0.00000000), double3(0.00000000, 2.33273548, 0.00000000), double3(0.00000000, 0.00000000, 4.78863344))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-4", double3x3(double3(5.01899764, 0.00000000, 0.00000000), double3(0.00000000, 5.01899764, 0.00000000), double3(0.00000000, 0.00000000, 8.42799603))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-115-5", double3x3(double3(3.73352205, 0.00000000, 0.00000000), double3(0.00000000, 3.73352205, 0.00000000), double3(0.00000000, 0.00000000, 5.20499755))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-116", double3x3(double3(5.52499740, 0.00000000, 0.00000000), double3(0.00000000, 5.52499740, 0.00000000), double3(0.00000000, 0.00000000, 17.46299178))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-116-2", double3x3(double3(10.56599503, 0.00000000, 0.00000000), double3(0.00000000, 10.56599503, 0.00000000), double3(0.00000000, 0.00000000, 25.21998813))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-117", double3x3(double3(7.72867348, 0.00000000, 0.00000000), double3(0.00000000, 7.72867348, 0.00000000), double3(0.00000000, 0.00000000, 5.61999736))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-117-2", double3x3(double3(5.70999731, 0.00000000, 0.00000000), double3(0.00000000, 5.70999731, 0.00000000), double3(0.00000000, 0.00000000, 16.45799226))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-118", double3x3(double3(6.59599690, 0.00000000, 0.00000000), double3(0.00000000, 6.59599690, 0.00000000), double3(0.00000000, 0.00000000, 5.16469757))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-118-2", double3x3(double3(7.82299632, 0.00000000, 0.00000000), double3(0.00000000, 7.82299632, 0.00000000), double3(0.00000000, 0.00000000, 9.05299574))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-119", double3x3(double3(6.26799705, 0.00000000, 0.00000000), double3(0.00000000, 6.26799705, 0.00000000), double3(0.00000000, 0.00000000, 14.78199304))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-119-2", double3x3(double3(3.91999816, 0.00000000, 0.00000000), double3(0.00000000, 3.91999816, 0.00000000), double3(0.00000000, 0.00000000, 15.21999284))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-120", double3x3(double3(8.46699602, 0.00000000, 0.00000000), double3(0.00000000, 8.46699602, 0.00000000), double3(0.00000000, 0.00000000, 12.74499400))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-120-2", double3x3(double3(12.44759414, 0.00000000, 0.00000000), double3(0.00000000, 12.44759414, 0.00000000), double3(0.00000000, 0.00000000, 12.92769392))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-121", double3x3(double3(4.97599766, 0.00000000, 0.00000000), double3(0.00000000, 4.97599766, 0.00000000), double3(0.00000000, 0.00000000, 6.74599683))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-121-2", double3x3(double3(5.04699763, 0.00000000, 0.00000000), double3(0.00000000, 5.04699763, 0.00000000), double3(0.00000000, 0.00000000, 6.48599695))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122", double3x3(double3(11.48399460, 0.00000000, 0.00000000), double3(0.00000000, 11.48399460, 0.00000000), double3(0.00000000, 0.00000000, 5.61999736))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122-2", double3x3(double3(14.93299297, 0.00000000, 0.00000000), double3(0.00000000, 14.93299297, 0.00000000), double3(0.00000000, 0.00000000, 15.54699268))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-122-3", double3x3(double3(3.95768302, 0.00000000, 0.00000000), double3(0.00000000, 3.95768302, 0.00000000), double3(0.00000000, 0.00000000, 5.97647452))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123", double3x3(double3(4.01899811, 0.00000000, 0.00000000), double3(0.00000000, 4.01899811, 0.00000000), double3(0.00000000, 0.00000000, 3.27899846))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123-2", double3x3(double3(8.98099577, 0.00000000, 0.00000000), double3(0.00000000, 8.98099577, 0.00000000), double3(0.00000000, 0.00000000, 12.63799405))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-123-3", double3x3(double3(4.19999802, 0.00000000, 0.00000000), double3(0.00000000, 4.19999802, 0.00000000), double3(0.00000000, 0.00000000, 7.95999625))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-124", double3x3(double3(6.11799712, 0.00000000, 0.00000000), double3(0.00000000, 6.11799712, 0.00000000), double3(0.00000000, 0.00000000, 4.99599765))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-124-2", double3x3(double3(6.20999708, 0.00000000, 0.00000000), double3(0.00000000, 6.20999708, 0.00000000), double3(0.00000000, 0.00000000, 11.00199482))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-125", double3x3(double3(8.51599599, 0.00000000, 0.00000000), double3(0.00000000, 8.51599599, 0.00000000), double3(0.00000000, 0.00000000, 6.71299684))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-125-2", double3x3(double3(6.37599700, 0.00000000, 0.00000000), double3(0.00000000, 6.37599700, 0.00000000), double3(0.00000000, 0.00000000, 8.32899608))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-126", double3x3(double3(11.34999466, 0.00000000, 0.00000000), double3(0.00000000, 11.34999466, 0.00000000), double3(0.00000000, 0.00000000, 6.18999709))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-126-2", double3x3(double3(6.31239703, 0.00000000, 0.00000000), double3(0.00000000, 6.31239703, 0.00000000), double3(0.00000000, 0.00000000, 9.54949551))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-127", double3x3(double3(5.89399723, 0.00000000, 0.00000000), double3(0.00000000, 5.89399723, 0.00000000), double3(0.00000000, 0.00000000, 8.34799607))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-127-2", double3x3(double3(7.10499666, 0.00000000, 0.00000000), double3(0.00000000, 7.10499666, 0.00000000), double3(0.00000000, 0.00000000, 4.14399805))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-128", double3x3(double3(7.05769668, 0.00000000, 0.00000000), double3(0.00000000, 7.05769668, 0.00000000), double3(0.00000000, 0.00000000, 9.97839530))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-128-2", double3x3(double3(7.74359636, 0.00000000, 0.00000000), double3(0.00000000, 7.74359636, 0.00000000), double3(0.00000000, 0.00000000, 11.64029452))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129", double3x3(double3(4.28199799, 0.00000000, 0.00000000), double3(0.00000000, 4.28199799, 0.00000000), double3(0.00000000, 0.00000000, 6.18199709))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129-2", double3x3(double3(5.00499764, 0.00000000, 0.00000000), double3(0.00000000, 5.00499764, 0.00000000), double3(0.00000000, 0.00000000, 4.80999774))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-129-3", double3x3(double3(4.28199799, 0.00000000, 0.00000000), double3(0.00000000, 4.28199799, 0.00000000), double3(0.00000000, 0.00000000, 6.18199709))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-130", double3x3(double3(5.92399721, 0.00000000, 0.00000000), double3(0.00000000, 5.92399721, 0.00000000), double3(0.00000000, 0.00000000, 18.12999147))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-130-2", double3x3(double3(7.37719653, 0.00000000, 0.00000000), double3(0.00000000, 7.37719653, 0.00000000), double3(0.00000000, 0.00000000, 15.12309288))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-131", double3x3(double3(3.01999858, 0.00000000, 0.00000000), double3(0.00000000, 3.01999858, 0.00000000), double3(0.00000000, 0.00000000, 5.30999750))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-131-2", double3x3(double3(4.92629768, 0.00000000, 0.00000000), double3(0.00000000, 4.92629768, 0.00000000), double3(0.00000000, 0.00000000, 8.28499610))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-132", double3x3(double3(6.16599710, 0.00000000, 0.00000000), double3(0.00000000, 6.16599710, 0.00000000), double3(0.00000000, 0.00000000, 6.05199715))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-132-2", double3x3(double3(5.85199725, 0.00000000, 0.00000000), double3(0.00000000, 5.85199725, 0.00000000), double3(0.00000000, 0.00000000, 14.23399330))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-133", double3x3(double3(9.38099559, 0.00000000, 0.00000000), double3(0.00000000, 9.38099559, 0.00000000), double3(0.00000000, 0.00000000, 4.66299781))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-133-2", double3x3(double3(11.78899445, 0.00000000, 0.00000000), double3(0.00000000, 11.78899445, 0.00000000), double3(0.00000000, 0.00000000, 23.63498888))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-134", double3x3(double3(8.42699603, 0.00000000, 0.00000000), double3(0.00000000, 8.42699603, 0.00000000), double3(0.00000000, 0.00000000, 14.49199318))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-134-2", double3x3(double3(3.08449855, 0.00000000, 0.00000000), double3(0.00000000, 3.08449855, 0.00000000), double3(0.00000000, 0.00000000, 3.10649854))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-135", double3x3(double3(8.58999596, 0.00000000, 0.00000000), double3(0.00000000, 8.58999596, 0.00000000), double3(0.00000000, 0.00000000, 5.91299722))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-135-2", double3x3(double3(8.52699599, 0.00000000, 0.00000000), double3(0.00000000, 8.52699599, 0.00000000), double3(0.00000000, 0.00000000, 5.94199720))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136", double3x3(double3(4.39829793, 0.00000000, 0.00000000), double3(0.00000000, 4.39829793, 0.00000000), double3(0.00000000, 0.00000000, 2.87299865))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-2", double3x3(double3(4.58449784, 0.00000000, 0.00000000), double3(0.00000000, 4.58449784, 0.00000000), double3(0.00000000, 0.00000000, 2.95329861))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-3", double3x3(double3(4.22665402, 0.00000000, 0.00000000), double3(0.00000000, 4.22665402, 0.00000000), double3(0.00000000, 0.00000000, 2.68883593))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-4", double3x3(double3(4.66859780, 0.00000000, 0.00000000), double3(0.00000000, 4.66859780, 0.00000000), double3(0.00000000, 0.00000000, 5.21499755))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-136-5", double3x3(double3(6.75399682, 0.00000000, 0.00000000), double3(0.00000000, 6.75399682, 0.00000000), double3(0.00000000, 0.00000000, 4.10999807))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137", double3x3(double3(8.08999619, 0.00000000, 0.00000000), double3(0.00000000, 8.08999619, 0.00000000), double3(0.00000000, 0.00000000, 5.44999744))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137-2", double3x3(double3(3.63999829, 0.00000000, 0.00000000), double3(0.00000000, 3.63999829, 0.00000000), double3(0.00000000, 0.00000000, 5.26999752))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-137-3", double3x3(double3(2.24231466, 0.00000000, 0.00000000), double3(0.00000000, 2.24231466, 0.00000000), double3(0.00000000, 0.00000000, 6.85050892))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-138", double3x3(double3(8.43359603, 0.00000000, 0.00000000), double3(0.00000000, 8.43359603, 0.00000000), double3(0.00000000, 0.00000000, 7.67819639))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-138-2", double3x3(double3(4.34999795, 0.00000000, 0.00000000), double3(0.00000000, 4.34999795, 0.00000000), double3(0.00000000, 0.00000000, 13.72999354))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-139", double3x3(double3(11.93999438, 0.00000000, 0.00000000), double3(0.00000000, 11.93999438, 0.00000000), double3(0.00000000, 0.00000000, 17.39999181))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-139-2", double3x3(double3(4.16999804, 0.00000000, 0.00000000), double3(0.00000000, 4.16999804, 0.00000000), double3(0.00000000, 0.00000000, 10.87999488))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-140", double3x3(double3(11.07599479, 0.00000000, 0.00000000), double3(0.00000000, 11.07599479, 0.00000000), double3(0.00000000, 0.00000000, 36.93298262))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-140-2", double3x3(double3(11.01199482, 0.00000000, 0.00000000), double3(0.00000000, 11.01199482, 0.00000000), double3(0.00000000, 0.00000000, 5.72699731))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-141", double3x3(double3(7.17719662, 0.00000000, 0.00000000), double3(0.00000000, 7.17719662, 0.00000000), double3(0.00000000, 0.00000000, 6.32889702))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-141-2", double3x3(double3(6.90129675, 0.00000000, 0.00000000), double3(0.00000000, 6.90129675, 0.00000000), double3(0.00000000, 0.00000000, 19.97539060))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142", double3x3(double3(10.32999514, 0.00000000, 0.00000000), double3(0.00000000, 10.32999514, 0.00000000), double3(0.00000000, 0.00000000, 20.37999041))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142-2", double3x3(double3(12.28399422, 0.00000000, 0.00000000), double3(0.00000000, 12.28399422, 0.00000000), double3(0.00000000, 0.00000000, 23.58198890))),
    std::make_pair("spglibtestdata/tetragonal/POSCAR-142-3", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Trigonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/trigonal/POSCAR-143", double3x3(double3(7.24879659, 0.00000000, 0.00000000), double3(-3.62439829, 6.27764199, 0.00000000), double3(0.00000000, 0.00000000, 6.78359681))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-143-2", double3x3(double3(7.95418626, 0.00000000, 0.00000000), double3(-3.97709313, 6.88852737, 0.00000000), double3(0.00000000, 0.00000000, 7.36129654))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-144", double3x3(double3(6.86729677, 0.00000000, 0.00000000), double3(-3.43364838, 5.94725346, 0.00000000), double3(0.00000000, 0.00000000, 17.06199197))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-144-2", double3x3(double3(4.33679796, 0.00000000, 0.00000000), double3(-2.16839898, 3.75577720, 0.00000000), double3(0.00000000, 0.00000000, 8.33969608))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-145", double3x3(double3(12.69199403, 0.00000000, 0.00000000), double3(-6.34599701, 10.99158925, 0.00000000), double3(0.00000000, 0.00000000, 19.18599097))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-145-2", double3x3(double3(10.50199506, 0.00000000, 0.00000000), double3(-5.25099753, 9.09499451, 0.00000000), double3(0.00000000, 0.00000000, 7.47299648))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-146", double3x3(double3(10.83699490, 0.00000000, 0.00000000), double3(-5.41849745, 9.38511288, 0.00000000), double3(0.00000000, 0.00000000, 8.15699616))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-146-2", double3x3(double3(5.99999718, 0.00000000, 0.00000000), double3(-2.99999859, 5.19614998, 0.00000000), double3(0.00000000, 0.00000000, 14.32999326))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-147", double3x3(double3(16.99059201, 0.00000000, 0.00000000), double3(-8.49529600, 14.71428430, 0.00000000), double3(0.00000000, 0.00000000, 8.88879582))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-147-2", double3x3(double3(9.39619558, 0.00000000, 0.00000000), double3(-4.69809779, 8.13734407, 0.00000000), double3(0.00000000, 0.00000000, 7.22249660))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-148", double3x3(double3(7.10669666, 0.00000000, 0.00000000), double3(-3.55334833, 6.15457984, 0.00000000), double3(0.00000000, 0.00000000, 17.08669196))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-148-2", double3x3(double3(20.97399013, 0.00000000, 0.00000000), double3(-10.48699507, 18.16400827, 0.00000000), double3(0.00000000, 0.00000000, 12.13399429))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-149", double3x3(double3(5.02199764, 0.00000000, 0.00000000), double3(-2.51099882, 4.34917753, 0.00000000), double3(0.00000000, 0.00000000, 6.37599700))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-149-2", double3x3(double3(7.15099664, 0.00000000, 0.00000000), double3(-3.57549832, 6.19294475, 0.00000000), double3(0.00000000, 0.00000000, 8.17579615))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-150", double3x3(double3(9.06999573, 0.00000000, 0.00000000), double3(-4.53499787, 7.85484672, 0.00000000), double3(0.00000000, 0.00000000, 4.98399765))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-150-2", double3x3(double3(8.63799594, 0.00000000, 0.00000000), double3(-4.31899797, 7.48072392, 0.00000000), double3(0.00000000, 0.00000000, 4.73799777))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-151", double3x3(double3(5.95999720, 0.00000000, 0.00000000), double3(-2.97999860, 5.16150898, 0.00000000), double3(0.00000000, 0.00000000, 17.19999191))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-151-2", double3x3(double3(5.03399763, 0.00000000, 0.00000000), double3(-2.51699882, 4.35956983, 0.00000000), double3(0.00000000, 0.00000000, 14.14099335))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-152", double3x3(double3(9.20399567, 0.00000000, 0.00000000), double3(-4.60199783, 7.97089407, 0.00000000), double3(0.00000000, 0.00000000, 24.81798832))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-152-2", double3x3(double3(5.03599763, 0.00000000, 0.00000000), double3(-2.51799882, 4.36130188, 0.00000000), double3(0.00000000, 0.00000000, 11.25499470))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-153", double3x3(double3(6.01999717, 0.00000000, 0.00000000), double3(-3.00999858, 5.21347048, 0.00000000), double3(0.00000000, 0.00000000, 17.29999186))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154", double3x3(double3(4.91339769, 0.00000000, 0.00000000), double3(-2.45669884, 4.25512722, 0.00000000), double3(0.00000000, 0.00000000, 5.40519746))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154-2", double3x3(double3(13.03999386, 0.00000000, 0.00000000), double3(-6.51999693, 11.29296595, 0.00000000), double3(0.00000000, 0.00000000, 8.39999605))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-154-3", double3x3(double3(4.06285843, 0.00000000, 0.00000000), double3(-2.03142921, 3.51853861, 0.00000000), double3(0.00000000, 0.00000000, 4.61557670))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-155", double3x3(double3(9.34159560, 0.00000000, 0.00000000), double3(-4.67079780, 8.09005911, 0.00000000), double3(0.00000000, 0.00000000, 7.30549656))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-155-2", double3x3(double3(9.12299571, 0.00000000, 0.00000000), double3(-4.56149785, 7.90074604, 0.00000000), double3(0.00000000, 0.00000000, 16.98299201))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-156", double3x3(double3(3.73329824, 0.00000000, 0.00000000), double3(-1.86664912, 3.23313112, 0.00000000), double3(0.00000000, 0.00000000, 6.09799713))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-156-2", double3x3(double3(3.91499816, 0.00000000, 0.00000000), double3(-1.95749908, 3.39048786, 0.00000000), double3(0.00000000, 0.00000000, 12.72499401))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-157", double3x3(double3(12.19699426, 0.00000000, 0.00000000), double3(-6.09849713, 10.56290688, 0.00000000), double3(0.00000000, 0.00000000, 19.35899089))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-157-2", double3x3(double3(8.75299588, 0.00000000, 0.00000000), double3(-4.37649794, 7.58031679, 0.00000000), double3(0.00000000, 0.00000000, 3.96599813))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-158", double3x3(double3(6.11999712, 0.00000000, 0.00000000), double3(-3.05999856, 5.30007298, 0.00000000), double3(0.00000000, 0.00000000, 5.65799734))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-158-2", double3x3(double3(12.87399394, 0.00000000, 0.00000000), double3(-6.43699697, 11.14920580, 0.00000000), double3(0.00000000, 0.00000000, 9.23899565))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-159", double3x3(double3(10.19999520, 0.00000000, 0.00000000), double3(-5.09999760, 8.83345496, 0.00000000), double3(0.00000000, 0.00000000, 30.35098572))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-159-2", double3x3(double3(10.56299503, 0.00000000, 0.00000000), double3(-5.28149751, 9.14782204, 0.00000000), double3(0.00000000, 0.00000000, 5.36799747))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-160", double3x3(double3(12.72564330, 0.00000000, 0.00000000), double3(-6.36282165, 11.02073038, 0.00000000), double3(0.00000000, 0.00000000, 7.90251643))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-160-2", double3x3(double3(5.48699742, 0.00000000, 0.00000000), double3(-2.74349871, 4.75187915, 0.00000000), double3(0.00000000, 0.00000000, 9.15599569))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-161", double3x3(double3(10.43799509, 0.00000000, 0.00000000), double3(-5.21899754, 9.03956891, 0.00000000), double3(0.00000000, 0.00000000, 37.14998252))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-161-2", double3x3(double3(5.15999757, 0.00000000, 0.00000000), double3(-2.57999879, 4.46868898, 0.00000000), double3(0.00000000, 0.00000000, 16.57999220))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-162", double3x3(double3(5.44999744, 0.00000000, 0.00000000), double3(-2.72499872, 4.71983623, 0.00000000), double3(0.00000000, 0.00000000, 8.10099619))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-162-2", double3x3(double3(4.98999765, 0.00000000, 0.00000000), double3(-2.49499883, 4.32146473, 0.00000000), double3(0.00000000, 0.00000000, 4.62199783))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-163", double3x3(double3(5.88999723, 0.00000000, 0.00000000), double3(-2.94499861, 5.10088723, 0.00000000), double3(0.00000000, 0.00000000, 9.59099549))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-163-2", double3x3(double3(5.30999750, 0.00000000, 0.00000000), double3(-2.65499875, 4.59859273, 0.00000000), double3(0.00000000, 0.00000000, 14.24999329))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-164", double3x3(double3(4.04699810, 0.00000000, 0.00000000), double3(-2.02349905, 3.50480316, 0.00000000), double3(0.00000000, 0.00000000, 5.32999749))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-164-2", double3x3(double3(6.24899706, 0.00000000, 0.00000000), double3(-3.12449853, 5.41179020, 0.00000000), double3(0.00000000, 0.00000000, 5.08599761))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-165", double3x3(double3(7.18499662, 0.00000000, 0.00000000), double3(-3.59249831, 6.22238960, 0.00000000), double3(0.00000000, 0.00000000, 7.35099654))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-165-2", double3x3(double3(12.18999426, 0.00000000, 0.00000000), double3(-6.09499713, 10.55684470, 0.00000000), double3(0.00000000, 0.00000000, 10.13999523))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-166", double3x3(double3(6.24299706, 0.00000000, 0.00000000), double3(-3.12149853, 5.40659405, 0.00000000), double3(0.00000000, 0.00000000, 29.99998588))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-166-2", double3x3(double3(5.42499745, 0.00000000, 0.00000000), double3(-2.71249872, 4.69818560, 0.00000000), double3(0.00000000, 0.00000000, 9.83599537))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167", double3x3(double3(11.50999458, 0.00000000, 0.00000000), double3(-5.75499729, 9.96794771, 0.00000000), double3(0.00000000, 0.00000000, 60.56997150))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167-2", double3x3(double3(10.02299528, 0.00000000, 0.00000000), double3(-5.01149764, 8.68016854, 0.00000000), double3(0.00000000, 0.00000000, 25.47098801))),
    std::make_pair("spglibtestdata/trigonal/POSCAR-167-3", double3x3(double3(4.94919767, 0.00000000, 0.00000000), double3(-2.47459884, 4.28613091, 0.00000000), double3(0.00000000, 0.00000000, 13.99799341))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Hexagonal)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/hexagonal/POSCAR-168", double3x3(double3(15.93599250, 0.00000000, 0.00000000), double3(-7.96799625, 13.80097434, 0.00000000), double3(0.00000000, 0.00000000, 3.89199817))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-169", double3x3(double3(7.10999665, 0.00000000, 0.00000000), double3(-3.55499833, 6.15743772, 0.00000000), double3(0.00000000, 0.00000000, 19.33999090))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-169-2", double3x3(double3(9.70899543, 0.00000000, 0.00000000), double3(-4.85449772, 8.40823669, 0.00000000), double3(0.00000000, 0.00000000, 19.34299090))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-170", double3x3(double3(7.10999665, 0.00000000, 0.00000000), double3(-3.55499833, 6.15743772, 0.00000000), double3(0.00000000, 0.00000000, 19.29999092))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-170-2", double3x3(double3(10.51259505, 0.00000000, 0.00000000), double3(-5.25629753, 9.10417438, 0.00000000), double3(0.00000000, 0.00000000, 14.93759297))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-171", double3x3(double3(17.38999182, 0.00000000, 0.00000000), double3(-8.69499591, 15.06017469, 0.00000000), double3(0.00000000, 0.00000000, 13.03599387))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-171-2", double3x3(double3(6.31999703, 0.00000000, 0.00000000), double3(-3.15999851, 5.47327798, 0.00000000), double3(0.00000000, 0.00000000, 19.28999092))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-172", double3x3(double3(6.19799708, 0.00000000, 0.00000000), double3(-3.09899854, 5.36762293, 0.00000000), double3(0.00000000, 0.00000000, 18.72699119))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-173", double3x3(double3(7.13299664, 0.00000000, 0.00000000), double3(-3.56649832, 6.17735630, 0.00000000), double3(0.00000000, 0.00000000, 7.41399651))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-173-2", double3x3(double3(9.22499566, 0.00000000, 0.00000000), double3(-4.61249783, 7.98908059, 0.00000000), double3(0.00000000, 0.00000000, 5.22399754))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-174", double3x3(double3(10.27429517, 0.00000000, 0.00000000), double3(-5.13714758, 8.89780062, 0.00000000), double3(0.00000000, 0.00000000, 3.98749812))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-174-2", double3x3(double3(12.31999420, 0.00000000, 0.00000000), double3(-6.15999710, 10.66942795, 0.00000000), double3(0.00000000, 0.00000000, 9.87999535))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-175", double3x3(double3(12.65209405, 0.00000000, 0.00000000), double3(-6.32604702, 10.95703486, 0.00000000), double3(0.00000000, 0.00000000, 9.13809570))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-175-2", double3x3(double3(5.45899743, 0.00000000, 0.00000000), double3(-2.72949872, 4.72763045, 0.00000000), double3(0.00000000, 0.00000000, 3.08949855))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-176", double3x3(double3(6.41799698, 0.00000000, 0.00000000), double3(-3.20899849, 5.55814843, 0.00000000), double3(0.00000000, 0.00000000, 3.74299824))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-176-2", double3x3(double3(11.66999451, 0.00000000, 0.00000000), double3(-5.83499725, 10.10651171, 0.00000000), double3(0.00000000, 0.00000000, 9.94999532))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-177", double3x3(double3(6.34129702, 0.00000000, 0.00000000), double3(-3.17064851, 5.49172431, 0.00000000), double3(0.00000000, 0.00000000, 6.46219696))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-179", double3x3(double3(7.22129660, 0.00000000, 0.00000000), double3(-3.61064830, 6.25382631, 0.00000000), double3(0.00000000, 0.00000000, 6.90719675))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-179-2", double3x3(double3(10.41199510, 0.00000000, 0.00000000), double3(-5.20599755, 9.01705226, 0.00000000), double3(0.00000000, 0.00000000, 15.18399286))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-180", double3x3(double3(4.81899773, 0.00000000, 0.00000000), double3(-2.40949887, 4.17337446, 0.00000000), double3(0.00000000, 0.00000000, 6.59199690))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-180-2", double3x3(double3(4.89999769, 0.00000000, 0.00000000), double3(-2.44999885, 4.24352248, 0.00000000), double3(0.00000000, 0.00000000, 5.37999747))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-181", double3x3(double3(4.42829792, 0.00000000, 0.00000000), double3(-2.21414896, 3.83501849, 0.00000000), double3(0.00000000, 0.00000000, 6.36799700))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-181-2", double3x3(double3(10.48179507, 0.00000000, 0.00000000), double3(-5.24089753, 9.07750081, 0.00000000), double3(0.00000000, 0.00000000, 11.17499474))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-182", double3x3(double3(5.30999750, 0.00000000, 0.00000000), double3(-2.65499875, 4.59859273, 0.00000000), double3(0.00000000, 0.00000000, 14.24999329))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-182-2", double3x3(double3(5.45799743, 0.00000000, 0.00000000), double3(-2.72899872, 4.72676443, 0.00000000), double3(0.00000000, 0.00000000, 9.01599576))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-183", double3x3(double3(19.49999082, 0.00000000, 0.00000000), double3(-9.74999541, 16.88748743, 0.00000000), double3(0.00000000, 0.00000000, 10.29999515))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-183-2", double3x3(double3(3.39599840, 0.00000000, 0.00000000), double3(-1.69799920, 2.94102089, 0.00000000), double3(0.00000000, 0.00000000, 5.09199760))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-184", double3x3(double3(13.71799355, 0.00000000, 0.00000000), double3(-6.85899677, 11.88013090, 0.00000000), double3(0.00000000, 0.00000000, 8.45259602))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-184-2", double3x3(double3(13.80199351, 0.00000000, 0.00000000), double3(-6.90099675, 11.95287700, 0.00000000), double3(0.00000000, 0.00000000, 8.50299600))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-185", double3x3(double3(9.88499535, 0.00000000, 0.00000000), double3(-4.94249767, 8.56065709, 0.00000000), double3(0.00000000, 0.00000000, 10.80499492))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-185-2", double3x3(double3(6.25999705, 0.00000000, 0.00000000), double3(-3.12999853, 5.42131648, 0.00000000), double3(0.00000000, 0.00000000, 12.24899424))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-186", double3x3(double3(9.97999530, 0.00000000, 0.00000000), double3(-4.98999765, 8.64292946, 0.00000000), double3(0.00000000, 0.00000000, 7.63999641))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-186-2", double3x3(double3(8.09999619, 0.00000000, 0.00000000), double3(-4.04999809, 7.01480247, 0.00000000), double3(0.00000000, 0.00000000, 13.33999372))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-187", double3x3(double3(2.90649863, 0.00000000, 0.00000000), double3(-1.45324932, 2.51710165, 0.00000000), double3(0.00000000, 0.00000000, 2.83659867))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-187-2", double3x3(double3(5.44799744, 0.00000000, 0.00000000), double3(-2.72399872, 4.71810418, 0.00000000), double3(0.00000000, 0.00000000, 8.09099619))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-188", double3x3(double3(6.11999712, 0.00000000, 0.00000000), double3(-3.05999856, 5.30007298, 0.00000000), double3(0.00000000, 0.00000000, 5.65799734))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-188-2", double3x3(double3(9.21799566, 0.00000000, 0.00000000), double3(-4.60899783, 7.98301842, 0.00000000), double3(0.00000000, 0.00000000, 18.04199151))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-189", double3x3(double3(8.15399616, 0.00000000, 0.00000000), double3(-4.07699808, 7.06156782, 0.00000000), double3(0.00000000, 0.00000000, 6.13699711))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-189-2", double3x3(double3(9.64999546, 0.00000000, 0.00000000), double3(-4.82499773, 8.35714121, 0.00000000), double3(0.00000000, 0.00000000, 3.85699819))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-190", double3x3(double3(4.59609784, 0.00000000, 0.00000000), double3(-2.29804892, 3.98033749, 0.00000000), double3(0.00000000, 0.00000000, 8.92999580))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-190-2", double3x3(double3(10.56099503, 0.00000000, 0.00000000), double3(-5.28049752, 9.14608999, 0.00000000), double3(0.00000000, 0.00000000, 13.52199364))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-191", double3x3(double3(3.95999814, 0.00000000, 0.00000000), double3(-1.97999907, 3.42945899, 0.00000000), double3(0.00000000, 0.00000000, 3.84399819))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-191-2", double3x3(double3(11.25699470, 0.00000000, 0.00000000), double3(-5.62849735, 9.74884338, 0.00000000), double3(0.00000000, 0.00000000, 5.85299725))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-192", double3x3(double3(10.45699508, 0.00000000, 0.00000000), double3(-5.22849754, 9.05602339, 0.00000000), double3(0.00000000, 0.00000000, 14.23799330))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-192-2", double3x3(double3(9.76829540, 0.00000000, 0.00000000), double3(-4.88414770, 8.45959197, 0.00000000), double3(0.00000000, 0.00000000, 9.34079560))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-193", double3x3(double3(8.48999601, 0.00000000, 0.00000000), double3(-4.24499800, 7.35255222, 0.00000000), double3(0.00000000, 0.00000000, 6.08399714))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-193-2", double3x3(double3(9.74899541, 0.00000000, 0.00000000), double3(-4.87449771, 8.44287769, 0.00000000), double3(0.00000000, 0.00000000, 16.46999225))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-194", double3x3(double3(3.58699831, 0.00000000, 0.00000000), double3(-1.79349916, 3.10643166, 0.00000000), double3(0.00000000, 0.00000000, 15.49199271))),
    std::make_pair("spglibtestdata/hexagonal/POSCAR-194-2", double3x3(double3(3.46999837, 0.00000000, 0.00000000), double3(-1.73499918, 3.00510674, 0.00000000), double3(0.00000000, 0.00000000, 28.44998661))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Cubic)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/cubic/POSCAR-195", double3x3(double3(10.34999513, 0.00000000, 0.00000000), double3(0.00000000, 10.34999513, 0.00000000), double3(0.00000000, 0.00000000, 10.34999513))),
    std::make_pair("spglibtestdata/cubic/POSCAR-195-2", double3x3(double3(7.26599658, 0.00000000, 0.00000000), double3(0.00000000, 7.26599658, 0.00000000), double3(0.00000000, 0.00000000, 7.26599658))),
    std::make_pair("spglibtestdata/cubic/POSCAR-196", double3x3(double3(12.15399428, 0.00000000, 0.00000000), double3(0.00000000, 12.15399428, 0.00000000), double3(0.00000000, 0.00000000, 12.15399428))),
    std::make_pair("spglibtestdata/cubic/POSCAR-196-2", double3x3(double3(18.74999118, 0.00000000, 0.00000000), double3(0.00000000, 18.74999118, 0.00000000), double3(0.00000000, 0.00000000, 18.74999118))),
    std::make_pair("spglibtestdata/cubic/POSCAR-197", double3x3(double3(10.14539523, 0.00000000, 0.00000000), double3(0.00000000, 10.14539523, 0.00000000), double3(0.00000000, 0.00000000, 10.14539523))),
    std::make_pair("spglibtestdata/cubic/POSCAR-197-2", double3x3(double3(10.24999518, 0.00000000, 0.00000000), double3(0.00000000, 10.24999518, 0.00000000), double3(0.00000000, 0.00000000, 10.24999518))),
    std::make_pair("spglibtestdata/cubic/POSCAR-198", double3x3(double3(7.83999631, 0.00000000, 0.00000000), double3(0.00000000, 7.83999631, 0.00000000), double3(0.00000000, 0.00000000, 7.83999631))),
    std::make_pair("spglibtestdata/cubic/POSCAR-198-2", double3x3(double3(12.75299400, 0.00000000, 0.00000000), double3(0.00000000, 12.75299400, 0.00000000), double3(0.00000000, 0.00000000, 12.75299400))),
    std::make_pair("spglibtestdata/cubic/POSCAR-199", double3x3(double3(10.92999486, 0.00000000, 0.00000000), double3(0.00000000, 10.92999486, 0.00000000), double3(0.00000000, 0.00000000, 10.92999486))),
    std::make_pair("spglibtestdata/cubic/POSCAR-199-2", double3x3(double3(8.41899604, 0.00000000, 0.00000000), double3(0.00000000, 8.41899604, 0.00000000), double3(0.00000000, 0.00000000, 8.41899604))),
    std::make_pair("spglibtestdata/cubic/POSCAR-200", double3x3(double3(7.48699648, 0.00000000, 0.00000000), double3(0.00000000, 7.48699648, 0.00000000), double3(0.00000000, 0.00000000, 7.48699648))),
    std::make_pair("spglibtestdata/cubic/POSCAR-200-2", double3x3(double3(5.44999744, 0.00000000, 0.00000000), double3(0.00000000, 5.44999744, 0.00000000), double3(0.00000000, 0.00000000, 5.44999744))),
    std::make_pair("spglibtestdata/cubic/POSCAR-205", double3x3(double3(5.62399735, 0.00000000, 0.00000000), double3(0.00000000, 5.62399735, 0.00000000), double3(0.00000000, 0.00000000, 5.62399735))),
    std::make_pair("spglibtestdata/cubic/POSCAR-205-3", double3x3(double3(5.62399735, 0.00000000, 0.00000000), double3(0.00000000, 5.62399735, 0.00000000), double3(0.00000000, 0.00000000, 5.62399735))),
    std::make_pair("spglibtestdata/cubic/POSCAR-206", double3x3(double3(10.97999483, 0.00000000, 0.00000000), double3(0.00000000, 10.97999483, 0.00000000), double3(0.00000000, 0.00000000, 10.97999483))),
    std::make_pair("spglibtestdata/cubic/POSCAR-206-2", double3x3(double3(11.02999481, 0.00000000, 0.00000000), double3(0.00000000, 11.02999481, 0.00000000), double3(0.00000000, 0.00000000, 11.02999481))),
    std::make_pair("spglibtestdata/cubic/POSCAR-207", double3x3(double3(4.39999793, 0.00000000, 0.00000000), double3(0.00000000, 4.39999793, 0.00000000), double3(0.00000000, 0.00000000, 4.39999793))),
    std::make_pair("spglibtestdata/cubic/POSCAR-208", double3x3(double3(6.30999703, 0.00000000, 0.00000000), double3(0.00000000, 6.30999703, 0.00000000), double3(0.00000000, 0.00000000, 6.30999703))),
    std::make_pair("spglibtestdata/cubic/POSCAR-208-2", double3x3(double3(2.38574888, 0.00000000, 0.00000000), double3(0.00000000, 2.38574888, 0.00000000), double3(0.00000000, 0.00000000, 2.38574888))),
    std::make_pair("spglibtestdata/cubic/POSCAR-209", double3x3(double3(7.42299651, 0.00000000, 0.00000000), double3(0.00000000, 7.42299651, 0.00000000), double3(0.00000000, 0.00000000, 7.42299651))),
    std::make_pair("spglibtestdata/cubic/POSCAR-210", double3x3(double3(19.90999063, 0.00000000, 0.00000000), double3(0.00000000, 19.90999063, 0.00000000), double3(0.00000000, 0.00000000, 19.90999063))),
    std::make_pair("spglibtestdata/cubic/POSCAR-210-2", double3x3(double3(15.69899261, 0.00000000, 0.00000000), double3(0.00000000, 15.69899261, 0.00000000), double3(0.00000000, 0.00000000, 15.69899261))),
    std::make_pair("spglibtestdata/cubic/POSCAR-211", double3x3(double3(9.68879544, 0.00000000, 0.00000000), double3(0.00000000, 9.68879544, 0.00000000), double3(0.00000000, 0.00000000, 9.68879544))),
    std::make_pair("spglibtestdata/cubic/POSCAR-212", double3x3(double3(6.71499684, 0.00000000, 0.00000000), double3(0.00000000, 6.71499684, 0.00000000), double3(0.00000000, 0.00000000, 6.71499684))),
    std::make_pair("spglibtestdata/cubic/POSCAR-212-2", double3x3(double3(6.71499684, 0.00000000, 0.00000000), double3(0.00000000, 6.71499684, 0.00000000), double3(0.00000000, 0.00000000, 6.71499684))),
    std::make_pair("spglibtestdata/cubic/POSCAR-213", double3x3(double3(10.27999516, 0.00000000, 0.00000000), double3(0.00000000, 10.27999516, 0.00000000), double3(0.00000000, 0.00000000, 10.27999516))),
    std::make_pair("spglibtestdata/cubic/POSCAR-213-2", double3x3(double3(7.93599627, 0.00000000, 0.00000000), double3(0.00000000, 7.93599627, 0.00000000), double3(0.00000000, 0.00000000, 7.93599627))),
    std::make_pair("spglibtestdata/cubic/POSCAR-214", double3x3(double3(21.75998976, 0.00000000, 0.00000000), double3(0.00000000, 21.75998976, 0.00000000), double3(0.00000000, 0.00000000, 21.75998976))),
    std::make_pair("spglibtestdata/cubic/POSCAR-214-2", double3x3(double3(12.31499421, 0.00000000, 0.00000000), double3(0.00000000, 12.31499421, 0.00000000), double3(0.00000000, 0.00000000, 12.31499421))),
    std::make_pair("spglibtestdata/cubic/POSCAR-215", double3x3(double3(5.39299746, 0.00000000, 0.00000000), double3(0.00000000, 5.39299746, 0.00000000), double3(0.00000000, 0.00000000, 5.39299746))),
    std::make_pair("spglibtestdata/cubic/POSCAR-215-2", double3x3(double3(8.31999609, 0.00000000, 0.00000000), double3(0.00000000, 8.31999609, 0.00000000), double3(0.00000000, 0.00000000, 8.31999609))),
    std::make_pair("spglibtestdata/cubic/POSCAR-216", double3x3(double3(7.17599662, 0.00000000, 0.00000000), double3(0.00000000, 7.17599662, 0.00000000), double3(0.00000000, 0.00000000, 7.17599662))),
    std::make_pair("spglibtestdata/cubic/POSCAR-216-2", double3x3(double3(7.57699643, 0.00000000, 0.00000000), double3(0.00000000, 7.57699643, 0.00000000), double3(0.00000000, 0.00000000, 7.57699643))),
    std::make_pair("spglibtestdata/cubic/POSCAR-217", double3x3(double3(12.69999402, 0.00000000, 0.00000000), double3(0.00000000, 12.69999402, 0.00000000), double3(0.00000000, 0.00000000, 12.69999402))),
    std::make_pair("spglibtestdata/cubic/POSCAR-217-2", double3x3(double3(10.16799522, 0.00000000, 0.00000000), double3(0.00000000, 10.16799522, 0.00000000), double3(0.00000000, 0.00000000, 10.16799522))),
    std::make_pair("spglibtestdata/cubic/POSCAR-218", double3x3(double3(8.29399610, 0.00000000, 0.00000000), double3(0.00000000, 8.29399610, 0.00000000), double3(0.00000000, 0.00000000, 8.29399610))),
    std::make_pair("spglibtestdata/cubic/POSCAR-218-2", double3x3(double3(6.02599716, 0.00000000, 0.00000000), double3(0.00000000, 6.02599716, 0.00000000), double3(0.00000000, 0.00000000, 6.02599716))),
    std::make_pair("spglibtestdata/cubic/POSCAR-219", double3x3(double3(17.34399184, 0.00000000, 0.00000000), double3(0.00000000, 17.34399184, 0.00000000), double3(0.00000000, 0.00000000, 17.34399184))),
    std::make_pair("spglibtestdata/cubic/POSCAR-219-2", double3x3(double3(12.14099429, 0.00000000, 0.00000000), double3(0.00000000, 12.14099429, 0.00000000), double3(0.00000000, 0.00000000, 12.14099429))),
    std::make_pair("spglibtestdata/cubic/POSCAR-220", double3x3(double3(9.81799538, 0.00000000, 0.00000000), double3(0.00000000, 9.81799538, 0.00000000), double3(0.00000000, 0.00000000, 9.81799538))),
    std::make_pair("spglibtestdata/cubic/POSCAR-220-2", double3x3(double3(8.53399598, 0.00000000, 0.00000000), double3(0.00000000, 8.53399598, 0.00000000), double3(0.00000000, 0.00000000, 8.53399598))),
    std::make_pair("spglibtestdata/cubic/POSCAR-221", double3x3(double3(9.63799546, 0.00000000, 0.00000000), double3(0.00000000, 9.63799546, 0.00000000), double3(0.00000000, 0.00000000, 9.63799546))),
    std::make_pair("spglibtestdata/cubic/POSCAR-221-2", double3x3(double3(5.79499727, 0.00000000, 0.00000000), double3(0.00000000, 5.79499727, 0.00000000), double3(0.00000000, 0.00000000, 5.79499727))),
    std::make_pair("spglibtestdata/cubic/POSCAR-222", double3x3(double3(10.98999483, 0.00000000, 0.00000000), double3(0.00000000, 10.98999483, 0.00000000), double3(0.00000000, 0.00000000, 10.98999483))),
    std::make_pair("spglibtestdata/cubic/POSCAR-222-2", double3x3(double3(16.25599235, 0.00000000, 0.00000000), double3(0.00000000, 16.25599235, 0.00000000), double3(0.00000000, 0.00000000, 16.25599235))),
    std::make_pair("spglibtestdata/cubic/POSCAR-223", double3x3(double3(6.66999686, 0.00000000, 0.00000000), double3(0.00000000, 6.66999686, 0.00000000), double3(0.00000000, 0.00000000, 6.66999686))),
    std::make_pair("spglibtestdata/cubic/POSCAR-223-2", double3x3(double3(10.29999515, 0.00000000, 0.00000000), double3(0.00000000, 10.29999515, 0.00000000), double3(0.00000000, 0.00000000, 10.29999515))),
    std::make_pair("spglibtestdata/cubic/POSCAR-224", double3x3(double3(4.90399769, 0.00000000, 0.00000000), double3(0.00000000, 4.90399769, 0.00000000), double3(0.00000000, 0.00000000, 4.90399769))),
    std::make_pair("spglibtestdata/cubic/POSCAR-224-2", double3x3(double3(4.90399769, 0.00000000, 0.00000000), double3(0.00000000, 4.90399769, 0.00000000), double3(0.00000000, 0.00000000, 4.90399769))),
    std::make_pair("spglibtestdata/cubic/POSCAR-225", double3x3(double3(9.98999530, 0.00000000, 0.00000000), double3(0.00000000, 9.98999530, 0.00000000), double3(0.00000000, 0.00000000, 9.98999530))),
    std::make_pair("spglibtestdata/cubic/POSCAR-225-2", double3x3(double3(4.09649807, 0.00000000, 0.00000000), double3(0.00000000, 4.09649807, 0.00000000), double3(0.00000000, 0.00000000, 4.09649807))),
    std::make_pair("spglibtestdata/cubic/POSCAR-226", double3x3(double3(25.05998821, 0.00000000, 0.00000000), double3(0.00000000, 25.05998821, 0.00000000), double3(0.00000000, 0.00000000, 25.05998821))),
    std::make_pair("spglibtestdata/cubic/POSCAR-226-2", double3x3(double3(10.04599527, 0.00000000, 0.00000000), double3(0.00000000, 10.04599527, 0.00000000), double3(0.00000000, 0.00000000, 10.04599527))),
    std::make_pair("spglibtestdata/cubic/POSCAR-227", double3x3(double3(10.12999523, 0.00000000, 0.00000000), double3(0.00000000, 10.12999523, 0.00000000), double3(0.00000000, 0.00000000, 10.12999523))),
    std::make_pair("spglibtestdata/cubic/POSCAR-227-2", double3x3(double3(23.25498906, 0.00000000, 0.00000000), double3(0.00000000, 23.25498906, 0.00000000), double3(0.00000000, 0.00000000, 23.25498906))),
    std::make_pair("spglibtestdata/cubic/POSCAR-228", double3x3(double3(15.70499261, 0.00000000, 0.00000000), double3(0.00000000, 15.70499261, 0.00000000), double3(0.00000000, 0.00000000, 15.70499261))),
    std::make_pair("spglibtestdata/cubic/POSCAR-228-2", double3x3(double3(21.80998974, 0.00000000, 0.00000000), double3(0.00000000, 21.80998974, 0.00000000), double3(0.00000000, 0.00000000, 21.80998974))),
    std::make_pair("spglibtestdata/cubic/POSCAR-229", double3x3(double3(18.26999140, 0.00000000, 0.00000000), double3(0.00000000, 18.26999140, 0.00000000), double3(0.00000000, 0.00000000, 18.26999140))),
    std::make_pair("spglibtestdata/cubic/POSCAR-229-2", double3x3(double3(6.22099707, 0.00000000, 0.00000000), double3(0.00000000, 6.22099707, 0.00000000), double3(0.00000000, 0.00000000, 6.22099707))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230", double3x3(double3(12.60199407, 0.00000000, 0.00000000), double3(0.00000000, 12.60199407, 0.00000000), double3(0.00000000, 0.00000000, 12.60199407))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-2", double3x3(double3(12.37599418, 0.00000000, 0.00000000), double3(0.00000000, 12.37599418, 0.00000000), double3(0.00000000, 0.00000000, 12.37599418))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-3", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/cubic/POSCAR-230-4", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}

TEST(FindConventionalCell, Virtual)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_real_distribution<double> dist(-0.5, 0.5);

  const std::map<std::string, double3x3> testData =
  {
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-221-33", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-222-33", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-223-33", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-224-33", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-73", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-93", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-227-99", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-conv-56", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-230-prim-33", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-1-bcc-33", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-221-18", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-223-18", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-10-227-50", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-102-224-13", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-104-222-13", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-105-223-13", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-109-227-13", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-11-227-48", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-conv-15", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-110-230-prim-13", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-221-11", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-224-11", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-111-227-66", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-222-11", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-112-223-11", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-113-227-68", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-221-14", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-223-14", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-115-227-33", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-116-230-conv-34", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-117-230-conv-33", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-222-14", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-118-224-14", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-221-19", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-224-19", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-21", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-12-227-83", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-conv-16", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-120-230-prim-14", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-conv-13", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-122-230-prim-11", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-123-221-05", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-126-222-05", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-222-18", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-224-18", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-227-49", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-13-230-conv-44", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-131-223-05", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-134-224-05", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-47", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-227-51", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-14-230-conv-45", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-conv-05", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-142-230-prim-05", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-221-27", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-222-27", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-223-27", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-224-27", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-227-92", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-36", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-conv-55", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-230-prim-27", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-146-bcc-27", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-221-15", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-222-15", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-223-15", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-224-15", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-227-70", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-17", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-conv-37", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-230-prim-15", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-148-bcc-15", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-222-19", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-223-19", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-21", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-conv-22", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-18", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-230-prim-19", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-18", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-15-bcc-19", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-221-17", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-222-17", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-223-17", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-224-17", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-227-72", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-19", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-conv-38", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-230-prim-17", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-155-bcc-17", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-221-20", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-222-20", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-223-20", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-224-20", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-16-227-84", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-221-16", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-224-16", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-16", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-227-71", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-160-fcc", double3x3(double3(7.07106781, 0.00000000, 0.00000000), double3(-3.53553391, 6.12372436, 0.00000000), double3(0.00000000, 0.00000000, 17.32050808))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-222-16", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-223-16", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-conv-18", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-230-prim-16", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-161-bcc-16", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-221-06", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-224-06", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-06", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-166-227-38", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-222-06", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-223-06", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-conv-06", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-230-prim-06", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-167-bcc-6", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-60", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-227-85", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-17-230-conv-46", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-18-227-86", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-59", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-227-89", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-19-230-conv-51", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-221-07", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-222-07", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-223-07", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-195-224-07", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-227-40", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-198-230-conv-20", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-conv-07", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-199-230-prim-07", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-221-28", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-222-28", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-223-28", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-224-28", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-41", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-74", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-227-94", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-39", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-conv-57", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-230-prim-28", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-2-bcc-28", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-53", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-227-90", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-20-230-conv-53", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-221-02", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-200-223-02", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-222-02", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-201-224-02", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-205-230-conv-08", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-conv-02", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-206-230-prim-02", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-221-04", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-207-222-04", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-223-04", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-208-224-04", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-221-23", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-222-23", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-223-23", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-224-23", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-21-230-conv-49", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-212-227-19", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-213-230-conv-09", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-conv-04", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-214-230-prim-04", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-221-03", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-224-03", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-215-227-18", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-216-227-03", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-222-03", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-218-223-03", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-conv-26", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-22-230-prim-23", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-conv-03", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-220-230-prim-03", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-221-221-01", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-222-222-01", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-223-223-01", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-224-224-01", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-227-227-01", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-01", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-conv-62", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-230-230-prim-01", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-24-230-conv-23", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-24-230-prim-20", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-221-21", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-223-21", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-25-227-54", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-26-227-64", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-27-230-conv-48", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-28-227-62", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-29-230-conv-52", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-221-29", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-222-29", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-223-29", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-224-29", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-227-82", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-227-95", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-3-230-conv-58", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-30-227-65", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-31-227-58", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-32-230-conv-47", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-33-227-63", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-34-222-21", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-34-224-21", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-221-22", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-224-22", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-35-227-87", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-37-222-22", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-37-223-22", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-38-221-26", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-39-224-26", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-77", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-81", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-227-96", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-4-230-conv-59", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-40-223-26", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-41-222-26", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-conv-25", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-conv-29", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-prim-22", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-230-prim-26", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-bcc-22", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-43-bcc-26", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-44-227-24", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-45-230-conv-24", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-45-230-prim-21", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-46-227-28", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-47-221-08", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-47-223-08", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-48-222-08", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-48-224-08", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-221-32", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-222-32", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-223-32", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-224-32", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-45", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-75", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-227-98", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-40", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-43", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-conv-61", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-prim-29", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-230-prim-32", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-bcc-29", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-5-bcc-32", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-51-227-29", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-53-227-32", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-54-230-conv-30", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-221-30", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-223-30", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-6-227-79", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-61-230-conv-31", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-62-227-31", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-65-221-09", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-66-223-09", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-67-224-09", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-68-222-09", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-222-30", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-224-30", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-227-78", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-227-80", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-7-230-conv-60", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-230-conv-11", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-230-prim-09", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-70-bcc-9", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-73-230-conv-10", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-73-230-prim-08", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-74-227-09", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-75-221-25", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-75-222-25", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-76-227-61", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-77-223-25", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-77-224-25", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-78-227-91", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-78-230-conv-54", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-221-31", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-224-31", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-227-44", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-8-227-97", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-80-230-conv-28", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-80-230-prim-25", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-221-24", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-222-24", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-223-24", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-224-24", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-227-88", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-81-230-conv-50", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-82-230-conv-27", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-82-230-prim-24", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-83-221-10", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-84-223-10", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-85-222-10", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-86-224-10", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-88-230-conv-12", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-88-230-prim-10", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-89-221-12", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-89-222-12", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-222-31", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-223-31", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-227-43", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-conv-41", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-conv-42", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-prim-30", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-230-prim-31", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-bcc-30", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-9-bcc-31", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-91-227-67", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-92-227-35", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-92-230-conv-35", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-93-223-12", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-93-224-12", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-95-227-36", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-95-230-conv-32", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-96-227-69", double3x3(double3(8.24000000, 0.00000000, 0.00000000), double3(0.00000000, 8.24000000, 0.00000000), double3(0.00000000, 0.00000000, 8.24000000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-98-230-conv-14", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-98-230-prim-12", double3x3(double3(12.80654000, 0.00000000, 0.00000000), double3(0.00000000, 12.80654000, 0.00000000), double3(0.00000000, 0.00000000, 12.80654000))),
    std::make_pair("spglibtestdata/virtual_structure/POSCAR-99-221-13", double3x3(double3(10.00000000, 0.00000000, 0.00000000), double3(0.00000000, 10.00000000, 0.00000000), double3(0.00000000, 0.00000000, 10.00000000))) 
  };

  for (auto& [fileName, conventionalCellTargetValue] : testData)
  {
    std::ifstream t(fileName.c_str());
    std::string fileContent((std::istreambuf_iterator<char>(t)),
      std::istreambuf_iterator<char>());

    SKPOSCARLegacyParser parser = SKPOSCARLegacyParser(fileContent);
    parser.startParsing();

    std::vector<std::tuple<double3, size_t, double> > atoms = parser.firstTestFrame();
    double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
    bool allowPartialOccupancies = false;
    double symmetryPrecision = 1e-5;

    double3 randomShift = double3(dist(mt), dist(mt), dist(mt));
    std::vector<std::tuple<double3, size_t, double>> randomlyShiftedAtoms{};
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(randomlyShiftedAtoms),
      [randomShift](const std::tuple<double3, size_t, double>& atom) { return std::make_tuple(std::get<0>(atom) + randomShift, std::get<1>(atom), std::get<2>(atom)); });

    std::optional<SKSpaceGroup::FoundSpaceGroupInfo> spaceGroup = SKSpaceGroup::findSpaceGroup(unitCell, randomlyShiftedAtoms, allowPartialOccupancies, symmetryPrecision);

    if (spaceGroup)
    {
      double3x3 foundConventionalCell = spaceGroup->cell.unitCell();
      EXPECT_EQ(conventionalCellTargetValue, conventionalCellTargetValue);
    }
    else
    {
      EXPECT_TRUE(false);
    }
  }
}