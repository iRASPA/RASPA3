#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import skposcarparser;

namespace
{

std::size_t totalAtomCount(const std::vector<std::tuple<double3, std::size_t, double>>& atoms)
{
  return atoms.size();
}

std::size_t uniqueTypeCount(const std::vector<std::tuple<double3, std::size_t, double>>& atoms)
{
  std::set<std::size_t> types{};
  for (const auto& atom : atoms)
  {
    types.insert(std::get<1>(atom));
  }
  return types.size();
}

}  // namespace

TEST(POSCARParser, ParsesVasp4SpeciesCounts)
{
  SKPOSCARParser parser = SKPOSCARParser::fromFile("spglibtestdata/cubic/POSCAR-225");
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  EXPECT_EQ(totalAtomCount(atoms), 36U);
  EXPECT_EQ(uniqueTypeCount(atoms), 3U);
  EXPECT_DOUBLE_EQ(parser.movies().front().front()->cell->unitCell().ax, 9.9899952992877097);
}

TEST(POSCARParser, ParsesVasp5ElementSymbols)
{
  const std::string poscar = R"(VASP5 example
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Sn Cl K
1 2 1
Direct
0.0 0.0 0.0
0.25 0.25 0.25
0.5 0.5 0.5
0.75 0.75 0.75
)";

  SKPOSCARParser parser(poscar);
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  EXPECT_EQ(totalAtomCount(atoms), 4U);
  EXPECT_EQ(uniqueTypeCount(atoms), 3U);
  EXPECT_EQ(std::get<1>(atoms[0]), 50U);
  EXPECT_EQ(std::get<1>(atoms[1]), 17U);
  EXPECT_EQ(std::get<1>(atoms[2]), 17U);
  EXPECT_EQ(std::get<1>(atoms[3]), 19U);
}

TEST(POSCARParser, ParsesSpglibVasp4WithCommentLineContainingElementLikeText)
{
  SKPOSCARParser parser =
      SKPOSCARParser::fromFile("spglibtestdata/virtual_structure/POSCAR-230-230-conv-62");
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  EXPECT_EQ(totalAtomCount(atoms), 208U);
  EXPECT_EQ(uniqueTypeCount(atoms), 208U);
}

TEST(POSCARParser, AppliesScaleFactor)
{
  const std::string poscar = R"(scaled cell
2.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
1
Direct
0.0 0.0 0.0
)";

  SKPOSCARParser parser(poscar);
  parser.startParsing();

  const double3x3 unitCell = parser.movies().front().front()->cell->unitCell();
  EXPECT_DOUBLE_EQ(unitCell.ax, 2.0);
  EXPECT_DOUBLE_EQ(unitCell.by, 2.0);
  EXPECT_DOUBLE_EQ(unitCell.cz, 2.0);
}

TEST(POSCARParser, ParsesSelectiveDynamicsAndCartesianMode)
{
  const std::string poscar = R"(selective dynamics example
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Fe
1
Selective dynamics
Cartesian
0.0 0.0 0.0 T T F
)";

  SKPOSCARParser parser(poscar);
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  ASSERT_EQ(totalAtomCount(atoms), 1U);
  EXPECT_EQ(std::get<1>(atoms.front()), 26U);
  EXPECT_NEAR(std::get<0>(atoms.front()).x, 0.0, 1e-12);
}

TEST(POSCARParser, StripsInlineCommentsOnAtomLines)
{
  SKPOSCARParser parser = SKPOSCARParser::fromFile("spglibtestdata/cubic/POSCAR-225");
  parser.startParsing();
  EXPECT_EQ(totalAtomCount(parser.firstTestFrame()), 36U);
}

TEST(POSCARParser, RejectsMismatchedVasp5SpeciesLines)
{
  const std::string poscar = R"(bad vasp5
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Sn Cl
1
Direct
0.0 0.0 0.0
)";

  SKPOSCARParser parser(poscar);
  EXPECT_THROW(parser.startParsing(), std::runtime_error);
}
