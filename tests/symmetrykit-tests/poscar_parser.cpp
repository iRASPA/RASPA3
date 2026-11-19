#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import skposcarparser;
import skstructure;
import skatom;
import skcell;
import bool3;

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
  SKPOSCARParser parser = SKPOSCARParser::fromFile("testdata/poscar/vasp4_counts.poscar");
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  EXPECT_EQ(totalAtomCount(atoms), 6U);
  EXPECT_EQ(uniqueTypeCount(atoms), 3U);
  EXPECT_DOUBLE_EQ(parser.movies().front().front()->cell->unitCell().ax, 10.0);
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
  SKPOSCARParser parser = SKPOSCARParser::fromFile("testdata/poscar/vasp4_element_like_comment.poscar");
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  EXPECT_EQ(totalAtomCount(atoms), 4U);
  EXPECT_EQ(uniqueTypeCount(atoms), 4U);
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

// VASP reads the first character of the coordinate-mode line alone: "Cart" is Cartesian as "Cartesian" is.
TEST(POSCARParser, CoordinateModeIsItsFirstCharacter)
{
  // (the frame is in fractional coordinates: 2 Angstrom along an axis of 4 is a half, and a direct 2 is 2)
  for (const auto& [mode, x] : {std::pair{"Cart", 0.5}, std::pair{"cartesian", 0.5}, std::pair{"K", 0.5},
                                std::pair{"Direct", 2.0}, std::pair{"dir", 2.0}})
  {
    const std::string poscar = std::format(R"(coordinate mode
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Fe
1
{}
2.0 0.0 0.0
)",
                                           mode);
    SKPOSCARParser parser(poscar);
    parser.startParsing();
    const auto atoms = parser.firstTestFrame();
    ASSERT_EQ(totalAtomCount(atoms), 1U) << mode;
    EXPECT_NEAR(std::get<0>(atoms.front()).x, x, 1e-12) << mode;
  }
}

TEST(POSCARParser, StripsInlineCommentsOnAtomLines)
{
  SKPOSCARParser parser = SKPOSCARParser::fromFile("testdata/poscar/vasp4_counts.poscar");
  parser.startParsing();
  EXPECT_EQ(totalAtomCount(parser.firstTestFrame()), 6U);
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

// VASP species labels with the POTCAR variant ('Fe_pv', 'O_s', and with the hash VASP 5.4 appends): the element is
// the part before it. The whole label missed the element table and took the species index as atomic number, so iron
// came out as hydrogen and oxygen as helium.
TEST(POSCARParser, ParsesSpeciesLabelsWithPotcarVariants)
{
  const std::string poscar = R"(FeO with POTCAR labels
1.0
4.3 0.0 0.0
0.0 4.3 0.0
0.0 0.0 4.3
Fe_pv O_s/5ae3c0
1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
)";

  SKPOSCARParser parser(poscar);
  parser.startParsing();

  const auto atoms = parser.firstTestFrame();
  ASSERT_EQ(totalAtomCount(atoms), 2U);
  EXPECT_EQ(std::get<1>(atoms[0]), 26U);
  EXPECT_EQ(std::get<1>(atoms[1]), 8U);
}

namespace
{
std::string poscarWith(const std::string& scale, const std::string& firstAtomLine,
                       const std::string& mode = "Direct")
{
  return "comment\n" + scale + "\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nSn Cl\n1 1\n" + mode + "\n" + firstAtomLine +
         "\n0.5 0.5 0.5 T T T\n";
}
}  // namespace

// A number is the whole token, with the 'D' exponent Fortran writes: "0.5D-01" was read as 0.5, and "nan" and
// "inf" as numbers.
TEST(POSCARParser, NumbersAreReadWhole)
{
  SKPOSCARParser parser(poscarWith("0.10000000D+01", "0.50000000D-01 +2.5d-1 1.0E-1"));
  parser.startParsing();
  const auto atoms = parser.firstTestFrame();
  ASSERT_EQ(atoms.size(), 2U);
  EXPECT_DOUBLE_EQ(std::get<0>(atoms[0]).x, 0.05);
  EXPECT_DOUBLE_EQ(std::get<0>(atoms[0]).y, 0.25);
  EXPECT_DOUBLE_EQ(std::get<0>(atoms[0]).z, 0.1);
  EXPECT_DOUBLE_EQ(parser.movies().front().front()->cell->unitCell().ax, 4.0);

  for (const auto& [scale, line] :
       {std::pair{"nan", "0 0 0"}, std::pair{"inf", "0 0 0"}, std::pair{"1.0", "nan 0 0"}, std::pair{"1.0", "0 inf 0"},
        std::pair{"1.0", "0.5x 0 0"}, std::pair{"1.0", "0 0 1e999"}, std::pair{"1.0abc", "0 0 0"}})
  {
    SKPOSCARParser refused(poscarWith(scale, line));
    EXPECT_THROW(refused.startParsing(), std::runtime_error) << scale << " / " << line;
  }
}

// The selective-dynamics flags are kept: 'F' holds the atom along that lattice vector.
TEST(POSCARParser, SelectiveDynamicsFlagsFixTheAtom)
{
  SKPOSCARParser parser(poscarWith("1.0", "0.0 0.0 0.0 F T .FALSE.", "Selective dynamics\nDirect"));
  parser.startParsing();
  const auto& atoms = parser.movies().front().front()->atoms;
  ASSERT_EQ(atoms.size(), 2U);
  EXPECT_TRUE(atoms[0]->isFixed().x);
  EXPECT_FALSE(atoms[0]->isFixed().y);
  EXPECT_TRUE(atoms[0]->isFixed().z);
  EXPECT_FALSE(atoms[1]->isFixed().x || atoms[1]->isFixed().y || atoms[1]->isFixed().z);

  SKPOSCARParser refused(poscarWith("1.0", "0.0 0.0 0.0 F T 7", "Selective dynamics\nDirect"));
  EXPECT_THROW(refused.startParsing(), std::runtime_error);
}

// Lines that end in a carriage return, with or without a line feed, are the lines of the file.
TEST(POSCARParser, CarriageReturnLineEnds)
{
  for (const std::string lineEnd : {"\r", "\r\n"})
  {
    std::string content = poscarWith("1.0", "0.25 0.5 0.75");
    std::string converted;
    for (const char character : content)
    {
      if (character == '\n')
        converted += lineEnd;
      else
        converted.push_back(character);
    }
    SKPOSCARParser parser(converted);
    ASSERT_NO_THROW(parser.startParsing());
    const auto atoms = parser.firstTestFrame();
    ASSERT_EQ(atoms.size(), 2U);
    EXPECT_DOUBLE_EQ(std::get<0>(atoms[0]).y, 0.5);
    EXPECT_EQ(std::get<1>(atoms[0]), 50U);
  }
}

// More unnamed species than there are elements: the ones past the table have no element (their number in the
// file was taken for one, and read past the end of the table of elements).
TEST(POSCARParser, UnnamedSpeciesPastTheTableOfElements)
{
  std::string content = "comment\n1.0\n40.0 0.0 0.0\n0.0 40.0 0.0\n0.0 0.0 40.0\n";
  const std::size_t numberOfSpecies = 125;
  for (std::size_t i = 0; i < numberOfSpecies; ++i) content += "1 ";
  content += "\nDirect\n";
  for (std::size_t i = 0; i < numberOfSpecies; ++i) content += std::format("{} 0.0 0.0\n", i / 200.0);

  SKPOSCARParser parser(content);
  ASSERT_NO_THROW(parser.startParsing());
  const auto atoms = parser.firstTestFrame();
  ASSERT_EQ(atoms.size(), numberOfSpecies);
  EXPECT_EQ(std::get<1>(atoms[0]), 1U);
  EXPECT_EQ(std::get<1>(atoms[117]), 118U);
  EXPECT_EQ(std::get<1>(atoms[118]), 0U);
  EXPECT_EQ(std::get<1>(atoms[124]), 0U);
}
