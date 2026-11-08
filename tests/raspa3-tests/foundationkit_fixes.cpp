#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import archive;
import characterset;
import getopt;
import indexpath;
import scanner;
import stringutils;

TEST(FOUNDATIONKIT_FIXES, scanner_reads_numbers_and_the_last_token_and_line)
{
  Scanner numbers("1.5 -7 2.25", CharacterSet::whitespaceAndNewlineCharacterSet());
  double x{};
  int n{};
  EXPECT_TRUE(numbers.scanDouble(x));
  EXPECT_DOUBLE_EQ(x, 1.5);
  EXPECT_TRUE(numbers.scanInt(n));
  EXPECT_EQ(n, -7);
  // the last token has no delimiter after it
  EXPECT_TRUE(numbers.scanDouble(x));
  EXPECT_DOUBLE_EQ(x, 2.25);
  EXPECT_FALSE(numbers.scanDouble(x));
  EXPECT_FALSE(Scanner("abc", CharacterSet::whitespaceAndNewlineCharacterSet()).scanDouble(x));

  Scanner tokens("alpha\tbeta", CharacterSet::whitespaceAndNewlineCharacterSet());
  std::string token;
  EXPECT_TRUE(tokens.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), token));
  EXPECT_EQ(token, "alpha");  // the tab is whitespace
  EXPECT_TRUE(tokens.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), token));
  EXPECT_EQ(token, "beta");
  EXPECT_TRUE(tokens.isAtEnd());

  Scanner lines("first\nlast", CharacterSet::whitespaceCharacterSet());
  std::string line;
  EXPECT_TRUE(lines.scanLine(line));
  EXPECT_EQ(line, "first");
  EXPECT_TRUE(lines.scanLine(line));
  EXPECT_EQ(line, "last");
  EXPECT_FALSE(lines.scanLine(line));
}

TEST(FOUNDATIONKIT_FIXES, whitespace_and_newline_set_holds_the_tab)
{
  const std::string set = CharacterSet::whitespaceAndNewlineCharacterSet().string();
  EXPECT_TRUE(set.contains('\t'));
  EXPECT_TRUE(set.contains(' '));
  EXPECT_TRUE(set.contains('\n'));
}

TEST(FOUNDATIONKIT_FIXES, equal_index_paths_compare_equal)
{
  EXPECT_TRUE(IndexPath(3) == IndexPath(3));
  EXPECT_FALSE(IndexPath(3) == IndexPath(4));
}

TEST(FOUNDATIONKIT_FIXES, an_optional_argument_takes_the_next_value_but_not_the_next_option)
{
  std::array<std::string, 5> storage{"raspa3", "-o", "value", "-o", "-v"};
  std::array<char*, 5> argv{};
  for (std::size_t i = 0; i < storage.size(); ++i) argv[i] = storage[i].data();

  std::vector<std::string> optionalValues;
  bool verbose = false;
  argparser::argparser parser(static_cast<int>(argv.size()), argv.data());
  parser.reg({"-o"}, "value", argparser::argparser::optional_argument, "",
             [&](const std::string& value) { optionalValues.push_back(value); });
  parser.reg({"-v"}, argparser::argparser::no_argument, "", [&](const std::string&) { verbose = true; });
  parser();

  ASSERT_EQ(optionalValues.size(), 2uz);
  EXPECT_EQ(optionalValues[0], "value");
  EXPECT_EQ(optionalValues[1], "");
  EXPECT_TRUE(verbose);
}

TEST(FOUNDATIONKIT_FIXES, trim_handles_empty_and_blank_strings)
{
  EXPECT_EQ(trim(""), "");
  EXPECT_EQ(trim("   "), "");
  EXPECT_EQ(trim("  a b \t"), "a b");
}

TEST(FOUNDATIONKIT_FIXES, reading_a_map_replaces_its_contents)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "map.bin";
  {
    std::ofstream stream(file, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << std::map<int, double>{{1, 1.5}};
    archive << std::unordered_map<int, double>{{2, 2.5}};
  }
  std::map<int, double> map{{7, 7.0}};
  std::unordered_map<int, double> unorderedMap{{8, 8.0}};
  {
    std::ifstream stream(file, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> map;
    archive >> unorderedMap;
  }
  EXPECT_EQ(map, (std::map<int, double>{{1, 1.5}}));
  EXPECT_EQ(unorderedMap, (std::unordered_map<int, double>{{2, 2.5}}));
}
