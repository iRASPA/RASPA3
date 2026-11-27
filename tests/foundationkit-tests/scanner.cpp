#include <gtest/gtest.h>

import std;

import characterset;
import scanner;

// A line of a file with Windows line ends is the text before its carriage return: the '\r' stayed on the line (on the
// last value of every line of such a file), the last line included.
TEST(scanner, a_line_with_a_windows_line_end_has_no_carriage_return)
{
  Scanner scanner("first line\r\nsecond\r\nlast\r", CharacterSet::whitespaceCharacterSet());
  std::string line;
  ASSERT_TRUE(scanner.scanLine(line));
  EXPECT_EQ(line, "first line");
  ASSERT_TRUE(scanner.scanLine(line));
  EXPECT_EQ(line, "second");
  ASSERT_TRUE(scanner.scanLine(line));
  EXPECT_EQ(line, "last");
  EXPECT_FALSE(scanner.scanLine(line));
}
