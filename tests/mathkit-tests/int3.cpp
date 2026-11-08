#include <gtest/gtest.h>

import std;

import int3;
import uint3;

TEST(int3, Test_less_than_is_lexicographic_and_usable_as_set_key)
{
  int3 a(1, 2, 3);
  int3 b(1, 3, -5);
  EXPECT_TRUE(a < b);
  EXPECT_FALSE(b < a);
  EXPECT_FALSE(a < a);

  std::set<int3> keys{b, a, int3(1, 2, 3), int3(-1, 0, 0)};
  EXPECT_EQ(keys.size(), 3uz);
  EXPECT_EQ(keys.begin()->x, -1);

  std::map<int3, int> counts;
  ++counts[a];
  ++counts[int3(1, 2, 3)];
  EXPECT_EQ(counts[a], 2);
}

TEST(uint3, Test_less_than_is_lexicographic_and_usable_as_set_key)
{
  uint3 a(1, 2, 3);
  uint3 b(1, 3, 0);
  EXPECT_TRUE(a < b);
  EXPECT_FALSE(b < a);
  EXPECT_FALSE(a < a);

  std::set<uint3> keys{b, a, uint3(1, 2, 3), uint3(0, 5, 5)};
  EXPECT_EQ(keys.size(), 3uz);
  EXPECT_EQ(keys.begin()->x, 0u);

  std::map<uint3, int> counts;
  ++counts[a];
  ++counts[uint3(1, 2, 3)];
  EXPECT_EQ(counts[a], 2);
}
