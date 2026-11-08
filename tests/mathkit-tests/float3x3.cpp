#include <gtest/gtest.h>

import std;

import float3;
import float3x3;

TEST(float3x3, Test_subscript_writes_into_the_matrix)
{
  float3x3 m{};
  m[0].x = 1.0f;
  m[1][2] = 2.0f;
  m[2] = float3(3.0f, 4.0f, 5.0f);

  EXPECT_EQ(m.m11, 1.0f);
  EXPECT_EQ(m.m32, 2.0f);
  EXPECT_EQ(m.m13, 3.0f);
  EXPECT_EQ(m.m23, 4.0f);
  EXPECT_EQ(m.m33, 5.0f);

  const float3x3& c = m;
  EXPECT_EQ(&c[1], &m.v[1]);
}
