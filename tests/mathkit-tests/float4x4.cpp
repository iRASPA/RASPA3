#include <gtest/gtest.h>

import std;

import float4;
import float4x4;

TEST(float4x4, Test_not_equal)
{
  // operator!= compares with another float4x4, not with a double4x4
  static_assert(std::is_same_v<decltype(&float4x4::operator!=), bool (float4x4::*)(const float4x4&) const>);

  float4x4 a(float4(1.0f, 0.0f, 0.0f, 0.0f), float4(0.0f, 1.0f, 0.0f, 0.0f), float4(0.0f, 0.0f, 1.0f, 0.0f),
             float4(0.0f, 0.0f, 0.0f, 1.0f));
  float4x4 b = a;
  EXPECT_FALSE(a != b);

  b.m[6] = 0.5f;
  EXPECT_TRUE(a != b);
}
