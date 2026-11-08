#include <gtest/gtest.h>

import std;

import float2;
import double2;
import float4;
import float4x4;
import double4x4;

TEST(free_operators, Test_float4_subtraction_and_double_scaling)
{
  float4 a(5.0f, 4.0f, 3.0f, 2.0f);
  float4 b(1.0f, 1.0f, 1.0f, 1.0f);

  float4 difference = a - b;
  EXPECT_EQ(difference.x, 4.0f);
  EXPECT_EQ(difference.w, 1.0f);

  float4 scaled = 2.0 * b;
  EXPECT_EQ(scaled.x, 2.0f);
  EXPECT_EQ(scaled.w, 2.0f);
}

TEST(free_operators, Test_stream_operators)
{
  std::ostringstream stream;
  stream << float2(1.5f, 2.0f) << " " << double2(-3.0, 4.25);
  EXPECT_EQ(stream.str(), "1.5,2 -3,4.25");

  std::istringstream input(stream.str());
  float2 f;
  double2 d;
  input >> f >> d;
  EXPECT_EQ(f.x, 1.5f);
  EXPECT_EQ(f.y, 2.0f);
  EXPECT_EQ(d.x, -3.0);
  EXPECT_EQ(d.y, 4.25);

  std::ostringstream matrices;
  matrices << double4x4::Identity() << " " << float4x4::Identity();
  EXPECT_EQ(matrices.str(), "1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1");
}
