#include <gtest/gtest.h>

import std;

import double3;
import float3;
import double4;
import double3x3;
import double4x4;
import simd_quatd;
import randomnumbers;

TEST(double3, Test_Perpendicular)
{
  std::random_device rd;
  RandomNumber random(rd());

  for (std::size_t i = 0; i < 100000; ++i)
  {
    double3 u{random.uniform(), random.uniform(), random.uniform()};
    double3 v{random.uniform(), random.uniform(), random.uniform()};

    double3 w = double3::perpendicular(u, v);

    EXPECT_TRUE(double3::is_perpendicular(w, u) && double3::is_perpendicular(w, v))
        << std::format("{} {} {}     {} {} {}\n", u.x, u.y, u.z, v.x, v.y, v.z);
  }
}

TEST(double3, Test_Perpendicular_small_norm_and_zero_vectors)
{
  // orthogonal vectors of small length are not collinear; the result is perpendicular to both
  double3 u{0.0, 0.0, 1e-6};
  double3 v{1e-6, 1e-6, 0.0};
  double3 w = double3::perpendicular(u, v);
  EXPECT_NEAR(w.length(), 1.0, 1e-12);
  EXPECT_NEAR(double3::dot(w, u.normalized()), 0.0, 1e-12);
  EXPECT_NEAR(double3::dot(w, v.normalized()), 0.0, 1e-12);

  // a zero vector imposes no constraint; the result is a unit vector perpendicular to the other one
  for (std::pair<double3, double3> pair : std::vector<std::pair<double3, double3>>{
           {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {{0.0, 0.0, 2.0}, {0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}})
  {
    w = double3::perpendicular(pair.first, pair.second);
    EXPECT_NEAR(w.length(), 1.0, 1e-12);
    EXPECT_NEAR(double3::dot(w, pair.first), 0.0, 1e-12);
    EXPECT_NEAR(double3::dot(w, pair.second), 0.0, 1e-12);
  }

  // parallel vectors of different lengths
  w = double3::perpendicular(double3(1e-6, 2e-6, 3e-6), double3(2.0, 4.0, 6.0));
  EXPECT_NEAR(w.length(), 1.0, 1e-12);
  EXPECT_NEAR(double3::dot(w, double3(1.0, 2.0, 3.0)), 0.0, 1e-12);
}

TEST(double3, Test_KeyHash_and_KeyEqual_are_consistent)
{
  std::unordered_set<double3, double3::KeyHash, double3::KeyEqual> keys;
  keys.insert(double3(1.0, 2.0, 3.0));
  keys.insert(double3(1.0, 2.0, 3.0 + 1e-4));
  keys.insert(double3(1.0, 2.0, 3.0));
  EXPECT_EQ(keys.size(), 2uz);

  // keys that compare equal must have equal hashes
  double3 a(1.0, 2.0, 3.0);
  double3 b(1.0, 2.0, 3.0 + 1e-4);
  EXPECT_TRUE(double3::KeyEqual()(a, double3(1.0, 2.0, 3.0)));
  EXPECT_EQ(double3::KeyHash()(a), double3::KeyHash()(double3(1.0, 2.0, 3.0)));
  EXPECT_FALSE(double3::KeyEqual()(a, b));

  EXPECT_TRUE(keys.contains(double3(1.0, 2.0, 3.0)));
  EXPECT_TRUE(keys.contains(double3(1.0, 2.0, 3.0 + 1e-4)));
  EXPECT_FALSE(keys.contains(double3(1.0, 2.0, 3.0 + 5e-4)));
}

TEST(double3, Test_fract_stays_below_one)
{
  // -5e-17 + 1.0 rounds to 1.0, which wraps to 0.0
  for (double value : {-5e-17, -1e-300, 5e-17, -0.25, 0.75, 3.0, -2.0, 0.999999})
  {
    double3 v(value, value, value);
    for (double3 f : {v.fract(), double3::fract(v)})
    {
      EXPECT_GE(f.x, 0.0) << value;
      EXPECT_LT(f.x, 1.0) << value;
      EXPECT_GE(f.y, 0.0) << value;
      EXPECT_LT(f.y, 1.0) << value;
      EXPECT_GE(f.z, 0.0) << value;
      EXPECT_LT(f.z, 1.0) << value;
    }
  }
  EXPECT_EQ(double3(-5e-17, 0.0, 0.0).fract().x, 0.0);
  EXPECT_EQ(double3::fract(double3(-5e-17, 0.0, 0.0)).x, 0.0);
  EXPECT_EQ(double3(-0.25, 0.0, 0.0).fract().x, 0.75);
}

TEST(float3, Test_fract_stays_below_one)
{
  for (float value : {-1e-9f, -0.25f, 0.75f, 3.0f})
  {
    float3 f = float3(value, value, value).fract();
    EXPECT_GE(f.x, 0.0f) << value;
    EXPECT_LT(f.x, 1.0f) << value;
    EXPECT_GE(f.z, 0.0f) << value;
    EXPECT_LT(f.z, 1.0f) << value;
  }
  EXPECT_EQ(float3(-1e-9f, 0.0f, 0.0f).fract().x, 0.0f);
}
