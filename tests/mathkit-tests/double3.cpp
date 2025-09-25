#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <numbers>
#include <print>
#include <random>

import double3;
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
