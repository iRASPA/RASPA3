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

TEST(random, Test_RandomVectorOnCone)
{
  std::random_device rd;
  RandomNumber random(rd());

  for (std::size_t i = 0; i < 100000; ++i)
  {
    double3 v{random.uniform(), random.uniform(), random.uniform()};
    double imposed_angle = random.uniform() * std::numbers::pi;

    double3 w = random.randomVectorOnCone(v, imposed_angle);

    double computed_angle = double3::angle(v, w);

    EXPECT_NEAR(imposed_angle, computed_angle, 1e-10)
        << std::format("wrong angle: {} vs {}\n", imposed_angle, computed_angle);
  }
}
