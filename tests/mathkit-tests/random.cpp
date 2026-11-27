#include <gtest/gtest.h>

import std;

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

// An infinite weight (an overflowed Boltzmann factor) outweighs the finite ones: the draw used to be infinite, past
// every running sum, and took the last category. Finite weights whose sum overflows are drawn in proportion, and a
// NaN weight is refused (the fast-math flags of the build fold the comparisons that were to catch it).
TEST(random, categorical_distribution_handles_infinite_overflowing_and_nan_weights)
{
  RandomNumber random(7);
  // from the bit patterns: the fast-math flags of the build let the compiler assume that no infinity or NaN occurs,
  // so std::numeric_limits does not reliably give one
  const double infinity = std::bit_cast<double>(0x7ff0000000000000ull);
  const double nan = std::bit_cast<double>(0x7ff8000000000000ull);
  const std::vector<double> oneInfinite{infinity, 1.0};
  for (std::size_t draw = 0; draw != 1000; ++draw) EXPECT_EQ(random.categoricalDistribution(oneInfinite), 0uz);

  const std::vector<double> twoInfinite{1.0, infinity, 0.0, infinity};
  std::array<std::size_t, 4> counts{};
  for (std::size_t draw = 0; draw != 20000; ++draw) ++counts[random.categoricalDistribution(twoInfinite)];
  EXPECT_EQ(counts[0] + counts[2], 0uz);
  EXPECT_NEAR(static_cast<double>(counts[1]) / 20000.0, 0.5, 0.02);

  const std::vector<double> overflowing{1.5e308, 0.5e308, 1.5e308};
  counts = {};
  for (std::size_t draw = 0; draw != 35000; ++draw) ++counts[random.categoricalDistribution(overflowing)];
  EXPECT_NEAR(static_cast<double>(counts[1]) / 35000.0, 1.0 / 7.0, 0.01);
  EXPECT_NEAR(static_cast<double>(counts[2]) / 35000.0, 3.0 / 7.0, 0.015);

  const std::vector<double> withNaN{1.0, nan};
  EXPECT_THROW(static_cast<void>(random.categoricalDistribution(withNaN)), std::invalid_argument);
}
