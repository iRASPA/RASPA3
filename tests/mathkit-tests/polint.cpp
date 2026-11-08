#include <gtest/gtest.h>

import std;

import polint;

TEST(polint, Test_reproduces_a_polynomial_of_degree_N_minus_1)
{
  std::array<double, 5> xa{-2.0, -1.0, 0.0, 1.0, 2.0};
  std::array<double, 5> ya{};
  for (std::size_t i = 0; i < xa.size(); ++i)
  {
    double x = xa[i];
    ya[i] = 3.0 - 2.0 * x + 0.5 * x * x - x * x * x + 0.25 * x * x * x * x;
  }

  // nearest to the first abscissa, to the last one, and to one in the middle
  for (double x : {-2.3, -1.9, 0.4, 1.6, 2.2})
  {
    double y{};
    double dy{};
    Interpolation::polint<5>(xa, ya, x, &y, &dy);
    EXPECT_NEAR(y, 3.0 - 2.0 * x + 0.5 * x * x - x * x * x + 0.25 * x * x * x * x, 1e-12) << x;
  }
}

TEST(polint, Test_equal_abscissae_throw)
{
  std::array<double, 3> xa{0.0, 1.0, 1.0};
  std::array<double, 3> ya{1.0, 2.0, 3.0};
  double y{};
  double dy{};
  EXPECT_THROW(Interpolation::polint<3>(xa, ya, 0.5, &y, &dy), std::invalid_argument);
}
