#include <gtest/gtest.h>

import std;

import special_functions;

// For z < -1 Pfaff's z / (z - 1) is near 1, and its distance from 1 was formed as 1 - z / (z - 1): it lost the digits
// of |z| and was exactly 1 for z below -2^53, where the function came out infinite (for b < a, where Gauss's value at
// 1 of the transformed function is). For b < a it decays as Gamma(c) Gamma(a - b) / (Gamma(a) Gamma(c - b)) (-z)^-b,
// up to a relative 1 / |z|.
TEST(special_functions, hypergeometric2F1_decays_for_large_negative_z)
{
  const double a = 1.5, b = 0.5, c = 2.0;
  const double leading = std::tgamma(c) * std::tgamma(a - b) / (std::tgamma(a) * std::tgamma(c - b));
  for (const double z : {-1.0e6, -1.0e12, -1.0e17, -1.0e20})
  {
    const double expected = leading * std::pow(-z, -b);
    const double value = hypergeometric2F1(a, b, c, z);
    EXPECT_EQ(std::to_string(value).find("inf"), std::string::npos) << z;
    EXPECT_NEAR(value, expected, 1.0e-5 * expected) << z;
  }
}

// A terminating series ends before a pole of c at a negative integer below a: 2F1(-1, 2; -2; z) = 1 + z. For
// |z| <= 1/2 the series was taken first and divided 0 by 0 at that pole.
TEST(special_functions, hypergeometric2F1_terminates_before_a_pole_of_c)
{
  for (const double z : {0.25, -0.4, 0.75, -3.0})
  {
    // (a NaN read off the printed value: the fast-math build folds the comparison of EXPECT_NEAR with one to a pass)
    const double value = hypergeometric2F1(-1.0, 2.0, -2.0, z);
    EXPECT_EQ(std::to_string(value).find("nan"), std::string::npos) << z;
    EXPECT_NEAR(value, 1.0 + z, 1.0e-14) << z;
  }
}
