#include <gtest/gtest.h>

import std;

import json;

// JSON has no NaN or infinity: the serializer writes null for them. Its test was std::isfinite, which the fast-math
// flags of the build compile to true, and the dump held the invalid tokens nan and inf. (The values are made from
// their bit patterns: the compiler may assume that a computation gives no NaN or infinity.)
TEST(json, a_nan_or_an_infinity_is_written_as_null)
{
  const double nan = std::bit_cast<double>(0x7ff8000000000000ull);
  const double infinity = std::bit_cast<double>(0x7ff0000000000000ull);
  nlohmann::json values = nlohmann::json::array({1.5, nan, infinity});
  EXPECT_EQ(values.dump(), "[1.5,null,null]");
}
