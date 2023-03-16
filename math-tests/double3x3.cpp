#include <gtest/gtest.h>

import double3;
import double3x3;

TEST(double3x3, Test_multiplication)
{
  double3x3 m{ double3{25.0, 2.5e-6, 0.0}, double3{0.0, 25.0, 0.0}, double3{0.0, 0.0, 25.0} };
  double3x3 n{ double3{0.04, -4.0e-9, 0.0}, double3{0.0, 0.04, 0.0}, double3{0.0, 0.0, 0.04} };
  
  double3x3 identity = m * n;
  EXPECT_NEAR(identity.ax, 1.0, 1e-10) << "wrong ax";
  EXPECT_NEAR(identity.bx, 0.0, 1e-10) << "wrong bx";
  EXPECT_NEAR(identity.cx, 0.0, 1e-10) << "wrong cx";
  EXPECT_NEAR(identity.ay, 0.0, 1e-10) << "wrong ay";
  EXPECT_NEAR(identity.by, 1.0, 1e-10) << "wrong by";
  EXPECT_NEAR(identity.cy, 0.0, 1e-10) << "wrong cy";
  EXPECT_NEAR(identity.az, 0.0, 1e-10) << "wrong az";
  EXPECT_NEAR(identity.bz, 0.0, 1e-10) << "wrong bz";
  EXPECT_NEAR(identity.cz, 1.0, 1e-10) << "wrong cz";
}

TEST(double3x3, Test_inverse)
{
  double3x3 m{double3{25.0, 2.5e-6, 0.0}, double3{0.0, 25.0, 0.0}, double3{0.0, 0.0, 25.0}};
  double3x3 inv = m.inverse();

  double3x3 identity = inv * m;

  EXPECT_NEAR(identity.ax, 1.0, 1e-10) << "wrong ax";
  EXPECT_NEAR(identity.bx, 0.0, 1e-10) << "wrong bx";
  EXPECT_NEAR(identity.cx, 0.0, 1e-10) << "wrong cx";
  EXPECT_NEAR(identity.bx, 0.0, 1e-10) << "wrong bx";
  EXPECT_NEAR(identity.by, 1.0, 1e-10) << "wrong by";
  EXPECT_NEAR(identity.bz, 0.0, 1e-10) << "wrong bz";
  EXPECT_NEAR(identity.cx, 0.0, 1e-10) << "wrong cx";
  EXPECT_NEAR(identity.cy, 0.0, 1e-10) << "wrong cy";
  EXPECT_NEAR(identity.cz, 1.0, 1e-10) << "wrong cz";
}