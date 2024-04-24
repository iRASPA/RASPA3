#include <gtest/gtest.h>

#include <cmath>

import double3;
import double3x3;
import simd_quatd;
import randomnumbers;

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

TEST(double3x3, Test_rotation_matrix_from_quaternion)
{
  RandomNumber random(10);

  for(size_t i = 0; i < 10000; ++i)
  {
    simd_quatd q = random.randomSimdQuatd();
    double3x3 rotationMatrix = double3x3::buildRotationMatrix(q);
    simd_quatd p = rotationMatrix.quaternion();
    EXPECT_NEAR(std::abs(p.ix), std::abs(q.ix), 1e-6) << "wrong ix";
    EXPECT_NEAR(std::abs(p.iy), std::abs(q.iy), 1e-6) << "wrong iy";
    EXPECT_NEAR(std::abs(p.iz), std::abs(q.iz), 1e-6) << "wrong iz";
    EXPECT_NEAR(std::abs(p.r),  std::abs(q.r),  1e-6) << "wrong r";
  }
}
