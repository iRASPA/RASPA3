#include <gtest/gtest.h>

#include <cmath>

import double3;
import double4;
import double3x3;
import double4x4;
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

TEST(double3x3, Test_quaternion)
{
  RandomNumber random(10);

  std::vector<double3> positions{double3{0.0, 0.0, 0.0}, double3{0.1,0.2,0.3}, double3{0.5, 0.5, 0.5}, double3{1.23,0.3,0.34}, double3{-0.2, 0.3, -0.1}};
  double3 com_position{10,11,12};
  std::vector<double> mass{1.0, 0.1, 4.5, 0.2, 1.5};
  double3 origin = positions[2];
  double3 neworigin = positions[1];

  double3 com = (mass[0] * positions[0] +
                 mass[1] * positions[1] +
                 mass[2] * positions[2] +
                 mass[3] * positions[3] +
                 mass[4] * positions[4]) / (1.0 + 0.1 + 4.5 + 0.2 + 1.5);

  std::cout << "center of mass: " << com.x << ", " << com.y << ", " << com.z << std::endl;


  // https://faculty.sites.iastate.edu/jia/files/inline-files/homogeneous-transform.pdf
    simd_quatd q = random.randomSimdQuatd();
    double3x3 rotationMatrix = double3x3::buildRotationMatrixInverse(q);
  {

    std::cout << std::endl;

    for(const double3 &pos : positions)
    {
      double3 pnew = rotationMatrix * (pos-origin) + origin+com_position;
      std::cout << "Correct generated positions " << pnew.x << ", " << pnew.y << ", " << pnew.z << std::endl;
    }
    std::cout << std::endl;

    double4x4 translation = double4x4(1.0, 0.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0, 0.0,
                                      0.0, 0.0, 1.0, 0.0,
                                      origin.x, origin.y, origin.z, 1.0);
    double4x4 m = translation * double4x4(rotationMatrix) * double4x4::inverse(translation);

    for(const double3 &pos : positions)
    {
      double4 pnew2 = m * double4(pos.x, pos.y, pos.z, 1.0);
      std::cout << "Matrix: " << pnew2.x << ", " << pnew2.y << ", " << pnew2.z << std::endl;
    }
    std::cout << std::endl;


    double3x3 newm = m.toDouble3x3();
    double3 trans = double3(m.m41, m.m42, m.m43);

    for(const double3 &pos : positions)
    {
      double3 pnew2 = rotationMatrix * (pos-neworigin) + neworigin
         + rotationMatrix * (neworigin - origin) - (neworigin - origin);
      std::cout << pnew2.x << ", " << pnew2.y << ", " << pnew2.z << std::endl;
    }
    std::cout << std::endl;
  }
  {
    double4x4 translation1 = double4x4(1.0, 0.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0, 0.0,
                                        0.0, 0.0, 1.0, 0.0,
                                        com.x, com.y, com.z, 1.0);
    double4x4 m1 = translation1 * double4x4(rotationMatrix) * double4x4::inverse(translation1);
    double3x3 newm1 = m1.toDouble3x3();

    double3 t = rotationMatrix * (com - origin) - (com - origin) + com_position;
    for(const double3 &pos : positions)
    {
      double3 pnew2 = rotationMatrix * (pos - com) + com + t;
      std::cout << "com: " << pnew2.x << ", " << pnew2.y << ", " << pnew2.z << std::endl;
    }
    std::cout << std::endl;

  }
}
