module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <utility>
#include <vector>
#pragma push_macro("__SSE3__")
#undef __SSE3__
#include <random>
#pragma pop_macro("__SSE3__")
#endif

#ifndef USE_LEGACY_HEADERS
#include <cstdlib>
#endif

module randomnumbers;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import simd_quatd;
import archive;

double3x3 RandomNumber::randomRotationMatrix()
{
  double X0 = uniform();
  double Y1 = 2.0 * std::numbers::pi * uniform();
  double Y2 = 2.0 * std::numbers::pi * uniform();
  double R1 = std::sqrt(1.0 - X0);
  double R2 = std::sqrt(X0);
  double U0 = std::cos(Y2) * R2;
  double U1 = std::sin(Y1) * R1;
  double U2 = std::cos(Y1) * R1;
  double U3 = std::sin(Y2) * R2;
  double COEFI = 2.0 * U0 * U0 - 1.0;
  double COEFUU = 2.0;
  double COEFE = 2.0 * U0;

  double3x3 R;
  R[0][0] = COEFI + COEFUU * U1 * U1;
  R[1][1] = COEFI + COEFUU * U2 * U2;
  R[2][2] = COEFI + COEFUU * U3 * U3;
  R[1][2] = COEFUU * U2 * U3 - COEFE * U1;
  R[2][0] = COEFUU * U3 * U1 - COEFE * U2;
  R[0][1] = COEFUU * U1 * U2 - COEFE * U3;
  R[2][1] = COEFUU * U3 * U2 + COEFE * U1;
  R[0][2] = COEFUU * U1 * U3 + COEFE * U2;
  R[1][0] = COEFUU * U2 * U1 + COEFE * U3;
  return R;
}

double3x3 RandomNumber::randomRotationAroundX(double angle)
{
  double c = std::cos(angle);
  double s = std::sin(angle);
  return double3x3(double3(1.0, 0.0, 0.0), double3(0.0, c, s), double3(0.0, -s, c));
}

double3x3 RandomNumber::randomRotationAroundY(double angle)
{
  double c = std::cos(angle);
  double s = std::sin(angle);
  return double3x3(double3(c, 0.0, -s), double3(0.0, 1.0, 0.0), double3(s, 0.0, c));
}

double3x3 RandomNumber::randomRotationAroundZ(double angle)
{
  double c = std::cos(angle);
  double s = std::sin(angle);
  return double3x3(double3(c, s, 0.0), double3(-s, c, 0.0), double3(0.0, 0.0, 1.0));
}

double3 RandomNumber::randomVectorOnUnitSphere()
{
  double ran1, ran2, ranh, ransq = 0.0;

  do
  {
    ran1 = 2.0 * uniform() - 1.0;
    ran2 = 2.0 * uniform() - 1.0;
    ransq = ran1 * ran1 + ran2 * ran2;
  } while (ransq >= 1.0);

  ranh = 2.0 * std::sqrt(1.0 - ransq);
  return double3(ran1 * ranh, ran2 * ranh, 1.0 - 2.0 * ransq);
}

// angle between 0..pi
double3 RandomNumber::randomVectorOnCone(double3 v, double angle)
{
  double3 normalized_v = v.normalized();

  double3 w = randomVectorOnUnitSphere();

  w -= double3::dot(normalized_v, w) * normalized_v;

  return (std::sin(angle) * w.normalized() + std::cos(angle) * normalized_v).normalized();
}

simd_quatd RandomNumber::randomSimdQuatd()
{
  double s = 0.0;
  double sigma1 = 0.0;
  double sigma2 = 0.0;
  double theta1 = 0.0;
  double theta2 = 0.0;

  s = uniform();
  sigma1 = std::sqrt(1.0 - s);
  sigma2 = std::sqrt(s);
  theta1 = 2.0 * std::numbers::pi * uniform();
  theta2 = 2.0 * std::numbers::pi * uniform();

  return simd_quatd(sigma2 * std::cos(theta2), sigma1 * std::sin(theta1), sigma1 * std::cos(theta1),
                    sigma2 * std::sin(theta2));
}

simd_quatd RandomNumber::smallRandomQuaternion(double angleRange)
{
  double3 randomDirection = randomVectorOnUnitSphere();
  double angle = angleRange * 2.0 * ((double(rand()) / RAND_MAX) - 0.5);
  return simd_quatd::fromAxisAngle(angle, randomDirection).normalized();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RandomNumber &r)
{
  archive << r.seed;
  archive << r.count;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RandomNumber &r)
{
  archive >> r.seed;
  archive >> r.count;
  r.mt = std::mt19937_64(r.seed);
  r.mt.discard(r.count);
  return archive;
}
