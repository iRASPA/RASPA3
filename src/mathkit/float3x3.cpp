module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <numbers>
#include <vector>
#endif

module float3x3;

#ifndef USE_LEGACY_HEADERS
import <numbers>;
import <vector>;
import <cmath>;
#endif

import simd_quatd;
import float3;

#define sqr(x) ((x) * (x))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// dlambda_limit, below which two lambdas are relatively equal
float flambda_limit = 1.0E-3f;
float fiszero_limit = 1.0E-20f;

float3x3::float3x3(simd_quatd q)
{
  double sqw = q.r * q.r;
  double sqx = q.ix * q.ix;
  double sqy = q.iy * q.iy;
  double sqz = q.iz * q.iz;

  // invs (inverse square length) is only required if quaternion is not already normalised
  double invs = 1 / (sqx + sqy + sqz + sqw);
  m11 = static_cast<float>((sqx - sqy - sqz + sqw) * invs);  // since sqw + sqx + sqy + sqz =1/invs*invs
  m22 = static_cast<float>((-sqx + sqy - sqz + sqw) * invs);
  m33 = static_cast<float>((-sqx - sqy + sqz + sqw) * invs);

  double tmp1 = q.ix * q.iy;
  double tmp2 = q.iz * q.r;
  m21 = static_cast<float>(2.0 * (tmp1 + tmp2) * invs);
  m12 = static_cast<float>(2.0 * (tmp1 - tmp2) * invs);

  tmp1 = q.ix * q.iz;
  tmp2 = q.iy * q.r;
  m31 = static_cast<float>(2.0 * (tmp1 - tmp2) * invs);
  m13 = static_cast<float>(2.0 * (tmp1 + tmp2) * invs);
  tmp1 = q.iy * q.iz;
  tmp2 = q.ix * q.r;
  m32 = static_cast<float>(2.0 * (tmp1 + tmp2) * invs);
  m23 = static_cast<float>(2.0 * (tmp1 - tmp2) * invs);
}

float float3x3::determinant(void)
{
  float determinant = +m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  return determinant;
}

float float3x3::trace(void) { return m11 + m22 + m33; }

float3x3 const float3x3::inverse()
{
  float determinant = +m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  float3x3 inverse;
  inverse.m11 = +(m22 * m33 - m32 * m23) / determinant;
  inverse.m21 = -(m21 * m33 - m31 * m23) / determinant;
  inverse.m31 = +(m21 * m32 - m31 * m22) / determinant;
  inverse.m12 = -(m12 * m33 - m32 * m13) / determinant;
  inverse.m22 = +(m11 * m33 - m31 * m13) / determinant;
  inverse.m32 = -(m11 * m32 - m31 * m12) / determinant;
  inverse.m13 = +(m12 * m23 - m22 * m13) / determinant;
  inverse.m23 = -(m11 * m23 - m21 * m13) / determinant;
  inverse.m33 = +(m11 * m22 - m21 * m12) / determinant;

  return inverse;
}

float3x3 const float3x3::inverse(const float3x3& a)
{
  float determinant = +a.m11 * (a.m22 * a.m33 - a.m23 * a.m32) - a.m12 * (a.m21 * a.m33 - a.m23 * a.m31) +
                      a.m13 * (a.m21 * a.m32 - a.m22 * a.m31);

  float3x3 inverse;
  inverse.m11 = +(a.m22 * a.m33 - a.m32 * a.m23) / determinant;
  inverse.m21 = -(a.m21 * a.m33 - a.m31 * a.m23) / determinant;
  inverse.m31 = +(a.m21 * a.m32 - a.m31 * a.m22) / determinant;
  inverse.m12 = -(a.m12 * a.m33 - a.m32 * a.m13) / determinant;
  inverse.m22 = +(a.m11 * a.m33 - a.m31 * a.m13) / determinant;
  inverse.m32 = -(a.m11 * a.m32 - a.m31 * a.m12) / determinant;
  inverse.m13 = +(a.m12 * a.m23 - a.m22 * a.m13) / determinant;
  inverse.m23 = -(a.m11 * a.m23 - a.m21 * a.m13) / determinant;
  inverse.m33 = +(a.m11 * a.m22 - a.m21 * a.m12) / determinant;

  return inverse;
}

float3x3 const float3x3::transpose(void)
{
  float3x3 res;

  res.m11 = m11;
  res.m12 = m21;
  res.m13 = m31;
  res.m21 = m12;
  res.m22 = m22;
  res.m23 = m32;
  res.m31 = m13;
  res.m32 = m23;
  res.m33 = m33;

  return res;
}

float3x3 float3x3::inversetranpose()
{
  float determinant = +m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  float3x3 inverse;
  inverse.m11 = +(m22 * m33 - m32 * m23) / determinant;
  inverse.m12 = -(m21 * m33 - m31 * m23) / determinant;
  inverse.m13 = +(m21 * m32 - m31 * m22) / determinant;
  inverse.m21 = -(m12 * m33 - m32 * m13) / determinant;
  inverse.m22 = +(m11 * m33 - m31 * m13) / determinant;
  inverse.m23 = -(m11 * m32 - m31 * m12) / determinant;
  inverse.m31 = +(m12 * m23 - m22 * m13) / determinant;
  inverse.m32 = -(m11 * m23 - m21 * m13) / determinant;
  inverse.m33 = +(m11 * m22 - m21 * m12) / determinant;

  return inverse;
}

float trunc_sqrt(float x) { return (x <= 0.0f ? 0.0f : std::sqrt(x)); }

float trunc_acos(float x)
{
  if (x >= 1.0f) return 0.0f;
  if (x <= -1.0f) return float(std::numbers::pi);
  return acos(x);
}

static double sign(double x) { return (x < 0.0 ? -1.0 : 1.0); }
