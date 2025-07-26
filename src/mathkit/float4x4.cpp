module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#endif

module float4x4;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import float3;
import float4;
import float3x3;
import double3x3;
import double4x4;
import simd_quatd;

float4x4::float4x4(const float4x4& a)  // copy constructor
{
  m11 = a.m11;
  m21 = a.m21;
  m31 = a.m31;
  m41 = a.m41;
  m12 = a.m12;
  m22 = a.m22;
  m32 = a.m32;
  m42 = a.m42;
  m13 = a.m13;
  m23 = a.m23;
  m33 = a.m33;
  m43 = a.m43;
  m14 = a.m14;
  m24 = a.m24;
  m34 = a.m34;
  m44 = a.m44;
}

float4x4::float4x4(const float* v)
{
  for (int i = 0; i < 16; i++)
  {
    m[i] = v[i];
  }
}

float4x4& float4x4::operator=(const float4x4& a)  // assignment operator
{
  if (this != &a)  // protect against self-assingment
  {
    m11 = a.m11;
    m21 = a.m21;
    m31 = a.m31;
    m41 = a.m41;
    m12 = a.m12;
    m22 = a.m22;
    m32 = a.m32;
    m42 = a.m42;
    m13 = a.m13;
    m23 = a.m23;
    m33 = a.m33;
    m43 = a.m43;
    m14 = a.m14;
    m24 = a.m24;
    m34 = a.m34;
    m44 = a.m44;
  }
  return *this;
}

float4x4::float4x4(const float3x3& a)
{
  m11 = a.m11;
  m21 = a.m21;
  m31 = a.m31;
  m41 = 0.0;
  m12 = a.m12;
  m22 = a.m22;
  m32 = a.m32;
  m42 = 0.0;
  m13 = a.m13;
  m23 = a.m23;
  m33 = a.m33;
  m43 = 0.0;
  m14 = 0.0;
  m24 = 0.0;
  m34 = 0.0;
  m44 = 1.0;
}

float4x4::float4x4(const double3x3& a)
{
  m11 = float(a.m11);
  m21 = float(a.m21);
  m31 = float(a.m31);
  m41 = float(0.0);
  m12 = float(a.m12);
  m22 = float(a.m22);
  m32 = float(a.m32);
  m42 = float(0.0);
  m13 = float(a.m13);
  m23 = float(a.m23);
  m33 = float(a.m33);
  m43 = float(0.0);
  m14 = float(0.0);
  m24 = float(0.0);
  m34 = float(0.0);
  m44 = float(1.0);
}

float4x4::float4x4(const double4x4& a)
{
  m11 = float(a.m11);
  m21 = float(a.m21);
  m31 = float(a.m31);
  m41 = float(a.m41);
  m12 = float(a.m12);
  m22 = float(a.m22);
  m32 = float(a.m32);
  m42 = float(a.m42);
  m13 = float(a.m13);
  m23 = float(a.m23);
  m33 = float(a.m33);
  m43 = float(a.m43);
  m14 = float(a.m14);
  m24 = float(a.m24);
  m34 = float(a.m34);
  m44 = float(a.m44);
}

float4x4 const float4x4::transpose(void)
{
  float4x4 res;

  res.m11 = m11;
  res.m12 = m21;
  res.m13 = m31;
  res.m14 = m41;
  res.m21 = m12;
  res.m22 = m22;
  res.m23 = m32;
  res.m24 = m42;
  res.m31 = m13;
  res.m32 = m32;
  res.m33 = m33;
  res.m34 = m43;
  res.m41 = m14;
  res.m42 = m42;
  res.m43 = m34;
  res.m44 = m44;

  return res;
}

float4x4 const float4x4::inverse(const float4x4& a)
{
  float4x4 inverse;

  float d01 = (a.m13 * a.m24 - a.m14 * a.m23);
  float d02 = (a.m13 * a.m34 - a.m14 * a.m33);
  float d12 = (a.m23 * a.m34 - a.m24 * a.m33);
  float d13 = (a.m23 * a.m44 - a.m24 * a.m43);
  float d23 = (a.m33 * a.m44 - a.m34 * a.m43);
  float d30 = (a.m43 * a.m14 - a.m44 * a.m13);

  const float cell0 = (a.m22 * d23 - a.m32 * d13 + a.m42 * d12);
  const float cell4 = -(a.m12 * d23 + a.m32 * d30 + a.m42 * d02);
  const float cell8 = (a.m12 * d13 + a.m22 * d30 + a.m42 * d01);
  const float cell12 = -(a.m12 * d12 - a.m22 * d02 + a.m32 * d01);

  const float determinant = a.m11 * cell0 + a.m21 * cell4 + a.m31 * cell8 + a.m41 * cell12;
  if (std::fabs(determinant) < 1e-8f)
  {
    std::cout << "Error";
  }

  inverse.m11 = cell0;
  inverse.m12 = cell4;
  inverse.m13 = cell8;
  inverse.m14 = cell12;

  inverse.m21 = -(a.m21 * d23 - a.m31 * d13 + a.m41 * d12);
  inverse.m22 = (a.m11 * d23 + a.m31 * d30 + a.m41 * d02);
  inverse.m23 = -(a.m11 * d13 + a.m21 * d30 + a.m41 * d01);
  inverse.m24 = (a.m11 * d12 - a.m21 * d02 + a.m31 * d01);

  d01 = a.m11 * a.m22 - a.m12 * a.m21;
  d02 = a.m11 * a.m32 - a.m12 * a.m31;
  d12 = a.m21 * a.m32 - a.m22 * a.m31;
  d13 = a.m21 * a.m42 - a.m22 * a.m41;
  d23 = a.m31 * a.m42 - a.m32 * a.m41;
  d30 = a.m41 * a.m12 - a.m42 * a.m11;

  inverse.m31 = (a.m24 * d23 - a.m34 * d13 + a.m44 * d12);
  inverse.m32 = -(a.m14 * d23 + a.m34 * d30 + a.m44 * d02);
  inverse.m33 = (a.m14 * d13 + a.m24 * d30 + a.m44 * d01);
  inverse.m34 = -(a.m14 * d12 - a.m24 * d02 + a.m34 * d01);

  inverse.m41 = -(a.m23 * d23 - a.m33 * d13 + a.m43 * d12);
  inverse.m42 = (a.m13 * d23 + a.m33 * d30 + a.m43 * d02);
  inverse.m43 = -(a.m13 * d13 + a.m23 * d30 + a.m43 * d01);
  inverse.m44 = (a.m13 * d12 - a.m23 * d02 + a.m33 * d01);

  const float invDeterminant = 1.0f / determinant;
  inverse.m11 *= invDeterminant;
  inverse.m21 *= invDeterminant;
  inverse.m31 *= invDeterminant;
  inverse.m41 *= invDeterminant;
  inverse.m12 *= invDeterminant;
  inverse.m22 *= invDeterminant;
  inverse.m32 *= invDeterminant;
  inverse.m42 *= invDeterminant;
  inverse.m13 *= invDeterminant;
  inverse.m23 *= invDeterminant;
  inverse.m33 *= invDeterminant;
  inverse.m43 *= invDeterminant;
  inverse.m14 *= invDeterminant;
  inverse.m24 *= invDeterminant;
  inverse.m34 *= invDeterminant;
  inverse.m44 *= invDeterminant;
  return inverse;
}

void float4x4::inverse()
{
  float4x4 inverse;

  float d01 = (m13 * m24 - m14 * m23);
  float d02 = (m13 * m34 - m14 * m33);
  float d12 = (m23 * m34 - m24 * m33);
  float d13 = (m23 * m44 - m24 * m43);
  float d23 = (m33 * m44 - m34 * m43);
  float d30 = (m43 * m14 - m44 * m13);

  const float cell0 = (m22 * d23 - m32 * d13 + m42 * d12);
  const float cell4 = -(m12 * d23 + m32 * d30 + m42 * d02);
  const float cell8 = (m12 * d13 + m22 * d30 + m42 * d01);
  const float cell12 = -(m12 * d12 - m22 * d02 + m32 * d01);

  const float determinant = m11 * cell0 + m21 * cell4 + m31 * cell8 + m41 * cell12;
  if (std::fabs(determinant) < 1e-8f)
  {
    std::cout << "Error";
  }

  inverse.m11 = cell0;
  inverse.m12 = cell4;
  inverse.m13 = cell8;
  inverse.m14 = cell12;

  inverse.m21 = -(m21 * d23 - m31 * d13 + m41 * d12);
  inverse.m22 = (m11 * d23 + m31 * d30 + m41 * d02);
  inverse.m23 = -(m11 * d13 + m21 * d30 + m41 * d01);
  inverse.m24 = (m11 * d12 - m21 * d02 + m31 * d01);

  d01 = m11 * m22 - m12 * m21;
  d02 = m11 * m32 - m12 * m31;
  d12 = m21 * m32 - m22 * m31;
  d13 = m21 * m42 - m22 * m41;
  d23 = m31 * m42 - m32 * m41;
  d30 = m41 * m12 - m42 * m11;

  inverse.m31 = (m24 * d23 - m34 * d13 + m44 * d12);
  inverse.m32 = -(m14 * d23 + m34 * d30 + m44 * d02);
  inverse.m33 = (m14 * d13 + m24 * d30 + m44 * d01);
  inverse.m34 = -(m14 * d12 - m24 * d02 + m34 * d01);

  inverse.m41 = -(m23 * d23 - m33 * d13 + m43 * d12);
  inverse.m42 = (m13 * d23 + m33 * d30 + m43 * d02);
  inverse.m43 = -(m13 * d13 + m23 * d30 + m43 * d01);
  inverse.m44 = (m13 * d12 - m23 * d02 + m33 * d01);

  const float invDeterminant = 1.0f / determinant;
  inverse.m11 *= invDeterminant;
  inverse.m21 *= invDeterminant;
  inverse.m31 *= invDeterminant;
  inverse.m41 *= invDeterminant;
  inverse.m12 *= invDeterminant;
  inverse.m22 *= invDeterminant;
  inverse.m32 *= invDeterminant;
  inverse.m42 *= invDeterminant;
  inverse.m13 *= invDeterminant;
  inverse.m23 *= invDeterminant;
  inverse.m33 *= invDeterminant;
  inverse.m43 *= invDeterminant;
  inverse.m14 *= invDeterminant;
  inverse.m24 *= invDeterminant;
  inverse.m34 *= invDeterminant;
  inverse.m44 *= invDeterminant;
  *this = inverse;
}

float4x4 float4x4::TransformationAroundArbitraryPoint(float4x4 m, float3 p)
{
  float3x3 R;

  R.m11 = 1.0f - m.m11;
  R.m12 = -m.m12;
  R.m13 = -m.m13;
  R.m21 = -m.m21;
  R.m22 = 1.0f - m.m22;
  R.m23 = -m.m23;
  R.m31 = -m.m31;
  R.m32 = -m.m32;
  R.m33 = 1.0f - m.m33;

  m.m14 += R.m11 * p.x + R.m12 * p.y + R.m13 * p.z;
  m.m24 += R.m21 * p.x + R.m22 * p.y + R.m23 * p.z;
  m.m34 += R.m31 * p.x + R.m32 * p.y + R.m33 * p.z;

  return m;
}

float4x4 float4x4::Identity(void)
{
  float4x4 R;

  R.m11 = 1.0;
  R.m12 = 0.0;
  R.m13 = 0.0;
  R.m14 = 0.0;
  R.m21 = 0.0;
  R.m22 = 1.0;
  R.m23 = 0.0;
  R.m24 = 0.0;
  R.m31 = 0.0;
  R.m32 = 0.0;
  R.m33 = 1.0;
  R.m34 = 0.0;
  R.m41 = 0.0;
  R.m42 = 0.0;
  R.m43 = 0.0;
  R.m44 = 1.0;

  return R;
}

// Book: "Essential Mathematics for Games & Interactive Applications", 2nd edition, page 157
float4x4 float4x4::AffinityMatrixToTransformationAroundArbitraryPoint(float4x4 m, float3 p)
{
  float4x4 R;

  R.m11 = 1.0f - m.m11;
  R.m12 = -m.m12;
  R.m13 = -m.m13;
  R.m21 = -m.m21;
  R.m22 = 1.0f - m.m22;
  R.m23 = -m.m23;
  R.m31 = -m.m31;
  R.m32 = -m.m32;
  R.m33 = 1.0f - m.m33;

  m.m14 += R.m11 * p.x + R.m12 * p.y + R.m13 * p.z;
  m.m24 += R.m21 * p.x + R.m22 * p.y + R.m23 * p.z;
  m.m34 += R.m31 * p.x + R.m32 * p.y + R.m33 * p.z;

  return m;
}

float4x4 float4x4::Rotationfloat4x4FromQuaternion(simd_quatd q)
{
  float4x4 m;
  double sqw = q.r * q.r;
  double sqx = q.ix * q.ix;
  double sqy = q.iy * q.iy;
  double sqz = q.iz * q.iz;

  // invs (inverse square length) is only required if quaternion is not already normalised
  double invs = 1 / (sqx + sqy + sqz + sqw);
  m.m11 = float((sqx - sqy - sqz + sqw) * invs);  // since sqw + sqx + sqy + sqz =1/invs*invs
  m.m22 = float((-sqx + sqy - sqz + sqw) * invs);
  m.m33 = float((-sqx - sqy + sqz + sqw) * invs);
  m.m44 = 1.0;

  double tmp1 = q.ix * q.iy;
  double tmp2 = q.iz * q.r;
  m.m21 = float(2.0 * (tmp1 + tmp2) * invs);
  m.m12 = float(2.0 * (tmp1 - tmp2) * invs);
  m.m14 = m.m41 = 0.0;

  tmp1 = q.ix * q.iz;
  tmp2 = q.iy * q.r;
  m.m31 = float(2.0 * (tmp1 - tmp2) * invs);
  m.m13 = float(2.0 * (tmp1 + tmp2) * invs);
  m.m24 = m.m42 = 0.0;

  tmp1 = q.iy * q.iz;
  tmp2 = q.ix * q.r;
  m.m32 = float(2.0 * (tmp1 + tmp2) * invs);
  m.m23 = float(2.0 * (tmp1 - tmp2) * invs);
  m.m34 = m.m43 = 0.0;

  return m;
}

std::ostream& operator<<(std::ostream& stream, const float4x4& v)  // output
{
  stream << v.m11 << v.m21 << v.m31 << v.m41;
  stream << v.m12 << v.m22 << v.m32 << v.m42;
  stream << v.m13 << v.m23 << v.m33 << v.m43;
  stream << v.m14 << v.m24 << v.m34 << v.m44;
  return stream;
}
