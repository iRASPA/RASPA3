module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#endif

export module float4x4;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import float3;
import float4;
import double3x3;
import double4x4;
import simd_quatd;
import float3x3;

export union float4x4
{
  float m[16];
  float mm[4][4];
  float4 v[4];
  struct
  {
    float m11, m21, m31, m41,  // 1st column
        m12, m22, m32, m42,    // 2nd column
        m13, m23, m33, m43,    // 3rd column
        m14, m24, m34, m44;    // 4th column, {tx,ty,tz,1}
  };
  struct
  {
    float ax, ay, az, aw,  // 1st column
        bx, by, bz, bw,    // 2nd column
        cx, cy, cz, cw,    // 3rd column
        wx, wy, wz, ww;    // 4th column, {tx,ty,tz,1}
  };

  inline float4& operator[](int i) { return v[i]; }
  inline const float4& operator[](int i) const { return v[i]; }

  float4x4(float m11 = 0.0, float m21 = 0.0, float m31 = 0.0, float m41 = 0.0, float m12 = 0.0, float m22 = 0.0,
           float m32 = 0.0, float m42 = 0.0, float m13 = 0.0, float m23 = 0.0, float m33 = 0.0, float m43 = 0.0,
           float m14 = 0.0, float m24 = 0.0, float m34 = 0.0, float m44 = 0.0)
      : m11(m11),
        m21(m21),
        m31(m31),
        m41(m41),
        m12(m12),
        m22(m22),
        m32(m32),
        m42(m42),
        m13(m13),
        m23(m23),
        m33(m33),
        m43(m43),
        m14(m14),
        m24(m24),
        m34(m34),
        m44(m44) {

        };

  float4x4(const float4x4&);             // copy constructor
  float4x4& operator=(const float4x4&);  // assignment operator
  float4x4(float4 v1, float4 v2, float4 v3, float4 v4) : v{v1, v2, v3, v4} {}
  float4x4(const float3x3&);
  float4x4(const double3x3&);
  float4x4(const double4x4&);
  float4x4(const float* v);

  inline bool operator==(const float4x4& b) const
  {
    bool equal = true;
    const float epsilon = 1e-8f;

    for (int i = 0; i < 16 && equal; i++) equal = std::fabs(this->m[i] - b.m[i]) <= epsilon;

    return equal;
  };
  inline bool operator!=(const double4x4& b) const { return !(*this == b); }

  float4x4 const transpose(void);
  static float4x4 const inverse(const float4x4& right);
  void inverse();
  static float4x4 TransformationAroundArbitraryPoint(float4x4 m, float3 p);

  static float4x4 Identity(void);
  static float4x4 AffinityMatrixToTransformationAroundArbitraryPoint(float4x4 m, float3 p);
  static float4x4 Rotationfloat4x4FromQuaternion(simd_quatd q);
};

export inline float4x4 operator*(const float4x4& a, const float4x4& b)
{
  float4x4 r;

  r.m11 = a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31 + a.m14 * b.m41;
  r.m21 = a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31 + a.m24 * b.m41;
  r.m31 = a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31 + a.m34 * b.m41;
  r.m41 = a.m41 * b.m11 + a.m42 * b.m21 + a.m43 * b.m31 + a.m44 * b.m41;

  r.m12 = a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32 + a.m14 * b.m42;
  r.m22 = a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32 + a.m24 * b.m42;
  r.m32 = a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32 + a.m34 * b.m42;
  r.m42 = a.m41 * b.m12 + a.m42 * b.m22 + a.m43 * b.m32 + a.m44 * b.m42;

  r.m13 = a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33 + a.m14 * b.m43;
  r.m23 = a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33 + a.m24 * b.m43;
  r.m33 = a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33 + a.m34 * b.m43;
  r.m43 = a.m41 * b.m13 + a.m42 * b.m23 + a.m43 * b.m33 + a.m44 * b.m43;

  r.m14 = a.m11 * b.m14 + a.m12 * b.m24 + a.m13 * b.m34 + a.m14 * b.m44;
  r.m24 = a.m21 * b.m14 + a.m22 * b.m24 + a.m23 * b.m34 + a.m24 * b.m44;
  r.m34 = a.m31 * b.m14 + a.m32 * b.m24 + a.m33 * b.m34 + a.m34 * b.m44;
  r.m44 = a.m41 * b.m14 + a.m42 * b.m24 + a.m43 * b.m34 + a.m44 * b.m44;

  return r;
}

export inline float4 operator*(const float4x4& a, const float4& b)
{
  float4 r;

  r.x = a.m11 * b.x + a.m12 * b.y + a.m13 * b.z + a.m14 * b.w;
  r.y = a.m21 * b.x + a.m22 * b.y + a.m23 * b.z + a.m24 * b.w;
  r.z = a.m31 * b.x + a.m32 * b.y + a.m33 * b.z + a.m34 * b.w;
  r.w = a.m41 * b.x + a.m42 * b.y + a.m43 * b.z + a.m44 * b.w;

  return r;
}
