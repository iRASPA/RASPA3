module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#endif

export module double4x4;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
#endif

import double3;
import double4;
import simd_quatd;
import double3x3;

export union double4x4
{
    double m[16];
    double mm[4][4];
    double4 v[4];
    struct {
        double m11, m21, m31, m41,   // 1st column
            m12, m22, m32, m42,   // 2nd column
            m13, m23, m33, m43,   // 3rd column
            m14, m24, m34, m44;   // 4th column, {tx,ty,tz,1}
    };
    struct {
        double ax, ay, az, aw,    // 1st column
            bx, by, bz, bw,    // 2nd column
            cx, cy, cz, cw,    // 3rd column
            wx, wy, wz, ww;    // 4th column, {tx,ty,tz,1}
    };

    inline double4& operator [] (int i) { return v[i]; }
    inline const double4& operator [] (int i) const { return v[i]; }

    double4x4()
      : m11(0.0), m21(0.0), m31(0.0), m41(0.0),
        m12(0.0), m22(0.0), m32(0.0), m42(0.0),
        m13(0.0), m23(0.0), m33(0.0), m43(0.0),
        m14(0.0), m24(0.0), m34(0.0), m44(0.0)
    {};
  
    explicit double4x4(double m11, double m22, double m33, double m44)
        : m11(m11), m21(0.0), m31(0.0), m41(0.0),
          m12(0.0), m22(m22), m32(0.0), m42(0.0),
          m13(0.0), m23(0.0), m33(m33), m43(0.0),
          m14(0.0), m24(0.0), m34(0.0), m44(m44)
    {};
  
    explicit double4x4(double4 v)
        : m11(v.x), m21(0.0), m31(0.0), m41(0.0),
          m12(0.0), m22(v.y), m32(0.0), m42(0.0),
          m13(0.0), m23(0.0), m33(v.z), m43(0.0),
          m14(0.0), m24(0.0), m34(0.0), m44(v.w)
    {};



    double4x4(double m11, double m21, double m31, double m41,
              double m12, double m22, double m32, double m42,
              double m13, double m23, double m33, double m43,
              double m14, double m24, double m34, double m44)
        : m11(m11), m21(m21), m31(m31), m41(m41),
        m12(m12), m22(m22), m32(m32), m42(m42),
        m13(m13), m23(m23), m33(m33), m43(m43),
        m14(m14), m24(m24), m34(m34), m44(m44)
    {

    }

    double4x4(double4 v1, double4 v2, double4 v3, double4 v4) : v{ v1,v2,v3,v4 } {}
    double4x4(const double3x3 &);

    double3x3 toDouble3x3() const
    {
        double3x3 w;
        w.m11 = this->m11; w.m21 = this->m21; w.m31 = this->m31;
        w.m12 = this->m12; w.m22 = this->m22; w.m32 = this->m32;
        w.m13 = this->m13; w.m23 = this->m23; w.m33 = this->m33;
        return w;
    }

    inline bool operator==(const double4x4& b) const
    {
        bool equal = true;
        const double epsilon = 1e-8;

        for (int i = 0; i < 16 && equal; i++)
            equal = fabs(this->m[i] - b.m[i]) <= epsilon;

        return equal;
    };
    inline bool operator!=(const double4x4& b) const { return !(*this == b); }

    double4x4 const transpose(void);
    static double4x4 const inverse(const double4x4& right);
    void inverse();
    static double4x4 TransformationAroundArbitraryPoint(double4x4 m, double3 p);

    static double4x4 Identity(void);
    static double4x4 AffinityMatrixToTransformationAroundArbitraryPoint(double4x4 m, double3 p);
    static double4x4 AffinityMatrixToTransformationAroundArbitraryPointWithTranslation(double4x4 m, double3 p, double3 t);
    static double4x4 Rotationdouble4x4FromQuaternion(simd_quatd q);
};


export inline double4x4 operator*(const double4x4& a, const double4x4& b)
{
    double4x4 r{};

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

export inline double4 operator*(const double4x4& a, const double4& b)
{
    double4 r{};

    r.x = a.m11 * b.x + a.m12 * b.y + a.m13 * b.z + a.m14 * b.w;
    r.y = a.m21 * b.x + a.m22 * b.y + a.m23 * b.z + a.m24 * b.w;
    r.z = a.m31 * b.x + a.m32 * b.y + a.m33 * b.z + a.m34 * b.w;
    r.w = a.m41 * b.x + a.m42 * b.y + a.m43 * b.z + a.m44 * b.w;

    return r;
}

export inline double4 operator*(const double4& a, const double4x4& b)
{
    double4 r{};

    r.x = a.x * b.m11 + a.y * b.m21 + a.z * b.m31 + a.w * b.m41;
    r.y = a.x * b.m12 + a.y * b.m22 + a.z * b.m32 + a.w * b.m42;
    r.z = a.x * b.m13 + a.y * b.m23 + a.z * b.m33 + a.w * b.m43;
    r.w = a.x * b.m14 + a.y * b.m24 + a.z * b.m34 + a.w * b.m44;

    return r;
}
