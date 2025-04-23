module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <ostream>
#endif

export module double3x3;

// #if defined(WIN32)
//     import <intrin.h>;
// #elif defined(__AVX__)
//     import <immintrin.h>;
// #endif

#ifndef USE_LEGACY_HEADERS
import <istream>;
import <ostream>;
import <fstream>;
import <cmath>;
#endif

import archive;
import double3;
import simd_quatd;
import int3x3;
import json;

export union double3x3
{
  double m[16];
  double mm[4][4];
  double3 v[4];
  struct
  {
    double m11, m21, m31, m41,  // 1st column
        m12, m22, m32, m42,     // 2nd column
        m13, m23, m33, m43,     // 3rd column
        m14, m24, m34, m44;     // 4th column, {tx,ty,tz,1}
  };
  struct
  {
    double ax, ay, az, aw,  // 1st column
        bx, by, bz, bw,     // 2nd column
        cx, cy, cz, cw,     // 3rd column
        wx, wy, wz, ww;     // 4th column, {tx,ty,tz,1}
  };

  double3x3()
      : m11(0.0),
        m21(0.0),
        m31(0.0),
        m41(0.0),
        m12(0.0),
        m22(0.0),
        m32(0.0),
        m42(0.0),
        m13(0.0),
        m23(0.0),
        m33(0.0),
        m43(0.0),
        m14(0.0),
        m24(0.0),
        m34(0.0),
        m44(0.0) {};

  explicit double3x3(double m11, double m22, double m33)
      : m11(m11),
        m21(0.0),
        m31(0.0),
        m41(0.0),
        m12(0.0),
        m22(m22),
        m32(0.0),
        m42(0.0),
        m13(0.0),
        m23(0.0),
        m33(m33),
        m43(0.0),
        m14(0.0),
        m24(0.0),
        m34(0.0),
        m44(0.0) {};

  explicit double3x3(double3 v)
      : m11(v.x),
        m21(0.0),
        m31(0.0),
        m41(0.0),
        m12(0.0),
        m22(v.y),
        m32(0.0),
        m42(0.0),
        m13(0.0),
        m23(0.0),
        m33(v.z),
        m43(0.0),
        m14(0.0),
        m24(0.0),
        m34(0.0),
        m44(0.0) {};

  explicit double3x3(double m11, double m21, double m31, double m12, double m22, double m32, double m13, double m23,
                     double m33)
      : m11(m11),
        m21(m21),
        m31(m31),
        m41(0.0),
        m12(m12),
        m22(m22),
        m32(m32),
        m42(0.0),
        m13(m13),
        m23(m23),
        m33(m33),
        m43(0.0),
        m14(0.0),
        m24(0.0),
        m34(0.0),
        m44(0.0) {

        };

  explicit double3x3(double3 v1, double3 v2, double3 v3)
      : m11(v1.x),
        m21(v1.y),
        m31(v1.z),
        m41(0.0),
        m12(v2.x),
        m22(v2.y),
        m32(v2.z),
        m42(0.0),
        m13(v3.x),
        m23(v3.y),
        m33(v3.z),
        m43(0.0),
        m14(0.0),
        m24(0.0),
        m34(0.0),
        m44(0.0) {

        };
  double3x3(simd_quatd q);

  double3x3(const int3x3& m);

  double3x3(double lattice[3][3]);

  bool operator==(double3x3 const& rhs)
  {
    return (ax == rhs.ax) && (ay == rhs.ay) && (az == rhs.az) && (bx == rhs.bx) && (by == rhs.by) && (bz == rhs.bz) &&
           (cx == rhs.cx) && (cy == rhs.cy) && (cz == rhs.cz);
  }

  inline double3& operator[](size_t i) { return v[i]; }
  inline const double3& operator[](size_t i) const { return v[i]; }

  static double3x3 identity();

  double determinant(void) const;
  double trace(void) const;
  double3x3 const inverse() const;
  static double3x3 inverse(const double3x3& right);
  static double3x3 transpose(const double3x3& right);
  double3x3 const transpose(void) const;
  double3x3 inversetranpose(void);
  void EigenSystemSymmetric(double3& eigenvalues, double3x3& eigenvectors);

  static double3x3 buildRotationMatrix(const simd_quatd& q);
  static double3x3 buildRotationMatrixInverse(const simd_quatd& q);
  simd_quatd quaternion();

  inline double3x3 operator-() const { return double3x3(-this->v[0], -this->v[1], -this->v[2]); }
  inline bool operator==(const double3x3& b) const
  {
    return ((std::fabs(this->m11 - b.m11) < 1e-5) && (std::fabs(this->m12 - b.m12) < 1e-5) &&
            (std::fabs(this->m13 - b.m13) < 1e-5) && (std::fabs(this->m21 - b.m21) < 1e-5) &&
            (std::fabs(this->m22 - b.m22) < 1e-5) && (std::fabs(this->m23 - b.m23) < 1e-5) &&
            (std::fabs(this->m31 - b.m31) < 1e-5) && (std::fabs(this->m32 - b.m32) < 1e-5) &&
            (std::fabs(this->m33 - b.m33) < 1e-5));
  }

  inline int3x3 toInt3x3()
  {
    int3x3 w;

    w.m11 = int(rint(this->m11));
    w.m21 = int(rint(this->m21));
    w.m31 = int(rint(this->m31));
    w.m12 = int(rint(this->m12));
    w.m22 = int(rint(this->m22));
    w.m32 = int(rint(this->m32));
    w.m13 = int(rint(this->m13));
    w.m23 = int(rint(this->m23));
    w.m33 = int(rint(this->m33));
    return w;
  }

  inline double3x3& operator+=(const double3x3& b)
  {
    this->m11 += b.m11;
    this->m21 += b.m21;
    this->m31 += b.m31;
    this->m12 += b.m12;
    this->m22 += b.m22;
    this->m32 += b.m32;
    this->m13 += b.m13;
    this->m23 += b.m23;
    this->m33 += b.m33;

    return *this;
  }

  inline double3x3& operator-=(const double3x3& b)
  {
    this->m11 -= b.m11;
    this->m21 -= b.m21;
    this->m31 -= b.m31;
    this->m12 -= b.m12;
    this->m22 -= b.m22;
    this->m32 -= b.m32;
    this->m13 -= b.m13;
    this->m23 -= b.m23;
    this->m33 -= b.m33;

    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const double3x3& vec);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, double3x3& vec);

  friend void to_json(nlohmann::json&, const double3x3&);
  friend void from_json(const nlohmann::json&, double3x3&);
};

export inline double3x3 operator+(const double3x3& a, const double3x3& b)
{
  double3x3 r{};

  r.m11 = a.m11 + b.m11;
  r.m12 = a.m12 + b.m12;
  r.m13 = a.m13 + b.m13;
  r.m21 = a.m21 + b.m21;
  r.m22 = a.m22 + b.m22;
  r.m23 = a.m23 + b.m23;
  r.m31 = a.m31 + b.m31;
  r.m32 = a.m32 + b.m32;
  r.m33 = a.m33 + b.m33;
  return r;
}

export inline double3x3 operator-(const double3x3& a, const double3x3& b)
{
  double3x3 r{};

  r.m11 = a.m11 - b.m11;
  r.m12 = a.m12 - b.m12;
  r.m13 = a.m13 - b.m13;
  r.m21 = a.m21 - b.m21;
  r.m22 = a.m22 - b.m22;
  r.m23 = a.m23 - b.m23;
  r.m31 = a.m31 - b.m31;
  r.m32 = a.m32 - b.m32;
  r.m33 = a.m33 - b.m33;
  return r;
}

export inline double3x3 operator*(const double3x3& a, const double3x3& b)
{
  double3x3 r{};

  r.m11 = a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31;
  r.m21 = a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31;
  r.m31 = a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31;

  r.m12 = a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32;
  r.m22 = a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32;
  r.m32 = a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32;

  r.m13 = a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33;
  r.m23 = a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33;
  r.m33 = a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33;

  return r;
}

export inline double3 operator*(const double3x3& a, const double3& b)
{
  double3 r{};

  r.x = a.m11 * b.x + a.m12 * b.y + a.m13 * b.z;
  r.y = a.m21 * b.x + a.m22 * b.y + a.m23 * b.z;
  r.z = a.m31 * b.x + a.m32 * b.y + a.m33 * b.z;

  return r;
}

export inline double3x3 sqr(const double3x3& a)
{
  double3x3 r{};

  r.m11 = a.m11 * a.m11;
  r.m21 = a.m21 * a.m21;
  r.m31 = a.m31 * a.m31;
  r.m12 = a.m12 * a.m12;
  r.m22 = a.m22 * a.m22;
  r.m32 = a.m32 * a.m32;
  r.m13 = a.m13 * a.m13;
  r.m23 = a.m23 * a.m23;
  r.m33 = a.m33 * a.m33;

  return r;
}

export inline double3 transposedMultiply(const double3x3& a, const double3& b)
{
  double3 r{};

  r.x = a.m11 * b.x + a.m21 * b.y + a.m31 * b.z;
  r.y = a.m12 * b.x + a.m22 * b.y + a.m32 * b.z;
  r.z = a.m13 * b.x + a.m23 * b.y + a.m33 * b.z;

  return r;
}

export inline double3 operator*(const double3x3& a, const std::array<double, 3>& b)
{
  double3 r{};

  r.x = a.m11 * b[0] + a.m12 * b[1] + a.m13 * b[2];
  r.y = a.m21 * b[0] + a.m22 * b[1] + a.m23 * b[2];
  r.z = a.m31 * b[0] + a.m32 * b[1] + a.m33 * b[2];

  return r;
}

// std::array<std::array<double,3>,3> is row-based
export inline double3x3 operator*(const double3x3& a, const std::array<std::array<double, 3>, 3>& b)
{
  double3x3 r{};

  r.m11 = a.m11 * b[0][0] + a.m12 * b[0][1] + a.m13 * b[0][2];
  r.m21 = a.m21 * b[0][0] + a.m22 * b[0][1] + a.m23 * b[0][2];
  r.m31 = a.m31 * b[0][0] + a.m32 * b[0][1] + a.m33 * b[0][2];

  r.m12 = a.m11 * b[1][0] + a.m12 * b[1][1] + a.m13 * b[1][2];
  r.m22 = a.m21 * b[1][0] + a.m22 * b[1][1] + a.m23 * b[1][2];
  r.m32 = a.m31 * b[1][0] + a.m32 * b[1][1] + a.m33 * b[1][2];

  r.m13 = a.m11 * b[2][0] + a.m12 * b[2][1] + a.m13 * b[2][2];
  r.m23 = a.m21 * b[2][0] + a.m22 * b[2][1] + a.m23 * b[2][2];
  r.m33 = a.m31 * b[2][0] + a.m32 * b[2][1] + a.m33 * b[2][2];

  return r;
}

// std::array<std::array<double,3>,3> is row-based
export inline double3x3 operator*(const std::array<std::array<double, 3>, 3>& a, const double3x3& b)
{
  double3x3 r{};

  r.m11 = a[0][0] * b.m11 + a[1][0] * b.m21 + a[2][0] * b.m31;
  r.m21 = a[0][1] * b.m11 + a[1][1] * b.m21 + a[2][1] * b.m31;
  r.m31 = a[0][2] * b.m11 + a[1][2] * b.m21 + a[2][2] * b.m31;

  r.m12 = a[0][0] * b.m12 + a[1][0] * b.m22 + a[2][0] * b.m32;
  r.m22 = a[0][1] * b.m12 + a[1][1] * b.m22 + a[2][1] * b.m32;
  r.m32 = a[0][2] * b.m12 + a[1][2] * b.m22 + a[2][2] * b.m32;

  r.m13 = a[0][0] * b.m13 + a[1][0] * b.m23 + a[2][0] * b.m33;
  r.m23 = a[0][1] * b.m13 + a[1][1] * b.m23 + a[2][1] * b.m33;
  r.m33 = a[0][2] * b.m13 + a[1][2] * b.m23 + a[2][2] * b.m33;

  return r;
}

export inline double3x3 operator*(const double& a, const double3x3& b)
{
  double3x3 r{};

  r.m11 = a * b.m11;
  r.m12 = a * b.m12;
  r.m13 = a * b.m13;
  r.m21 = a * b.m21;
  r.m22 = a * b.m22;
  r.m23 = a * b.m23;
  r.m31 = a * b.m31;
  r.m32 = a * b.m32;
  r.m33 = a * b.m33;

  return r;
}

export inline double3x3 operator*(const double3x3& a, const double& b)
{
  double3x3 r{};

  r.m11 = a.m11 * b;
  r.m12 = a.m12 * b;
  r.m13 = a.m13 * b;
  r.m21 = a.m21 * b;
  r.m22 = a.m22 * b;
  r.m23 = a.m23 * b;
  r.m31 = a.m31 * b;
  r.m32 = a.m32 * b;
  r.m33 = a.m33 * b;

  return r;
}

export inline double3x3 operator/(const double3x3& a, const double& b)
{
  double3x3 r{};
  r.m11 = a.m11 / b;
  r.m12 = a.m12 / b;
  r.m13 = a.m13 / b;
  r.m21 = a.m21 / b;
  r.m22 = a.m22 / b;
  r.m23 = a.m23 / b;
  r.m31 = a.m31 / b;
  r.m32 = a.m32 / b;
  r.m33 = a.m33 / b;
  return r;
}

export inline double3x3 sqrt(const double3x3& b)
{
  double3x3 r{};
  r.m11 = std::sqrt(b.m11);
  r.m12 = std::sqrt(b.m12);
  r.m13 = std::sqrt(b.m13);
  r.m21 = std::sqrt(b.m21);
  r.m22 = std::sqrt(b.m22);
  r.m23 = std::sqrt(b.m23);
  r.m31 = std::sqrt(b.m31);
  r.m32 = std::sqrt(b.m32);
  r.m33 = std::sqrt(b.m33);
  return r;
}

export void to_json(nlohmann::json& j, const double3x3& a)
{
  j = nlohmann::json{{a.ax, a.ay, a.az}, {a.bx, a.by, a.bz}, {a.cx, a.cy, a.cz}};
}

export void from_json(const nlohmann::json& j, double3x3& a) { j.get_to(a); }
