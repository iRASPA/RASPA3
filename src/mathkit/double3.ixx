module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#endif

export module double3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import bool3;
import hashcombine;
import archive;
import json;

export union double3
{
  //    #ifdef __AVX__
  //      __m256d vec;
  //    #endif
  double v[4];
  struct
  {
    double x, y, z, w;
  };

  //    #ifdef __AVX__
  //      double3(__m256d const& v): vec(v) {}
  //      double3& operator = (__m256d const& x) { vec = x; return *this; }
  //    #endif
  double3() : x(0.0), y(0.0), z(0.0), w(0.0) {}
  double3(double v) : x(v), y(v), z(v), w(0.0) {}
  double3(double x, double y, double z) : x(x), y(y), z(z), w(0.0) {}
  double3(int3 a) : x(double(a.x)), y(double(a.y)), z(double(a.z)), w(0.0) {}
  double3(double a[3]) : x(double(a[0])), y(double(a[1])), z(double(a[2])), w(0.0) {}

  bool operator==(double3 const& rhs) const { return (x == rhs.x) && (y == rhs.y) && (z == rhs.z); }

  inline double& operator[](std::size_t i) { return v[i]; }
  inline const double& operator[](std::size_t i) const { return v[i]; }

  double3 normalized();
  double3 fract() const;

  inline double length() { return std::sqrt(x * x + y * y + z * z); }
  inline double length_squared() const { return (x * x + y * y + z * z); }
  inline double length() const { return std::sqrt(x * x + y * y + z * z); }
  inline static double3 abs(double3 v1) { return double3(std::abs(v1.x), std::abs(v1.y), std::abs(v1.z)); }
  inline static double3 rint(double3 const& v) { return double3(std::rint(v.x), std::rint(v.y), std::rint(v.z)); }
  inline static double3 fract(double3 const& v)
  {
    return double3(v.x - std::floor(v.x), v.y - std::floor(v.y), v.z - std::floor(v.z));
  }
  inline static double dot(const double3& v1, const double3& v2)
  {
    // #ifdef __AVX__
    //         __m256d dp = _mm256_mul_pd(v1.vec, v2.vec);
    //         __m128d a = _mm256_extractf128_pd(dp, 0);
    //         __m128d b = _mm256_extractf128_pd(dp, 1);
    //         __m128d c = _mm_add_pd(a, b);
    //         __m128d yy = _mm_unpackhi_pd(c, c);
    //         __m128d dotproduct = _mm_add_pd(c, yy);
    //         return _mm_cvtsd_f64(dotproduct);
    // #else
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    // #endif
  }
  inline static double3 max(const double3& v1, const double3& v2)
  {
    return double3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
  }
  inline static double3 min(const double3& v1, const double3& v2)
  {
    return double3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
  }
  inline static double3 cross(const double3& v1, const double3& v2)
  {
    return double3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
  }
  inline static double3 normalize(const double3& v)
  {
    double f = 1.0 / std::sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    return double3(f * v.x, f * v.y, f * v.z);
  }
  inline static double3 flip(double3 v, bool3 flip, double3 boundary)
  {
    return double3(flip.x ? boundary.x - v.x : v.x, flip.y ? boundary.y - v.y : v.y, flip.z ? boundary.z - v.z : v.z);
  }
  inline void clamp(double3 low, double3 high)
  {
    x = std::clamp(x, low.x, high.x);
    y = std::clamp(y, low.y, high.y);
    z = std::clamp(z, low.z, high.z);
  }
  inline void clamp(double low, double high)
  {
    x = std::clamp(x, low, high);
    y = std::clamp(y, low, high);
    z = std::clamp(z, low, high);
  }

  inline static bool is_perpendicular(double3 u, double3 v, double precision = 1e-10)
  {
    return std::fabs(double3::dot(u, v)) < precision;
  }

  inline static double3 perpendicular(double3 u, double3 v)
  {
    double3 w;

    std::array<double, 3> denominator = {std::fabs(u.z * v.y - u.y * v.z), std::fabs(u.z * v.x - u.x * v.z),
                                         std::fabs(u.y * v.x - u.x * v.y)};

    std::array<double, 3>::iterator maximum = std::max_element(denominator.begin(), denominator.end());

    if (*maximum < 1e-10)
    {
      // u and v are parallel
      // https://math.stackexchange.com/questions/137362/how-to-find-perpendicular-vector-to-another-vector
      w = double3(std::copysign(u[2], u[0]), std::copysign(u[2], u[1]),
                  -std::copysign(u[0], u[2]) - std::copysign(u[1], u[2]));
      return w.normalized();
    }

    auto choise = std::distance(denominator.begin(), maximum);

    switch (choise)
    {
      case 0:
        // u * z == 0 && v * z == 0, solving for y and z (and x can be chosen freely)
        w.x = 1.0;
        w.y = (u.x * v.z - u.z * v.x) / (u.z * v.y - u.y * v.z);
        w.z = (u.y * v.x - u.x * v.y) / (u.z * v.y - u.y * v.z);
        break;
      case 1:
        // u * z == 0 && v * z == 0, solving for x and z (and y can be chosen freely)
        w.x = (u.y * v.z - u.z * v.y) / (u.z * v.x - u.x * v.z);
        w.y = 1.0;
        w.z = (u.x * v.y - u.y * v.x) / (u.z * v.x - u.x * v.z);
        break;
      case 2:
        // u * z == 0 && v * z == 0, solving for x and y (and z can be chosen freely)
        w.x = (u.z * v.y - u.y * v.z) / (u.y * v.x - u.x * v.y);
        w.y = (u.x * v.z - u.z * v.x) / (u.y * v.x - u.x * v.y);
        w.z = 1.0;
        break;
      default:
        std::unreachable();
    }

    return w.normalized();
  }

  inline static double angle(double3 v, double3 w)
  {
    return std::atan2(double3::cross(v, w).length(), double3::dot(v, w));
  }

  inline static double angle(double3 posA, double3 posB, double3 posC)
  {
    double3 Rab = posA - posB;
    double3 Rcb = posC - posB;
    return std::atan2(double3::cross(Rab, Rcb).length(), double3::dot(Rab, Rcb));
  }

  inline double3 rotateAroundAxis(const double3& coord, double angle) const
  {
    double c = std::cos(angle);
    double w = 1.0 - c;
    double s = std::sqrt(1.0 - c * c);

    if (angle < 0.0) s = -s;

    double r1 = coord.x * (this->x * this->x * w + c) + coord.y * (this->x * this->y * w + this->z * s) +
                coord.z * (this->x * this->z * w - this->y * s);
    double r2 = coord.x * (this->x * this->y * w - this->z * s) + coord.y * (this->y * this->y * w + c) +
                coord.z * (this->y * this->z * w + this->x * s);
    double r3 = coord.x * (this->x * this->z * w + this->y * s) + coord.y * (this->y * this->z * w - this->x * s) +
                coord.z * (this->z * this->z * w + c);

    return {r1, r2, r3};
  }

  inline double3 operator+(const double3& b) const { return double3(this->x + b.x, this->y + b.y, this->z + b.z); }

  inline double3 operator-(const double3& b) const { return double3(this->x - b.x, this->y - b.y, this->z - b.z); }

  double3 operator-() const { return double3(-this->x, -this->y, -this->z); }
  double3& operator+=(const double3& b)
  {
    this->x += b.x, this->y += b.y, this->z += b.z;
    return *this;
  }
  double3& operator-=(const double3& b)
  {
    this->x -= b.x, this->y -= b.y, this->z -= b.z;
    return *this;
  }
  double3& operator*=(const double& b)
  {
    this->x *= b, this->y *= b, this->z *= b;
    return *this;
  }
  double3& operator/=(const double& b)
  {
    this->x /= b, this->y /= b, this->z /= b;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const double3& vec);

  struct KeyHash
  {
    std::size_t operator()(const double3& k) const
    {
      std::size_t h1 = std::hash<double>()(k.x);
      std::size_t h2 = std::hash<double>()(k.y);
      std::size_t h3 = std::hash<double>()(k.z);
      return (h1 ^ (h2 << 1)) ^ h3;
    }
  };

  struct KeyEqual
  {
    bool operator()(const double3& lhs, const double3& rhs) const
    {
      return (std::fabs(lhs.x - rhs.x) < 1e-3) && std::fabs(lhs.y - rhs.y) < 1e-3 && std::fabs(lhs.z == rhs.z) < 1e-3;
    }
  };

  // const Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive) const
  //{
  //   archive << x << y << z;
  //   return archive;
  // }

  // Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive)
  //{
  //   archive >> x >> y >> z;
  //   return archive;
  // }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const double3& vec);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, double3& vec);

  std::string to_string();
  friend void to_json(nlohmann::json&, const double3&);
  friend void from_json(const nlohmann::json&, double3&);
};

export inline double3 operator*(const double3& a, const double3& b) { return double3(a.x * b.x, a.y * b.y, a.z * b.z); }

export inline double3 operator/(const double3& a, const double3& b) { return double3(a.x / b.x, a.y / b.y, a.z / b.z); }

export inline double3 operator*(const double3& a, double b) { return double3(a.x * b, a.y * b, a.z * b); }

export inline double3 operator*(const double& a, const double3& b) { return double3(a * b.x, a * b.y, a * b.z); }

export inline double3 operator/(const double3& a, double b) { return double3(a.x / b, a.y / b, a.z / b); }

export inline double3 sqrt(const double3& a) { return double3(std::sqrt(a.x), std::sqrt(a.y), std::sqrt(a.z)); }

export double3 clamp(double3 value, double3 low, double3 high)
{
  return double3(std::clamp(value.x, low.x, high.x), std::clamp(value.y, low.y, high.y),
                 std::clamp(value.z, low.z, high.z));
}

export template <>
struct std::formatter<double3> : std::formatter<std::string_view>
{
  auto format(const double3& v, std::format_context& ctx) const
  {
    std::string temp{};
    std::format_to(std::back_inserter(temp), "({}, {}, {})", v.x, v.y, v.z);
    return std::formatter<std::string_view>::format(temp, ctx);
  }
};

export void to_json(nlohmann::json& j, const double3& a) { j = nlohmann::json{a.x, a.y, a.z}; }

export void from_json(const nlohmann::json& j, double3& a)
{
  std::array<double, 3> data;
  j.get_to(data);
  a = double3(data[0], data[1], data[2]);
}
