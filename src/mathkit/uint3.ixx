module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <functional>
#include <tuple>
#endif

export module uint3;

#ifdef USE_STD_IMPORT
import std;
#endif

import hashcombine;
import archive;

export union uint3
{
  std::size_t v[3];
  struct
  {
    std::size_t x, y, z;
  };

  uint3() : x(0), y(0), z(0) {}
  uint3(std::size_t x, std::size_t y, std::size_t z) : x(x), y(y), z(z) {}
  inline std::size_t& operator[](int i) { return v[i]; }
  inline const std::size_t& operator[](int i) const { return v[i]; }

  inline std::size_t length_squared() const { return (x * x + y * y + z * z); }

  inline uint3& operator+=(const uint3& v1)
  {
    this->x += v1.x;
    this->y += v1.y;
    this->z += v1.z;
    return *this;
  }
  inline uint3& operator-=(const uint3& v1)
  {
    this->x -= v1.x;
    this->y -= v1.y;
    this->z -= v1.z;
    return *this;
  }

  inline uint3& operator*=(const std::size_t& v1)
  {
    this->x *= v1;
    this->y *= v1;
    this->z *= v1;
    return *this;
  }
  inline uint3& operator/=(const std::size_t& v1)
  {
    this->x /= v1;
    this->y /= v1;
    this->z /= v1;
    return *this;
  }

  inline uint3 operator-() const { return uint3(-this->x, -this->y, -this->z); }

  inline bool operator==(const uint3& rhs) const
  {
    return (this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z);
  }

  struct hashFunction
  {
    std::size_t operator()(const uint3& k) const
    {
      std::size_t h = 0;
      hash_combine(h, k.x);
      hash_combine(h, k.y);
      hash_combine(h, k.z);
      return h;
    }
  };
  // friend QDebug operator<<(QDebug debug, const uint3& m);
  // friend QDataStream& operator<<(QDataStream&, const uint3&);
  // friend QDataStream& operator>>(QDataStream&, uint3&);

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const uint3& vec);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, uint3& vec);

};

export inline uint3 operator+(uint3 lhs, const uint3& rhs) { return uint3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z); }
export inline uint3 operator-(uint3 lhs, const uint3& rhs) { return uint3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z); }
export inline uint3 operator*(uint3 lhs, const std::size_t& rhs) { return uint3(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs); }
export inline uint3 operator/(uint3 lhs, const std::size_t& rhs) { return uint3(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs); }

export inline bool operator<(const uint3& lhs, const uint3& rhs)
{
  return std::tie(lhs.x, lhs.y, lhs.z) < std::tie(rhs.x, rhs.y, rhs.z);
}

export inline bool operator<(const uint3 lhs, const uint3 rhs)
{
  return std::tie(lhs.x, lhs.y, lhs.z) < std::tie(rhs.x, rhs.y, rhs.z);
}

template <>
struct std::formatter<uint3> : std::formatter<std::string_view>
{
  auto format(const uint3& v, std::format_context& ctx) const
  {
    std::string temp{};
    std::format_to(std::back_inserter(temp), "({}, {}, {})", v.x, v.y, v.z);
    return std::formatter<std::string_view>::format(temp, ctx);
  }
};

/*
namespace std
{
    export template <> struct hash<int3>
    {
        size_t operator()(const uint3& k) const
        {
            std::size_t h = 0;
            hash_combine(h, k.x);
            hash_combine(h, k.y);
            hash_combine(h, k.z);
            return h;
        }
    };
}*/
