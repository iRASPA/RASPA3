module;

#ifdef USE_LEGACY_HEADERS
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <tuple>
#endif

export module int3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import hashcombine;
import archive;

export union int3
{
  std::int32_t v[3];
  struct
  {
    std::int32_t x, y, z;
  };

  int3() : x(0), y(0), z(0) {}
  int3(int x, int y, int z) : x(x), y(y), z(z) {}
  inline int& operator[](int i) { return v[i]; }
  inline const int& operator[](int i) const { return v[i]; }

  inline int length_squared() const { return (x * x + y * y + z * z); }

  static int3 greatestCommonDivisor(int3 a, int b);

  inline int3& operator+=(const int3& v1)
  {
    this->x += v1.x;
    this->y += v1.y;
    this->z += v1.z;
    return *this;
  }
  inline int3& operator-=(const int3& v1)
  {
    this->x -= v1.x;
    this->y -= v1.y;
    this->z -= v1.z;
    return *this;
  }

  inline int3& operator*=(const int& v1)
  {
    this->x *= v1;
    this->y *= v1;
    this->z *= v1;
    return *this;
  }
  inline int3& operator/=(const int& v1)
  {
    this->x /= v1;
    this->y /= v1;
    this->z /= v1;
    return *this;
  }

  inline int3 operator-() const { return int3(-this->x, -this->y, -this->z); }

  inline bool operator==(const int3& rhs) const
  {
    return (this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z);
  }

  struct hashFunction
  {
    std::size_t operator()(const int3& k) const
    {
      std::size_t h = 0;
      hash_combine(h, k.x);
      hash_combine(h, k.y);
      hash_combine(h, k.z);
      return h;
    }
  };
  // friend QDebug operator<<(QDebug debug, const int3& m);
  // friend QDataStream& operator<<(QDataStream&, const int3&);
  // friend QDataStream& operator>>(QDataStream&, int3&);

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const int3& vec);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, int3& vec);
};

export inline int3 operator+(int3 lhs, const int3& rhs) { return int3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z); }
export inline int3 operator-(int3 lhs, const int3& rhs) { return int3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z); }
export inline int3 operator*(int3 lhs, const int& rhs) { return int3(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs); }
export inline int3 operator/(int3 lhs, const int& rhs) { return int3(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs); }

export inline bool operator<(const int3& lhs, const int3& rhs)
{
  return std::tie(lhs.x, lhs.y, lhs.z) < std::tie(rhs.x, rhs.y, rhs.z);
}

export inline bool operator<(const int3 lhs, const int3 rhs)
{
  return std::tie(lhs.x, lhs.y, lhs.z) < std::tie(rhs.x, rhs.y, rhs.z);
}

/*
namespace std
{
    export template <> struct hash<int3>
    {
        size_t operator()(const int3& k) const
        {
            std::size_t h = 0;
            hash_combine(h, k.x);
            hash_combine(h, k.y);
            hash_combine(h, k.z);
            return h;
        }
    };
}*/
