module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <ostream>
#endif

export module float3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;

export union float3
{
  float v[3];
  struct
  {
    float x, y, z;
  };

  float3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}
  float3(int3 a) : x(float(a.x)), y(float(a.y)), z(float(a.z)) {}

  inline float& operator[](int i) { return v[i]; }
  inline const float& operator[](int i) const { return v[i]; }

  inline float length() { return std::sqrt(x * x + y * y + z * z); }
  inline float length() const { return std::sqrt(x * x + y * y + z * z); }
  float3 normalise();
  float3 fract();

  inline static float dot(const float3& v1, const float3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
  inline static float3 max(const float3& v1, const float3& v2)
  {
    return float3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
  }
  inline static float3 min(const float3& v1, const float3& v2)
  {
    return float3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
  }
  inline static float3 cross(const float3& v1, const float3& v2)
  {
    return float3(v1.y * v2.z - v1.z * v2.y, -v1.x * v2.z + v1.z * v2.x, v1.x * v2.y - v1.y * v2.x);
  }

  float3 operator-() const { return float3(-this->x, -this->y, -this->z); }
  float3& operator+=(const float3& b)
  {
    this->x += b.x, this->y += b.y, this->z += b.z;
    return *this;
  }
  float3& operator-=(const float3& b)
  {
    this->x -= b.x, this->y -= b.y, this->z -= b.z;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const float3& vec);
};

export inline float3 operator+(const float3& a, const float3& b) { return float3(a.x + b.x, a.y + b.y, a.z + b.z); }

export inline float3 operator-(const float3& a, const float3& b) { return float3(a.x - b.x, a.y - b.y, a.z - b.z); }

export inline float3 operator*(const float3& a, const float3& b) { return float3(a.x * b.x, a.y * b.y, a.z * b.z); }

export inline float3 operator*(const float3& a, float b) { return float3(a.x * b, a.y * b, a.z * b); }

export inline float3 operator*(float a, const float3& b) { return float3(a * b.x, a * b.y, a * b.z); }
