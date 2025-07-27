module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#endif

export module double4x3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double4;

export union double4x3
{
  double m[12];
  double mm[4][3];
  double3 v[4];
  struct
  {
    double m11, m21, m31, m12, m22, m32, m13, m23, m33, m14, m24, m34;
  };
  struct
  {
    double ax, ay, az, bx, by, bz, cx, cy, cz, wx, wy, wz;
  };

  double4x3(double3 v1, double3 v2, double3 v3, double3 v4) : v{v1, v2, v3, v4} {}

  inline double3& operator[](std::size_t i) { return v[i]; }
  inline const double3& operator[](std::size_t i) const { return v[i]; }
};
