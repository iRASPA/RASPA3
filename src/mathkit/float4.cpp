module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#endif

module float4;

#ifdef USE_STD_IMPORT
import std;
#endif

void float4::normalise()
{
  float magnitude = std::sqrt((x * x) + (y * y) + (z * z) * (w * w));

  if (magnitude != 0)
  {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
    w /= magnitude;
  }
}
