module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cmath>
#endif

module float4;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
#endif

void float4::normalise()
{
  float magnitude = sqrt((x * x) + (y * y) + (z * z) * (w * w));

  if (magnitude != 0)
  {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
    w /= magnitude;
  }
}
