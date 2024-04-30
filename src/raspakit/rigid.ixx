module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <tuple>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <tuple>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;

export module rigid;

export namespace Rigid
{
  std::pair<simd_quatd, simd_quatd> NoSquishFreeRotorOrderTwo(double dt, std::pair<simd_quatd, simd_quatd> q, double3 inverseInertiaVector);
}

namespace Rigid
{
  std::pair<simd_quatd, simd_quatd> NoSquishRotate(size_t k, double dt, std::pair<simd_quatd, simd_quatd> q, double3 inverseInertiaVector);
}

