module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <print>
#include <tuple>
#endif

export module rigid;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <tuple>;
import <print>;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;

export namespace Rigid
{
std::pair<simd_quatd, simd_quatd> NoSquishFreeRotorOrderTwo(double dt, std::pair<simd_quatd, simd_quatd> q,
                                                            double3 inverseInertiaVector);
std::pair<simd_quatd, simd_quatd> NoSquishRotate(size_t k, double dt, std::pair<simd_quatd, simd_quatd> q,
                                                 double3 inverseInertiaVector);
}  // namespace Rigid
