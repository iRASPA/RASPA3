module;

#ifdef USE_LEGACY_HEADERS
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <vector>
#include <tuple>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module rigid;

#ifndef USE_LEGACY_HEADERS
import <vector>;
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


std::pair<simd_quatd, simd_quatd> Rigid::NoSquishRotate(size_t k, double dt, std::pair<simd_quatd, simd_quatd> quat, double3 inverseInertiaVector)
{
  simd_quatd pn, qn;
  double zeta;

  auto& [p, q] = quat;
  switch(k)
  {
    case 1:
      zeta = dt * (-p.r * q.ix + p.ix * q.r + p.iy * q.iz - p.iz * q.iy) * inverseInertiaVector.x / 4.0;
      pn.r  = std::cos(zeta) * p.r  - std::sin(zeta) * p.ix;
      pn.ix = std::cos(zeta) * p.ix + std::sin(zeta) * p.r;
      pn.iy = std::cos(zeta) * p.iy + std::sin(zeta) * p.iz;
      pn.iz = std::cos(zeta) * p.iz - std::sin(zeta) * p.iy;

      qn.r  = std::cos(zeta) * q.r  - std::sin(zeta) * q.ix;
      qn.ix = std::cos(zeta) * q.ix + std::sin(zeta) * q.r;
      qn.iy = std::cos(zeta) * q.iy + std::sin(zeta) * q.iz;
      qn.iz = std::cos(zeta) * q.iz - std::sin(zeta) * q.iy;
      return {pn, qn};
    case 2:
      zeta = dt * (-p.r * q.iy - p.ix * q.iz + p.iy * q.r + p.iz * q.ix) * inverseInertiaVector.y / 4.0;
      pn.r  = std::cos(zeta) * p.r  - std::sin(zeta) * p.iy;
      pn.ix = std::cos(zeta) * p.ix - std::sin(zeta) * p.iz;
      pn.iy = std::cos(zeta) * p.iy + std::sin(zeta) * p.r;
      pn.iz = std::cos(zeta) * p.iz + std::sin(zeta) * p.ix;

      qn.r  = std::cos(zeta) * q.r  - std::sin(zeta) * q.iy;
      qn.ix = std::cos(zeta) * q.ix - std::sin(zeta) * q.iz;
      qn.iy = std::cos(zeta) * q.iy + std::sin(zeta) * q.r;
      qn.iz = std::cos(zeta) * q.iz + std::sin(zeta) * q.ix;
      return {pn, qn};
    case 3:
      zeta = dt * (-p.r * q.iz + p.ix * q.iy - p.iy * q.ix + p.iz * q.r) * inverseInertiaVector.z / 4.0;
      pn.r  = std::cos(zeta) * p.r  - std::sin(zeta) * p.iz;
      pn.ix = std::cos(zeta) * p.ix + std::sin(zeta) * p.iy;
      pn.iy = std::cos(zeta) * p.iy - std::sin(zeta) * p.ix;
      pn.iz = std::cos(zeta) * p.iz + std::sin(zeta) * p.r;

      qn.r  = std::cos(zeta) * q.r  - std::sin(zeta) * q.iz;
      qn.ix = std::cos(zeta) * q.ix + std::sin(zeta) * q.iy;
      qn.iy = std::cos(zeta) * q.iy - std::sin(zeta) * q.ix;
      qn.iz = std::cos(zeta) * q.iz + std::sin(zeta) * q.r;
      return {pn, qn};
      break;
    default:
      break;
  }
  return {};
}

std::pair<simd_quatd, simd_quatd>  Rigid::NoSquishFreeRotorOrderTwo(double dt, std::pair<simd_quatd, simd_quatd> q, double3 inverseInertiaVector)
{
  // second order
  for(size_t i = 0;i != 5; ++i)
  {
    q = NoSquishRotate(3, 0.5 * dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(2, 0.5 * dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(1, dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(2, 0.5 * dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(3, 0.5 * dt / 5.0, q, inverseInertiaVector);
  }
  return q;
}

