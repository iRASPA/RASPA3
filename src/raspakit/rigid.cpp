module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <format>
#include <print>
#include <tuple>
#include <vector>
#endif

module rigid;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;

std::pair<simd_quatd, simd_quatd> Rigid::NoSquishRotate(std::size_t k, double dt,
                                                        std::pair<simd_quatd, simd_quatd> quat,
                                                        double3 inverseInertiaVector)
{
  simd_quatd pn, qn;
  double zeta, cos_zeta, sin_zeta;

  auto& [p, q] = quat;
  switch (k)
  {
    case 1:
      zeta = dt * (-p.r * q.ix + p.ix * q.r + p.iy * q.iz - p.iz * q.iy) * inverseInertiaVector.x / 4.0;
      cos_zeta = std::cos(zeta);
      sin_zeta = std::sin(zeta);

      pn.r = cos_zeta * p.r - sin_zeta * p.ix;
      pn.ix = cos_zeta * p.ix + sin_zeta * p.r;
      pn.iy = cos_zeta * p.iy + sin_zeta * p.iz;
      pn.iz = cos_zeta * p.iz - sin_zeta * p.iy;

      qn.r = cos_zeta * q.r - sin_zeta * q.ix;
      qn.ix = cos_zeta * q.ix + sin_zeta * q.r;
      qn.iy = cos_zeta * q.iy + sin_zeta * q.iz;
      qn.iz = cos_zeta * q.iz - sin_zeta * q.iy;
      return {pn, qn};
    case 2:
      zeta = dt * (-p.r * q.iy - p.ix * q.iz + p.iy * q.r + p.iz * q.ix) * inverseInertiaVector.y / 4.0;
      cos_zeta = std::cos(zeta);
      sin_zeta = std::sin(zeta);

      pn.r = cos_zeta * p.r - sin_zeta * p.iy;
      pn.ix = cos_zeta * p.ix - sin_zeta * p.iz;
      pn.iy = cos_zeta * p.iy + sin_zeta * p.r;
      pn.iz = cos_zeta * p.iz + sin_zeta * p.ix;

      qn.r = cos_zeta * q.r - sin_zeta * q.iy;
      qn.ix = cos_zeta * q.ix - sin_zeta * q.iz;
      qn.iy = cos_zeta * q.iy + sin_zeta * q.r;
      qn.iz = cos_zeta * q.iz + sin_zeta * q.ix;
      return {pn, qn};
    case 3:
      zeta = dt * (-p.r * q.iz + p.ix * q.iy - p.iy * q.ix + p.iz * q.r) * inverseInertiaVector.z / 4.0;
      cos_zeta = std::cos(zeta);
      sin_zeta = std::sin(zeta);

      pn.r = cos_zeta * p.r - sin_zeta * p.iz;
      pn.ix = cos_zeta * p.ix + sin_zeta * p.iy;
      pn.iy = cos_zeta * p.iy - sin_zeta * p.ix;
      pn.iz = cos_zeta * p.iz + sin_zeta * p.r;

      qn.r = cos_zeta * q.r - sin_zeta * q.iz;
      qn.ix = cos_zeta * q.ix + sin_zeta * q.iy;
      qn.iy = cos_zeta * q.iy - sin_zeta * q.ix;
      qn.iz = cos_zeta * q.iz + sin_zeta * q.r;
      return {pn, qn};
      break;
    default:
      break;
  }
  return {};
}

std::pair<simd_quatd, simd_quatd> Rigid::NoSquishFreeRotorOrderTwo(double dt, std::pair<simd_quatd, simd_quatd> q,
                                                                   double3 inverseInertiaVector)
{
  // second order
  for (std::size_t i = 0; i != 5; ++i)
  {
    q = NoSquishRotate(3, 0.5 * dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(2, 0.5 * dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(1, dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(2, 0.5 * dt / 5.0, q, inverseInertiaVector);
    q = NoSquishRotate(3, 0.5 * dt / 5.0, q, inverseInertiaVector);
  }
  return q;
}
