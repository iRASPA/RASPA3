module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>
#include <print>
#endif

#ifdef USE_STD_IMPORT
#include <stdio.h>
#endif

export module polint;

#ifdef USE_STD_IMPORT
import std;
#endif

export namespace Interpolation
{
template <std::size_t N>
void polint(const std::array<double, N> &xa, const std::array<double, N> &ya, double x, double *y, double *dy)
{
  std::size_t i, m, ns = 0;
  double den, dif, dift, ho, hp, w;
  std::array<double, N> c, d;

  dif = std::fabs(x - xa[0]);
  for (i = 0; i < N; i++)
  {
    if ((dift = std::fabs(x - xa[i])) < dif)
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[ns];
  --ns;
  for (m = 1; m < N; m++)
  {
    for (i = 0; i < N - m; i++)
    {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      if ((den = ho - hp) == 0.0)
      {
        fprintf(stderr, "Error in routine polint\n");
        std::exit(0);
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    if (2 * (ns + 1) < (N - m))
    {
      *dy = c[ns + 1];
    }
    else
    {
      *dy = d[ns];
      --ns;
    }
    *y += *dy;
  }
}

}  // namespace Interpolation
