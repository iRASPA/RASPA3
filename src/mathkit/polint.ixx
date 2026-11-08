module;

export module polint;

import std;

export namespace Interpolation
{
template <std::size_t N>
void polint(const std::array<double, N> &xa, const std::array<double, N> &ya, double x, double *y, double *dy)
{
  std::size_t i, m, nearest = 0;
  double den, dif, dift, ho, hp, w;
  std::array<double, N> c, d;

  dif = std::fabs(x - xa[0]);
  for (i = 0; i < N; i++)
  {
    if ((dift = std::fabs(x - xa[i])) < dif)
    {
      nearest = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[nearest];

  // 'next' is one past the index of the last correction taken along the tableau (the 1-based 'ns' of Numerical
  // Recipes); the 'c'-branch is always taken when it is zero, since 2 * 0 < N - m, so it never goes below zero
  std::size_t next = nearest;
  for (m = 1; m < N; m++)
  {
    for (i = 0; i < N - m; i++)
    {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      if ((den = ho - hp) == 0.0)
      {
        throw std::invalid_argument("Interpolation::polint: two abscissae are equal");
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    if (2 * next < (N - m))
    {
      *dy = c[next];
    }
    else
    {
      *dy = d[next - 1];
      --next;
    }
    *y += *dy;
  }
}

}  // namespace Interpolation
