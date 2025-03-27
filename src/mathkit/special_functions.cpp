module;

#ifdef USE_LEGACY_HEADERS
#include <bitset>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <limits>
#endif

module special_functions;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
import <cstring>;
import <bitset>;
import <limits>;
#endif

// routine by Alexander Voigt
// https://arxiv.org/abs/2201.01678
// "Comparison of methods for the calculation of the real dilogarithm regarding instruction-level parallelism"
// https://arxiv.org/src/2201.01678v1/anc/Li2.c
double li2(double x)
{
  const double PI = 3.1415926535897932;
  const double P[] = {0.9999999999999999502e+0,  -2.6883926818565423430e+0, 2.6477222699473109692e+0,
                      -1.1538559607887416355e+0, 2.0886077795020607837e-1,  -1.0859777134152463084e-2};
  const double Q[] = {1.0000000000000000000e+0,  -2.9383926818565635485e+0, 3.2712093293018635389e+0,
                      -1.7076702173954289421e+0, 4.1596017228400603836e-1,  -3.9801343754084482956e-2,
                      8.2743668974466659035e-4};

  double y = 0.0, r = 0.0, s = 1.0;

  if (x < -1.0)
  {
    const double l = std::log(1.0 - x);
    y = 1.0 / (1.0 - x);
    r = -PI * PI / 6.0 + l * (0.5 * l - std::log(-x));
    s = 1;
  }
  else if (x == -1.0)
  {
    return -PI * PI / 12.0;
  }
  else if (x < 0.0)
  {
    const double l = std::log1p(-x);
    y = x / (x - 1);
    r = -0.5 * l * l;
    s = -1;
  }
  else if (x == 0.0)
  {
    return 0;
  }
  else if (x < 0.5)
  {
    y = x;
    r = 0.0;
    s = 1.0;
  }
  else if (x < 1.0)
  {
    y = 1.0 - x;
    r = PI * PI / 6.0 - std::log(x) * std::log(y);
    s = -1.0;
  }
  else if (x == 1.0)
  {
    return PI * PI / 6.0;
  }
  else if (x < 2.0)
  {
    const double l = std::log(x);
    y = 1.0 - 1.0 / x;
    r = PI * PI / 6.0 - l * (std::log(y) + 0.5 * l);
    s = 1.0;
  }
  else
  {
    const double l = std::log(x);
    y = 1.0 / x;
    r = PI * PI / 3.0 - 0.5 * l * l;
    s = -1.0;
  }

  const double y2 = y * y;
  const double y4 = y2 * y2;
  const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]) + y4 * (P[4] + y * P[5]);
  const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]) + y4 * (Q[4] + y * Q[5] + y2 * Q[6]);

  return r + s * y * p / q;
}

// Note convergence restrictions: abs(x) < 1 and c not a negative integer or zero
double hypergeometric(double a, double b, double c, double x)
{
  const double TOLERANCE = 1.0e-3;
  double term = a * b * x / c;
  double value = 1.0 + term;
  int n = 1;

  while (std::fabs(term) > TOLERANCE)
  {
    a++, b++, c++, n++;
    term *= a * b * x / c / n;
    value += term;
  }

  return value;
}

int areClose(const double x1, const double x2, const double epsilon)
{
  int exponent;
  double delta, difference, maxXY;

  /* Find exponent of largest absolute value */
  maxXY = (std::fabs(x1) > std::fabs(x2)) ? x1 : x2;
  std::frexp(maxXY, &exponent);

  /* Form a neighborhood of size  2 * delta */
  delta = std::ldexp(epsilon, exponent);
  difference = x1 - x2;

  if (difference > delta) /* x1 > x2 */
    return 0;
  else if (difference < -delta) /* x1 < x2 */
    return 0;
  else        /* -delta <= difference <= delta */
    return 1; /* x1 ~=~ x2 */
}

int areEqual(double a, double b)
{
  if (a == b) return 1;

  if (areClose(a, b, std::numeric_limits<double>::epsilon())) return 1;

  return 0;
}

double maxNORM(double z)  // = maximum norm
{
  return std::fabs(z);
}

// abs(z) <= 1, actually only used for abs(z) < 1/2
double series_2F1(double a, double b, double c, double z)
{
  int n;
  int nMax;
  double r;
  double s;
  double sNew;
  double t;

  nMax = 2000;

  if (areEqual(a, 0.0) == 1 || areEqual(b, 0.0) == 1) return (1.0);
  if (areEqual(z, 0.0) == 1) return (1.0);

  n = 0;
  s = 1.0;
  t = 1.0;
  while (n < nMax)
  {
    n = n + 1;
    r = (a * b / c / static_cast<double>(n) * z);
    t = (t * r);
    sNew = (t + s);
    if (areEqual(sNew, s) == 1 && 3 < n && maxNORM(t) < std::numeric_limits<double>::epsilon())
      break;
    else
      s = (sNew);

    a = (a + 1.0);
    b = (b + 1.0);
    c = (c + 1.0);
  }

  return (s);
}

// Gosper's method
// http://www.math.utexas.edu/pipermail/maxima/2006/000126.html
// https://www.mapleprimes.com/posts/41570-c-code-for-the-hypergeometric
double hypergeometric2F1(double a, double b, double c, double z)
{
  int k;
  // double K;
  int kMax;
  double d_k;
  double e_k;
  double f_k;
  double d_k1;
  double e_k1;
  double f_k1, t;
  double xi;
  double lambda;
  double mu;
  double result;

  if (std::fabs(z) <= 0.5)
  {
    result = series_2F1(a, b, c, z);
    return (result);
  }

  kMax = 2000;
  xi = z / (z - 1.0);
  mu = c - a - b;
  d_k = 0.0;
  e_k = 1.0;
  f_k = 0.0;
  for (k = 0; k <= kMax; k++)
  {
    lambda = (static_cast<double>(k) + a) * (static_cast<double>(k) + b) / (static_cast<double>(k) + 1.0) /
             (2.0 * static_cast<double>(k) + c) / (2.0 * static_cast<double>(k) + c + 1.0);
    d_k1 = (lambda * z * (xi * d_k * (mu + static_cast<double>(k)) + e_k));
    e_k1 = (lambda * z * (-a * b * xi * d_k + (static_cast<double>(k) + c) * e_k));
    t = d_k * (-static_cast<double>(k) / (-1.0 + z) - (b * a - (mu + static_cast<double>(k)) * static_cast<double>(k)) /
                                                          (2.0 * static_cast<double>(k) + c) * xi) +
        e_k;
    f_k1 = (f_k + t);
    if (areEqual(f_k1, f_k) == 1 && maxNORM(t) <= std::numeric_limits<double>::epsilon() && 3 < k) break;
    d_k = (d_k1);
    e_k = (e_k1);
    f_k = (f_k1);
  }
  result = f_k;

  return (result);
}
