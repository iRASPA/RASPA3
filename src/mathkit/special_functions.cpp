module;

module special_functions;

import std;

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
// the direct series sum_n (a)_n (b)_n / (c)_n x^n / n!, which converges for |x| < 1 (and for a terminating series)
double hypergeometric(double a, double b, double c, double x)
{
  const int nMax = 2000;
  double term = 1.0;
  double value = 1.0;

  for (int n = 0; n < nMax; ++n)
  {
    // the ratio of successive terms tends to x, but for large a or b the early terms grow before they decay
    const double ratio = (a + n) * (b + n) * x / ((c + n) * (n + 1.0));
    term *= ratio;
    if (term == 0.0) break;
    value += term;
    if (std::fabs(ratio) < 1.0 && std::fabs(term) <= std::numeric_limits<double>::epsilon() * std::fabs(value)) break;
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

namespace
{
bool isNonPositiveInteger(double x) { return x <= 0.0 && x == std::floor(x); }

// 1 / Gamma(x), zero at the poles of Gamma
double reciprocalGamma(double x) { return isNonPositiveInteger(x) ? 0.0 : 1.0 / std::tgamma(x); }

// psi(x) = Gamma'(x) / Gamma(x): reflection for x <= 0, the recurrence up to x >= 10 and the asymptotic series
double digamma(double x)
{
  if (x <= 0.0) return digamma(1.0 - x) - std::numbers::pi / std::tan(std::numbers::pi * x);
  double result = 0.0;
  while (x < 10.0)
  {
    result -= 1.0 / x;
    x += 1.0;
  }
  const double f = 1.0 / (x * x);
  return result + std::log(x) - 0.5 / x -
         f * (1.0 / 12.0 - f * (1.0 / 120.0 - f * (1.0 / 252.0 - f * (1.0 / 240.0 - f / 132.0))));
}

// the logarithmic sum of Abramowitz & Stegun 15.3.10-12 in w = 1 - z,
//   sum_n (a + s)_n (b + s)_n / (n! (n + m)!) w^n [ln w - psi(n + 1) - psi(n + m + 1) + psi(a + s + n) + psi(b + s +
//   n)]
double logarithmicSum(double a, double b, int m, double w)
{
  const double logW = std::log(w);
  double term = 1.0 / std::tgamma(static_cast<double>(m) + 1.0);
  double sum{};
  for (int n = 0; n < 100000; ++n)
  {
    if (n > 0) term *= (a + n - 1) * (b + n - 1) / (static_cast<double>(n) * static_cast<double>(n + m)) * w;
    const double value = term * (logW - digamma(n + 1.0) - digamma(n + m + 1.0) + digamma(a + n) + digamma(b + n));
    sum += value;
    if (n > 3 && std::fabs(value) <= std::numeric_limits<double>::epsilon() * std::fabs(sum)) break;
  }
  return sum;
}

// the finite sum sum_{n < m} (a)_n (b)_n / (n! (1 - m)_n) w^n of A&S 15.3.11-12
double finiteSum(double a, double b, int m, double w)
{
  double term = 1.0;
  double sum{};
  for (int n = 0; n < m; ++n)
  {
    if (n > 0) term *= (a + n - 1) * (b + n - 1) / (static_cast<double>(n) * (1.0 - m + n - 1)) * w;
    sum += term;
  }
  return sum;
}

// 2F1(a, b; c; z) for 1/2 < z < 1 from the connection formula in 1 - z (A&S 15.3.6), and its logarithmic forms
// (15.3.10-12) when c - a - b is an integer (within 1e-8, where the two terms of 15.3.6 cancel to the rounding of
// Gamma(c - a - b))
double hypergeometric2F1NearOne(double a, double b, double c, double z)
{
  const double w = 1.0 - z;
  const double m = c - a - b;
  const double integerM = std::round(m);
  if (std::fabs(m - integerM) > 1.0e-8)
  {
    return std::tgamma(c) * std::tgamma(m) * reciprocalGamma(c - a) * reciprocalGamma(c - b) *
               hypergeometric(a, b, 1.0 - m, w) +
           std::pow(w, m) * std::tgamma(c) * std::tgamma(-m) * reciprocalGamma(a) * reciprocalGamma(b) *
               hypergeometric(c - a, c - b, m + 1.0, w);
  }
  if (integerM >= 0.0)
  {
    const int k = static_cast<int>(integerM);
    const double finite = k > 0 ? std::tgamma(static_cast<double>(k)) * std::tgamma(c) * reciprocalGamma(a + k) *
                                      reciprocalGamma(b + k) * finiteSum(a, b, k, w)
                                : 0.0;
    return finite - std::pow(-w, k) * std::tgamma(c) * reciprocalGamma(a) * reciprocalGamma(b) *
                        logarithmicSum(a + k, b + k, k, w);
  }
  const int k = static_cast<int>(-integerM);
  const double finite = std::tgamma(static_cast<double>(k)) * std::tgamma(c) * reciprocalGamma(a) * reciprocalGamma(b) *
                        std::pow(w, -k) * finiteSum(a - k, b - k, k, w);
  return finite - (k % 2 == 0 ? 1.0 : -1.0) * std::tgamma(c) * reciprocalGamma(a - k) * reciprocalGamma(b - k) *
                      logarithmicSum(a, b, k, w);
}
}  // namespace

// 2F1(a, b; c; z) for real z <= 1: the power series for |z| <= 1/2 (and for a terminating series), Pfaff's
// transformation 2F1(a, b; c; z) = (1 - z)^(-a) 2F1(a, c - b; c; z / (z - 1)) for z < -1/2, and the connection
// formula in 1 - z for 1/2 < z < 1. (Gosper's recurrence, used before for |z| > 1/2, diverges for
// 2 sqrt(2) - 2 < z < 1 and z <= -4.)
double hypergeometric2F1(double a, double b, double c, double z)
{
  if (std::fabs(z) <= 0.5)
  {
    return series_2F1(a, b, c, z);
  }
  if (isNonPositiveInteger(a) || isNonPositiveInteger(b))
  {
    // a polynomial
    return hypergeometric(a, b, c, z);
  }

  const double mu = c - a - b;
  if (z == 1.0)
  {
    // the connection formula below takes log(1 - z); at z = 1 Gauss's theorem gives Gamma(c) Gamma(c-a-b) / (Gamma(c-a)
    // Gamma(c-b)) for c - a - b > 0 and the function is unbounded for c - a - b <= 0, with the sign of its prefactor
    // Gamma(c) Gamma(a+b-c) / (Gamma(a) Gamma(b)) (of -Gamma(c) log(1-z) / (Gamma(a) Gamma(b)) for c - a - b = 0)
    if (mu > 0.0) return std::tgamma(c) * std::tgamma(mu) / (std::tgamma(c - a) * std::tgamma(c - b));
    const double prefactor = std::tgamma(c) * (mu < 0.0 ? std::tgamma(-mu) : 1.0) / (std::tgamma(a) * std::tgamma(b));
    return std::copysign(std::numeric_limits<double>::infinity(), prefactor);
  }
  if (z > 1.0) return std::numeric_limits<double>::quiet_NaN();
  if (z < 0.0) return std::pow(1.0 - z, -a) * hypergeometric2F1(a, c - b, c, z / (z - 1.0));
  return hypergeometric2F1NearOne(a, b, c, z);
}
