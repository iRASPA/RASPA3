module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#endif

module ring;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

Ring::Ring() {}

int Ring::floorDivision(int a, int b) { return int(std::floor(double(a) / double(b))); }

int Ring::modulo(int a, int b) { return a - b * int(std::floor(double(a) / double(b))); }

int Ring::greatestCommonDivisor(int arg1, int arg2)
{
  int a = arg1;
  int b = arg2;
  while (b != 0)
  {
    int stored_a = a;
    a = b;
    b = stored_a % b;
  }
  return std::abs(a);
}

std::tuple<int, int, int> Ring::extendedGreatestCommonDivisor(int a, int b)
{
  int ai = b;    // ai stands for: a with index i
  int aim1 = a;  // aim1 stands for: a with index i-1
  int bim1 = 0;
  int cim1 = 0;

  // We can accelerate the first step
  if (ai != 0)
  {
    // compute both quotient and remainder
    int q = aim1 / ai;
    int r = aim1 % ai;

    aim1 = ai;
    ai = r;
    bim1 = 0;  // before: bi = 0, bim1 = 1
    int bi = 1;
    cim1 = 1;  // before: ci = 1, cim1 = 0
    int ci = -q;
    // Now continue
    while (ai != 0)
    {
      // compute both quotient and remainder
      q = aim1 / ai;
      r = aim1 % ai;

      aim1 = ai;
      ai = r;

      int stored_bim1 = bim1;
      bim1 = bi;
      bi = stored_bim1 - q * bi;

      int stored_cim1 = cim1;
      cim1 = ci;
      ci = stored_cim1 - q * ci;
    }
  }
  else
  {
    bim1 = 1;
    cim1 = 0;
  }
  if (aim1 < 0)
  {
    // Make sure that the GCD is non-negative
    aim1 = -aim1;
    bim1 = -bim1;
    cim1 = -cim1;
  }

  return std::make_tuple(aim1, bim1, cim1);
}

std::pair<int, int> Ring::divisionModulo(int a, int b)
{
  // guard (b != 0) else {throw NumericalError.divisionByZero}
  if (b == 0)
  {
    return std::make_pair(0, 0);
  }
  int temp = int(std::floor(double(a) / double(b)));
  return std::make_pair(temp, a - b * temp);
}
