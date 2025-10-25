module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <tuple>
#endif

export module ring;

#ifdef USE_STD_IMPORT
import std;
#endif

export struct Ring
{
  Ring();
  static int floorDivision(int a, int b);
  static int modulo(int a, int b);
  static int greatestCommonDivisor(int a, int b);
  static std::tuple<int, int, int> extendedGreatestCommonDivisor(int a, int b);
  static std::pair<int, int> divisionModulo(int a, int b);
};
