module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <tuple>
#include <utility>
#endif

export module skasymmetricunit;

#ifndef USE_LEGACY_HEADERS
import <tuple>;
import <utility>;
#endif

import int3;
import int3x3;
import double3;
import double3x3;

export struct SKAsymmetricUnit
{
  std::pair<int, int> a;
  std::pair<int, int> b;
  std::pair<int, int> c;

  SKAsymmetricUnit();
  SKAsymmetricUnit(std::pair<int, int> p1, std::pair<int, int> p2, std::pair<int, int> p3) : a(p1), b(p2), c(p3) {};

  SKAsymmetricUnit(const SKAsymmetricUnit& t)
  {
    this->a = t.a;
    this->b = t.b;
    this->c = t.c;
  }

  SKAsymmetricUnit& operator=(const SKAsymmetricUnit& other)
  {
    this->a = other.a;
    this->b = other.b;
    this->c = other.c;
    return *this;
  }

  bool contains(double3 point);
  bool isInsideRange(double point, double leftBoundary, int equalLeft, double rightBoundary, int equalRight);
};
