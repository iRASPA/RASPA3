module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <tuple>
#endif

module skasymmetricunit;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;

SKAsymmetricUnit::SKAsymmetricUnit() {}

bool SKAsymmetricUnit::contains(double3 point)
{
  return (isInsideRange(point.x, double(a.first) / 48.0, a.first % 2, double(a.second) / 48.0, a.second % 2) &&
          isInsideRange(point.y, double(b.first) / 48.0, b.first % 2, double(b.second) / 48.0, b.second % 2) &&
          isInsideRange(point.z, double(c.first) / 48.0, c.first % 2, double(c.second) / 48.0, c.second % 2));
}

bool SKAsymmetricUnit::isInsideRange(double point, double leftBoundary, int equalLeft, double rightBoundary,
                                     int equalRight)
{
  double centeredPoint = point - std::rint(point);

  if (equalLeft == 0 && equalRight == 0)
  {
    if (leftBoundary <= centeredPoint && centeredPoint <= rightBoundary)
    {
      return true;
    }
    if (leftBoundary <= (centeredPoint + 1.0) && (centeredPoint + 1.0) <= rightBoundary)
    {
      return true;
    }
  }

  if (equalLeft != 0 && equalRight == 0)
  {
    if (leftBoundary < centeredPoint && centeredPoint <= rightBoundary)
    {
      return true;
    }
    if (leftBoundary < (centeredPoint + 1.0) && (centeredPoint + 1.0) <= rightBoundary)
    {
      return true;
    }
  }

  if (equalLeft == 0 && equalRight != 0)
  {
    if (leftBoundary <= centeredPoint && centeredPoint < rightBoundary)
    {
      return true;
    }
    if (leftBoundary <= (centeredPoint + 1.0) && (centeredPoint + 1.0) < rightBoundary)
    {
      return true;
    }
  }

  if (equalLeft != 0 && equalRight != 0)
  {
    if (leftBoundary < centeredPoint && centeredPoint < rightBoundary)
    {
      return true;
    }
    if (leftBoundary < (centeredPoint + 1.0) && (centeredPoint + 1.0) < rightBoundary)
    {
      return true;
    }
  }

  return false;
}
