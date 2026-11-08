module;

module skasymmetricunit;

import std;

import int3;
import int3x3;
import double3;
import double3x3;

SKAsymmetricUnit::SKAsymmetricUnit() {}

bool SKAsymmetricUnit::contains(double3 point)
{
  // the bounds are in units of 1/48; an odd value marks an open boundary one unit above the even bound it encodes
  // (49 is the open boundary at 48/48 = 1), so the parity is both the flag and an offset to strip
  return (isInsideRange(point.x, double(a.first - std::abs(a.first % 2)) / 48.0, a.first % 2,
                        double(a.second - std::abs(a.second % 2)) / 48.0, a.second % 2) &&
          isInsideRange(point.y, double(b.first - std::abs(b.first % 2)) / 48.0, b.first % 2,
                        double(b.second - std::abs(b.second % 2)) / 48.0, b.second % 2) &&
          isInsideRange(point.z, double(c.first - std::abs(c.first % 2)) / 48.0, c.first % 2,
                        double(c.second - std::abs(c.second % 2)) / 48.0, c.second % 2));
}

bool SKAsymmetricUnit::isInsideRange(double point, double leftBoundary, int equalLeft, double rightBoundary,
                                     int equalRight)
{
  // the range may lie anywhere within [-1, 1] (42 settings have a lower bound below -1/2), so the point is tested
  // in all three of its images that can fall there; the image below used to be left out
  double centeredPoint = point - std::rint(point);
  for (const double shift : {-1.0, 0.0, 1.0})
  {
    const double x = centeredPoint + shift;
    const bool aboveLeft = equalLeft == 0 ? leftBoundary <= x : leftBoundary < x;
    const bool belowRight = equalRight == 0 ? x <= rightBoundary : x < rightBoundary;
    if (aboveLeft && belowRight) return true;
  }
  return false;
}
