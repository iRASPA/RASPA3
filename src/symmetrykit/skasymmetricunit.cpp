module;

module skasymmetricunit;

import std;

import int3;
import int3x3;
import double3;
import double3x3;

SKAsymmetricUnit::SKAsymmetricUnit() {}

std::optional<double> SKAsymmetricUnit::coordinateInside(double point, std::pair<int, int> range,
                                                         double tolerance) const
{
  // the bounds are in units of 1/48; an odd value marks an open boundary one unit above the even bound it encodes
  // (49 is the open boundary at 48/48 = 1), so the parity is both the flag and an offset to strip
  const bool openLeft = range.first % 2 != 0;
  const bool openRight = range.second % 2 != 0;
  const double left = double(range.first - std::abs(range.first % 2)) / 48.0;
  const double right = double(range.second - std::abs(range.second % 2)) / 48.0;

  // the range may lie anywhere within [-1, 1], so the point is tested in all three of its images that can fall there
  const double centred = point - std::rint(point);
  for (const double shift : {-1.0, 0.0, 1.0})
  {
    const double x = centred + shift;
    const bool aboveLeft = openLeft ? left < x : left - tolerance <= x;
    const bool belowRight = openRight ? x < right : x <= right + tolerance;
    if (aboveLeft && belowRight) return x;
  }
  return std::nullopt;
}

std::optional<double3> SKAsymmetricUnit::imageInside(double3 point, double tolerance) const
{
  const std::optional<double> x = coordinateInside(point.x, a, tolerance);
  const std::optional<double> y = coordinateInside(point.y, b, tolerance);
  const std::optional<double> z = coordinateInside(point.z, c, tolerance);
  if (!x || !y || !z) return std::nullopt;
  return double3(*x, *y, *z);
}
