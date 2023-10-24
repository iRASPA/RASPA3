module;

#include <cmath>

module double4;

import <fstream>;
import <complex>;

import archive;

void double4::normalise()
{
  double magnitude = sqrt((x * x) + (y * y) + (z * z) * (w * w));

  if (magnitude != 0)
  {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
    w /= magnitude;
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const double4 &vec)
{
  archive << vec.x << vec.y << vec.z << vec.w;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, double4 &vec)
{
  archive >> vec.x >> vec.y >> vec.z >> vec.w;
  return archive;
}

