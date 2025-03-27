module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <fstream>
#include <map>
#include <vector>
#endif

module int3;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
#endif

import ring;
import archive;

int3 int3::greatestCommonDivisor(int3 a, int b)
{
  return int3(Ring::greatestCommonDivisor(a.x, b), Ring::greatestCommonDivisor(a.y, b),
              Ring::greatestCommonDivisor(a.z, b));
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const int3 &vec)
{
  archive << vec.x << vec.y << vec.z;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, int3 &vec)
{
  archive >> vec.x >> vec.y >> vec.z;
  return archive;
}
