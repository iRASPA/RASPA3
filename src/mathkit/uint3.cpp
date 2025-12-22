module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <fstream>
#include <map>
#include <vector>
#endif

module uint3;

#ifdef USE_STD_IMPORT
import std;
#endif

import ring;
import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const uint3 &vec)
{
  archive << vec.x << vec.y << vec.z;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, uint3 &vec)
{
  archive >> vec.x >> vec.y >> vec.z;
  return archive;
}
