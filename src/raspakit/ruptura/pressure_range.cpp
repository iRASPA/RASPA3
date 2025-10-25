module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <utility>
#include <vector>
#endif

module pressure_range;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PressureRange &r)
{
  archive << r.pressureStart;
  archive << r.pressureEnd;
  archive << r.numberOfPoints;
  archive << r.scale;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PressureRange &r)
{
  archive >> r.pressureStart;
  archive >> r.pressureEnd;
  archive >> r.numberOfPoints;
  archive >> r.scale;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PressureRange: Error in binary restart\n"));
  }
#endif

  return archive;
}
