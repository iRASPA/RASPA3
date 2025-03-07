module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <utility>
#include <vector>
#endif

module pressure_range;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <utility>;
import <print>;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PressureRange &r)
{
  archive << r.pressureStart;
  archive << r.pressureEnd;
  archive << r.numberOfPoints;
  archive << r.scale;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PressureRange &r)
{
  archive >> r.pressureStart;
  archive >> r.pressureEnd;
  archive >> r.numberOfPoints;
  archive >> r.scale;

  return archive;
}
