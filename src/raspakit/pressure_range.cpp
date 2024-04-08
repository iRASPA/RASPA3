module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <exception>
#include <source_location>
#include <complex>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module pressure_range;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
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

