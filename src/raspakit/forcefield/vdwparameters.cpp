module;

#ifdef USE_LEGACY_HEADERS
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <string>
#include <string_view>
#include <optional>
#include <numbers>
#include <algorithm>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <exception>
#include <source_location>
#include <complex>
#include <type_traits>
#include <iterator>
#include <functional>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module vdwparameters;

#ifndef USE_LEGACY_HEADERS
import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <vector>;
import <array>;
import <map>;
import <cmath>;
import <string>;
import <string_view>;
import <optional>;
import <numbers>;
import <algorithm>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <type_traits>;
import <iterator>;
import <functional>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;

bool VDWParameters::operator==(const VDWParameters &other) const
{
  return (parameters == other.parameters && shift == other.shift &&
          tailCorrectionEnergy == other.tailCorrectionEnergy && type == other.type);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p)
{
  archive << p.parameters;
  archive << p.shift;
  archive << p.tailCorrectionEnergy;
  archive << p.type;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p)
{
  archive >> p.parameters;
  archive >> p.shift;
  archive >> p.tailCorrectionEnergy;
  archive >> p.type;

  return archive;
}

