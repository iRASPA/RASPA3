module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <optional>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
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
import <print>;
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
