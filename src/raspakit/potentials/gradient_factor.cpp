module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <vector>
#endif

module gradient_factor;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <print>;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const GradientFactor &e)
{
  archive << e.energy;
  archive << e.dUdlambda;
  archive << e.gradientFactor;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, GradientFactor &e)
{
  archive >> e.energy;
  archive >> e.dUdlambda;
  archive >> e.gradientFactor;

  return archive;
}
