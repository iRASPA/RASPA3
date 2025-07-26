module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <vector>
#endif

module triquintic_derivative_factor;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Potentials::TriquinticDerivativeFactor &e)
{
  archive << e.energy;
  archive << e.firstDerivativeFactor;
  archive << e.secondDerivativeFactor;
  archive << e.thirdDerivativeFactor;
  archive << e.fourthDerivativeFactor;
  archive << e.fifthDerivativeFactor;
  archive << e.sixthDerivativeFactor;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Potentials::TriquinticDerivativeFactor &e)
{
  archive >> e.energy;
  archive >> e.firstDerivativeFactor;
  archive >> e.secondDerivativeFactor;
  archive >> e.thirdDerivativeFactor;
  archive >> e.fourthDerivativeFactor;
  archive >> e.fifthDerivativeFactor;
  archive >> e.sixthDerivativeFactor;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Potentials::TriquinticDerivativeFactor: Error in binary restart\n"));
  }
#endif

  return archive;
}
