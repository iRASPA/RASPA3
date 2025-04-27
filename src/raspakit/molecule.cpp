module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <istream>
#include <map>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <utility>
#include <vector>
#endif

module molecule;

#ifndef USE_LEGACY_HEADERS
import <istream>;
import <ostream>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <algorithm>;
import <print>;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Molecule &atom)
{
  archive << atom.centerOfMassPosition;
  archive << atom.velocity;
  archive << atom.gradient;
  archive << atom.orientation;
  archive << atom.orientationMomentum;
  archive << atom.orientationGradient;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
};

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Molecule &atom)
{
  archive >> atom.centerOfMassPosition;
  archive >> atom.velocity;
  archive >> atom.gradient;
  archive >> atom.orientation;
  archive >> atom.orientationMomentum;
  archive >> atom.orientationGradient;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Molecule: Error in binary restart\n"));
  }
#endif

  return archive;
}
