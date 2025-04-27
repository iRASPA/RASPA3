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

module atom;

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
import stringutils;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Atom &atom)
{
  archive << atom.position;
  archive << atom.velocity;
  archive << atom.gradient;
  archive << atom.charge;
  archive << atom.scalingVDW;
  archive << atom.scalingCoulomb;
  archive << atom.moleculeId;
  archive << atom.type;
  archive << atom.componentId;
  archive << atom.groupId;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
};

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Atom &atom)
{
  archive >> atom.position;
  archive >> atom.velocity;
  archive >> atom.gradient;
  archive >> atom.charge;
  archive >> atom.scalingVDW;
  archive >> atom.scalingCoulomb;
  archive >> atom.moleculeId;
  archive >> atom.type;
  archive >> atom.componentId;
  archive >> atom.groupId;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Atom: Error in binary restart\n"));
  }
#endif

  return archive;
}
