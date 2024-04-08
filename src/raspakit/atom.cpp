module;

#ifdef USE_LEGACY_HEADERS
#include <istream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
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
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
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

  return archive;
}


