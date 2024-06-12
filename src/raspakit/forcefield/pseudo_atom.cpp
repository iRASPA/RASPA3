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
#include <algorithm>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#include <type_traits>
#include <iterator>
#include <functional>
#include <print>
#endif

module pseudo_atom;

#ifndef USE_LEGACY_HEADERS
import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <vector>;
import <array>;
import <map>;
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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PseudoAtom &a)
{
  archive << a.versionNumber;
  archive << a.name;
  archive << a.mass;
  archive << a.charge;
  archive << a.atomicNumber;
  archive << a.printToPDB;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PseudoAtom &a)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > a.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PseudoAtom' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> a.name;
  archive >> a.mass;
  archive >> a.charge;
  archive >> a.atomicNumber;
  archive >> a.printToPDB;

  return archive;
}

bool PseudoAtom::operator==(const PseudoAtom &other) const
{
  return (versionNumber == other.versionNumber && name == other.name && mass == other.mass && charge == other.charge &&
          atomicNumber == other.atomicNumber && printToPDB == other.printToPDB && source == other.source);
}
