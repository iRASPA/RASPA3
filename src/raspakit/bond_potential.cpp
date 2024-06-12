module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <exception>
#include <source_location>
#include <complex>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <algorithm>
#include <print>
#endif
  
module bond_potential;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondPotential &b)
{
  archive << b.versionNumber;

  archive << b.bondType;
  archive << b.bondIds;
  archive << b.parameters;

  return archive;
}


Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondPotential &b)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > b.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.bondType;
  archive >> b.bondIds;
  archive >> b.parameters;

  return archive;
}

