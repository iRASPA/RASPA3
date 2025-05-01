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

module enthalpy_of_adsorption;

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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnthalpyOfAdsorption &p)
{
  archive << p.size;
  archive << p.values;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnthalpyOfAdsorption &p)
{
  archive >> p.size;
  archive >> p.values;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyEnthalpy: Error in binary restart\n"));
  }
#endif

  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnthalpyOfAdsorptionTerms &p)
{
  archive << p.size;
  archive << p.swappableComponents;
  archive << p.totalEnergyTimesNumberOfMolecules;
  archive << p.numberOfMoleculesSquared;
  archive << p.numberOfMolecules;
  archive << p.temperature;
  archive << p.totalEnergy;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnthalpyOfAdsorptionTerms &p)
{
  archive >> p.size;
  archive >> p.swappableComponents;
  archive >> p.totalEnergyTimesNumberOfMolecules;
  archive >> p.numberOfMoleculesSquared;
  archive >> p.numberOfMolecules;
  archive >> p.temperature;
  archive >> p.totalEnergy;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("EnthalpyOfAdsorptionTerms: Error in binary restart\n"));
  }
#endif

  return archive;
}
