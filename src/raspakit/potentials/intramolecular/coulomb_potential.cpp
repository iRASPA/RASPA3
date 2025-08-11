module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <utility>
#include <vector>
#endif

module coulomb_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

CoulombPotential::CoulombPotential(std::array<std::size_t, 2> identifiers, CoulombType type,
                                   std::vector<double> vector_parameters, double scaling)
    : identifiers(identifiers), type(type), scaling(scaling)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfCoulombParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case CoulombType::Coulomb:
      parameters[0] *= Units::KelvinToEnergy;
      break;
  }
}

std::string CoulombPotential::print() const
{
  switch (type)
  {
    case CoulombType::Coulomb:
      return std::format("{} - {} : COULOMB p_0/k_B={:g} [K], p_1={:g} [Å]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
  }
}

double CoulombPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  switch (type)
  {
    case CoulombType::Coulomb:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      return 1.0 / r;
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CoulombPotential &b)
{
  archive << b.versionNumber;

  archive << b.type;
  archive << b.identifiers;
  archive << b.parameters;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CoulombPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'CoulombPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.type;
  archive >> b.identifiers;
  archive >> b.parameters;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("CoulombPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
