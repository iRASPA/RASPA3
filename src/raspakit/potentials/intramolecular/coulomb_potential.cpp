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
import units;
import double3;

CoulombPotential::CoulombPotential(std::array<std::size_t, 2> identifiers, CoulombType type,
                                   std::vector<double> vector_parameters, double scaling)
    : identifiers(identifiers), type(type), scaling(scaling)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfCoulombParameters); ++i)
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
      return std::format("{} - {} : COULOMB p_0={:g} [{}], p_1={:g} [{}]\n", identifiers[0], identifiers[1],
                         parameters[0], "e", 
                         parameters[1], "e");
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
      return Units::CoulombicConversionFactor * parameters[0] * parameters[1] / r;
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
