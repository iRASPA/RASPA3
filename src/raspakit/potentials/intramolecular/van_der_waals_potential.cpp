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

module van_der_waals_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

VanDerWaalsPotential::VanDerWaalsPotential(std::array<std::size_t, 2> identifiers, VanDerWaalsType type,
                                           std::vector<double> vector_parameters, double shift, double scaling)
    : identifiers(identifiers), type(type), shift(shift), scaling(scaling)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfVanDerWaalsParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      parameters[0] *= Units::KelvinToEnergy;
      break;
  }
}

std::string VanDerWaalsPotential::print() const
{
  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      return std::format("{} - {} : LENNARD_JONES p_0/k_B={:g} [K], p_1={:g} [Å], shift/k_B={:g}\n", identifiers[0],
                         identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         shift * Units::EnergyToKelvin);
  }
}

double VanDerWaalsPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  double temp;
  double rri;

  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);

  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      rri = (parameters[1] * parameters[1]) / rr;
      temp = rri * rri * rri;
      return 4.0 * parameters[0] * (temp * (temp - 1.0)) - shift;
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VanDerWaalsPotential &b)
{
  archive << b.versionNumber;

  archive << b.type;
  archive << b.identifiers;
  archive << b.shift;
  archive << b.scaling;
  archive << b.parameters;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VanDerWaalsPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'VanDerWaalsPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.type;
  archive >> b.identifiers;
  archive >> b.shift;
  archive >> b.scaling;
  archive >> b.parameters;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("VanDerWaalsPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
