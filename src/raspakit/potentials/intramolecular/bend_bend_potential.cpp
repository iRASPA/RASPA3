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
#include <numbers>
#include <print>
#include <source_location>
#include <tuple>
#include <utility>
#include <vector>
#endif

module bend_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BendBendPotential::BendBendPotential(std::array<std::size_t, 4> identifiers, BendBendType type,
                                     std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfBendBendParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BendBendType::CVFF:
    case BendBendType::CFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      parameters[0] *= Units::EnergyToKelvin;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::DegreesToRadians;
      break;
    case BendBendType::MM3:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      parameters[0] *= 0.02191418 * Units::KCalPerMolToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::DegreesToRadians;
      break;
    default:
      std::unreachable();
  }
}

std::string BendBendPotential::print() const
{
  switch (type)
  {
    case BendBendType::CVFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : CFF p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::RadiansToDegrees,
                         parameters[2] * Units::RadiansToDegrees);
    case BendBendType::CFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : CFF p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::RadiansToDegrees,
                         parameters[2] * Units::RadiansToDegrees);
    case BendBendType::MM3:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : MM3 p_0={:g} [mdyne A/rad^2], p_1={:g} [degrees], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKCalPerMol / 0.02191418,
                         parameters[1] * Units::RadiansToDegrees, parameters[2] * Units::RadiansToDegrees);
    default:
      std::unreachable();
  }
}

double BendBendPotential::calculateEnergy([[maybe_unused]] const double3 &posA, [[maybe_unused]] const double3 &posB,
                                          [[maybe_unused]] const double3 &posC,
                                          [[maybe_unused]] const double3 &posD) const
{
  switch (type)
  {
    case BendBendType::CVFF:
      return 0.0;
    case BendBendType::CFF:
      return 0.0;
    case BendBendType::MM3:
      return 0.0;
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendBendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendBendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BendBendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BendBendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
