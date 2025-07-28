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

module bond_bond_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BondBondPotential::BondBondPotential(std::array<std::size_t, 3> identifiers, BondBondType type,
                                     std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBondBondParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondBondType::CVFF:
    case BondBondType::CFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      parameters[0] *= Units::KelvinToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string BondBondPotential::print() const
{
  switch (type)
  {
    case BondBondType::CVFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return std::format("{} - {} - {} : CVFF p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2]);
    case BondBondType::CFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return std::format("{} - {} - {} : CFF p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2]);
    default:
      std::unreachable();
  }
}

double BondBondPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC) const
{
  double3 dr_ab = posA - posB;
  double rr_ab = double3::dot(dr_ab, dr_ab);
  double r_ab = std::sqrt(rr_ab);

  double3 dr_cb = posC - posB;
  double rr_cb = double3::dot(dr_cb, dr_cb);
  double r_cb = std::sqrt(rr_cb);

  switch (type)
  {
    case BondBondType::CVFF:
    case BondBondType::CFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return parameters[0] * (r_ab - parameters[1]) * (r_cb - parameters[2]);
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondBondPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondBondPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondBondPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BondBondPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
