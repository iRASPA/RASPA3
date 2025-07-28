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

module bond_torsion_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BondTorsionPotential::BondTorsionPotential(std::array<std::size_t, 4> identifiers, BondTorsionType type,
                                           std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBondTorsionParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondTorsionType::MM3:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      parameters[0] *= Units::KCalPerMolToEnergy;
      parameters[1] *= Units::KCalPerMolToEnergy;
      parameters[2] *= Units::KCalPerMolToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string BondTorsionPotential::print() const
{
  switch (type)
  {
    case BondTorsionType::MM3:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      return std::format(
          "{} - {} - {} - {} : MM3 p_0={:g} [kcal/A mole], p_1={:g} [kcal/A mole], p_2={:g} [kcal/A mole], "
          "p_3/k_B={:g} [A]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKCalPerMol,
          parameters[1] * Units::EnergyToKCalPerMol, parameters[2] * Units::EnergyToKCalPerMol, parameters[3]);
      break;
    default:
      std::unreachable();
  }
}

double BondTorsionPotential::calculateEnergy([[maybe_unused]] const double3 &posA, [[maybe_unused]] const double3 &posB,
                                             [[maybe_unused]] const double3 &posC,
                                             [[maybe_unused]] const double3 &posD) const
{
  switch (type)
  {
    case BondTorsionType::MM3:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      return 0.0;
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondTorsionPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondTorsionPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondTorsionPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BondTorsionPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
