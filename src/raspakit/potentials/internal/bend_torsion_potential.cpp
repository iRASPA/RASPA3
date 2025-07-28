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

module bend_torsion_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BendTorsionPotential::BendTorsionPotential(std::array<std::size_t, 4> identifiers, BendTorsionType type,
                                           std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBendTorsionParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BendTorsionType::Smoothed:
      break;
    case BendTorsionType::SmoothedThreeCosine:
      break;
    case BendTorsionType::Nicholas:
      break;
    case BendTorsionType::CFF:
      break;
    case BendTorsionType::SmoothedCFF:
      break;
    case BendTorsionType::SmoothedCFF2:
      break;
    case BendTorsionType::SmoothedCFF3:
      break;
    case BendTorsionType::CVFF:
      break;
    default:
      std::unreachable();
  }
}

std::string BendTorsionPotential::print() const
{
  switch (type)
  {
    case BendTorsionType::Smoothed:
      return "not implemente yet";
    case BendTorsionType::SmoothedThreeCosine:
      return "not implemente yet";
    case BendTorsionType::Nicholas:
      return "not implemente yet";
    case BendTorsionType::CFF:
      return "not implemente yet";
    case BendTorsionType::SmoothedCFF:
      return "not implemente yet";
    case BendTorsionType::SmoothedCFF2:
      return "not implemente yet";
    case BendTorsionType::SmoothedCFF3:
      return "not implemente yet";
    case BendTorsionType::CVFF:
      return "not implemente yet";
    default:
      std::unreachable();
  }
}

double BendTorsionPotential::calculateEnergy([[maybe_unused]] const double3 &posA, [[maybe_unused]] const double3 &posB,
                                             [[maybe_unused]] const double3 &posC,
                                             [[maybe_unused]] const double3 &posD) const
{
  switch (type)
  {
    case BendTorsionType::Smoothed:
      return 0.0;
    case BendTorsionType::SmoothedThreeCosine:
      return 0.0;
    case BendTorsionType::Nicholas:
      return 0.0;
    case BendTorsionType::CFF:
      return 0.0;
    case BendTorsionType::SmoothedCFF:
      return 0.0;
    case BendTorsionType::SmoothedCFF2:
      return 0.0;
    case BendTorsionType::SmoothedCFF3:
      return 0.0;
    case BendTorsionType::CVFF:
      return 0.0;
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendTorsionPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendTorsionPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BendTorsionPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BendTorsionPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
