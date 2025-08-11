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

module out_of_plane_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

OutOfPlaneBendPotential::OutOfPlaneBendPotential(std::array<std::size_t, 4> identifiers, OutOfPlaneBendType type,
                                                 std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < parameters.size(); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case OutOfPlaneBendType::Harmonic:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    default:
      std::unreachable();
  }
}

std::string OutOfPlaneBendPotential::print() const
{
  switch (type)
  {
    case OutOfPlaneBendType::Harmonic:
      return std::string("not implemented yet");
      break;
    default:
      std::unreachable();
  }
}

double OutOfPlaneBendPotential::calculateEnergy([[maybe_unused]] const double3 &posA,
                                                [[maybe_unused]] const double3 &posB,
                                                [[maybe_unused]] const double3 &posC,
                                                [[maybe_unused]] const double3 &posD) const
{
  switch (type)
  {
    case OutOfPlaneBendType::Harmonic:
      return 0.0;
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const OutOfPlaneBendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, OutOfPlaneBendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'OutOfPlaneBendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("OutOfPlaneBendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
