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

module bond_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BondBendPotential::BondBendPotential(std::array<std::size_t, 4> identifiers, BondBendType type,
                                     std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBondBendParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondBendType::CVFF:
    case BondBendType::CFF:
      // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
      // =========================================
      // p_0     [degrees]
      // p_1/k_B [K/A/rad]
      // p_2     [A]
      // p_3/k_B [K/A/rad]
      // p_4     [A]
      parameters[0] *= Units::DegreesToRadians;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case BondBendType::MM3:
      // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
      // =====================================
      // p_0     [mdyne/rad]
      // p_1     [A]
      // p_2     [A]
      // p_3     [degrees]
      parameters[0] *= Units::KCalPerMolToEnergy;
      parameters[3] *= Units::DegreesToRadians;
      break;
    case BondBendType::TruncatedHarmonic:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
      // ================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case BondBendType::ScreenedHarmonic:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case BondBendType::ScreenedVessal:
      // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
      // ============================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case BondBendType::TruncatedVessal:
      // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
      //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
      // ============================================================================
      // p_0/k_B [K/rad^(4+p_2)]
      // p_1     [degrees]
      // p_2     [-]
      // p_3     [A]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
  }
}

std::string BondBendPotential::print() const
{
  switch (type)
  {
    case BondBendType::CVFF:
      // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
      // =========================================
      // p_0     [degrees]
      // p_1/k_B [K/A/rad]
      // p_2     [A]
      // p_3/k_B [K/A/rad]
      // p_4     [A]
      return std::format(
          "{} - {} - {} : CVFF p_0={:g} [degrees], p_1/k_B={:g} [K/A/rad], p_2={:g} [A], p_3/k_B={:g} [K/A/rad], "
          "p_4={:g} [A]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::RadiansToDegrees,
          parameters[1] * Units::EnergyToKelvin, parameters[2], parameters[3] * Units::EnergyToKelvin, parameters[4]);
    case BondBendType::CFF:
      // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
      // =========================================
      // p_0     [degrees]
      // p_1/k_B [K/A/rad]
      // p_2     [A]
      // p_3/k_B [K/A/rad]
      // p_4     [A]
      return std::format(
          "{} - {} - {} : CVFF p_0={:g} [degrees], p_1/k_B={:g} [K/A/rad], p_2={:g} [A], p_3/k_B={:g} [K/A/rad], "
          "p_4={:g} [A]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::RadiansToDegrees,
          parameters[1] * Units::EnergyToKelvin, parameters[2], parameters[3] * Units::EnergyToKelvin, parameters[4]);
    case BondBendType::MM3:
      // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
      // =====================================
      // p_0     [mdyne/rad]
      // p_1     [A]
      // p_2     [A]
      // p_3     [degrees]
      return std::format("{} - {} - {} : MM3 p_0={:g} [mdyne/rad], p_1={:g} [A], p_2={:g} [A], p_3={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2],
                         parameters[0] / (2.51117 * Units::KCalPerMolToEnergy), parameters[1], parameters[2],
                         parameters[3] * Units::RadiansToDegrees);
    case BondBendType::TruncatedHarmonic:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
      // ================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      return std::format("{} - {} - {} : TRUNCATED_HARMONIC p_0/k_B={:g} [K/A/rad], p_2={:g} [degrees], p_3={:g} [A]\n",
                         identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees, parameters[2]);
    case BondBendType::ScreenedHarmonic:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      return std::format(
          "{} - {} - {} : SCREENED_HARMONIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2={:g} [A], p_3={:g} [A]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2], parameters[3]);
      break;
    case BondBendType::ScreenedVessal:
      // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
      // ============================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      return std::format(
          "{} - {} - {} : SCREENED_VESSAL p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2={:g} [A], p_3={:g} [A]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2], parameters[3]);
    case BondBendType::TruncatedVessal:
      // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
      //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
      // ============================================================================
      // p_0/k_B [K/rad^(4+p_2)]
      // p_1     [degrees]
      // p_2     [-]
      // p_3     [A]
      return std::format(
          "{} - {} - {} : TRUNCATED_VESSAL p_0/k_B={:g} [K/rad^(4+p_2)], p_1={:g} [degrees], p_2={:g} [-], p_3={:g} "
          "[A]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2], parameters[3]);
  }
}

double BondBendPotential::calculateEnergy([[maybe_unused]] const double3 &posA, [[maybe_unused]] const double3 &posB,
                                          [[maybe_unused]] const double3 &posC,
                                          [[maybe_unused]] const double3 &posD) const
{
  switch (type)
  {
    case BondBendType::CVFF:
      return 0.0;
    case BondBendType::CFF:
      return 0.0;
    case BondBendType::MM3:
      return 0.0;
    case BondBendType::TruncatedHarmonic:
      return 0.0;
    case BondBendType::ScreenedHarmonic:
      return 0.0;
    case BondBendType::ScreenedVessal:
      return 0.0;
    case BondBendType::TruncatedVessal:
      return 0.0;
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondBendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondBendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondBendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BondBendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
