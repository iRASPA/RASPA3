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

module inversion_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

InversionBendPotential::InversionBendPotential(std::array<std::size_t, 4> identifiers, InversionBendType type,
                                               std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfInversionBendParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::Harmonic2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case InversionBendType::HarmonicCosine:
    case InversionBendType::HarmonicCosine2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] = std::cos(parameters[1] * Units::DegreesToRadians);
      break;
    case InversionBendType::Planar:
    case InversionBendType::Planar2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case InversionBendType::MM3:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      parameters[0] *= 0.021914 * Units::KCalPerMolToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    default:
      std::unreachable();
  }
}

std::string InversionBendPotential::print() const
{
  switch (type)
  {
    case InversionBendType::Harmonic:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case InversionBendType::HarmonicCosine:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC_COSINE p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, std::acos(parameters[1]) * Units::RadiansToDegrees);
    case InversionBendType::Planar:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      return std::format("{} - {} - {} - {} : PLANAR p_0/k_B={:g} [K]\n", identifiers[0], identifiers[1],
                         identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin);
    case InversionBendType::Harmonic2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC2 p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case InversionBendType::HarmonicCosine2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC_COSINE2 p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, std::acos(parameters[1]) * Units::RadiansToDegrees);
    case InversionBendType::Planar2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      return std::format("{} - {} - {} - {} : PLANAR2 p_0/k_B={:g} [K]\n", identifiers[0], identifiers[1],
                         identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin);
    case InversionBendType::MM3:
      // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : MM3 p_0/k_B={:g} [mdyne A/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    default:
      std::unreachable();
  }
}

double InversionBendPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                               const double3 &posD) const
{
  double c, chi;
  double temp, temp2;

  double3 Rab = posA - posB;
  double rab2 = double3::dot(Rab, Rab);
  double rrab = std::sqrt(rab2);

  double3 Rbc = posC - posB;
  double3 Rbd = posD - posB;
  double3 Rcd = posD - posC;
  double3 Rad = posD - posA;

  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::HarmonicCosine:
    case InversionBendType::Planar:
      // w is a vector perpendicular to the B-C-D plane
      // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
      temp = double3::dot(Rbc, Rbd);
      c = double3::dot(Rbc, Rbc) * double3::dot(Rbd, Rbd) - temp * temp;
      break;
    case InversionBendType::Harmonic2:
    case InversionBendType::HarmonicCosine2:
    case InversionBendType::Planar2:
    case InversionBendType::MM3:
      // w is a vector perpendicular to the A-C-D plane
      // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
      temp = double3::dot(Rad, Rcd);
      c = double3::dot(Rcd, Rcd) * double3::dot(Rad, Rad) - temp * temp;
      break;
    default:
      std::unreachable();
  }

  double e = double3::dot(Rab, double3::cross(Rbd, Rbc));
  double cos_chi = std::sqrt(double3::dot(Rab, Rab) - e * e / c) / rrab;
  cos_chi = std::clamp(cos_chi, -1.0, 1.0);

  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::Harmonic2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      chi = std::acos(cos_chi);
      temp = chi - parameters[1];
      return 0.5 * parameters[0] * temp * temp;
    case InversionBendType::HarmonicCosine:
    case InversionBendType::HarmonicCosine2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      temp = cos_chi - parameters[1];
      return 0.5 * parameters[0] * temp * temp;
      break;
    case InversionBendType::Planar:
    case InversionBendType::Planar2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      return parameters[0] * (1.0 - cos_chi);
    case InversionBendType::MM3:
      // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      chi = std::acos(cos_chi);
      temp = (chi - parameters[1]) * Units::RadiansToDegrees;
      temp2 = temp * temp;
      return parameters[0] * temp2 *
             (1.0 - 0.014 * temp + 5.6e-5 * temp2 - 7.0e-7 * temp * temp2 + 2.2e-8 * temp2 * temp2);
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const InversionBendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, InversionBendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'InversionBendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("InversionBendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
