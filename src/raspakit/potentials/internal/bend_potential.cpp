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
#include <optional>
#endif

module bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BendPotential::BendPotential(std::array<std::size_t, 3> identifiers, BendType type,
                             std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBendParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BendType::Fixed:
      parameters[0] *= Units::DegreesToRadians;
      break;
    case BendType::Rigid:
      break;
    case BendType::Harmonic:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case BendType::CoreShell:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case BendType::Quartic:
      // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
      // ======================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case BendType::CFF_Quartic:
      // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
      // =====================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case BendType::HarmonicCosine:
      // p_0*(1+cos(p_1*theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::DegreesToRadians;
      break;
    case BendType::Cosine:
      // p_0*(1+cos(p_1*theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::DegreesToRadians;
      break;
    case BendType::Tafipolsky:
      // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
      // ===============================================
      // p_0/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case BendType::MM3:
    case BendType::MM3_inplane:
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

std::string BendPotential::print() const
{
  switch (type)
  {
    case BendType::Fixed:
      return std::format("{} - {} - {} : FIXED\n", identifiers[0], identifiers[1], identifiers[2]);
    case BendType::Rigid:
      return std::format("{} - {} - {} : RIGID\n", identifiers[0], identifiers[1], identifiers[2]);
    case BendType::Harmonic:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} : HARMONIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case BendType::CoreShell:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} : CORE_SHELL p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case BendType::Quartic:
      // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
      // ======================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      return std::format(
          "{} - {} - {} : QUARTIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2/k_B={:g} [K/rad^3], p_3/k_B={:g} "
          "[K/rad^4]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1],
          parameters[2] * Units::EnergyToKelvin, parameters[3] * Units::EnergyToKelvin);
    case BendType::CFF_Quartic:
      // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
      // =====================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      return std::format(
          "{} - {} - {} : CFF_QUARTIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2/k_B={:g} [K/rad^3], p_3/k_B={:g} "
          "[K/rad^4]\n",
          identifiers[0], identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin);
    case BendType::HarmonicCosine:
      // (1/2)*p_0*(cos(theta)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} : HARMONIC_COSINE p_0/k_B={:g} [K], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case BendType::Cosine:
      // p_0*(1+cos(p_1*theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      return std::format("{} - {} - {} : COSINE p_0/k_B={:g} [K], p_1={:g} [-], p_2={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2] * Units::RadiansToDegrees);
    case BendType::Tafipolsky:
      // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
      // ===============================================
      // p_0/k_B [K]
      return std::format("{} - {} - {} : TAFIPOLSKY p_0/k_B={:g} [K]\n", identifiers[0], identifiers[1], identifiers[2],
                         parameters[0] * Units::EnergyToKelvin);
    case BendType::MM3:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} : MM3 p_0/k_B={:g} [mdyne A/rad^2], p_1={:g} [degrees]]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKCalPerMol,
                         parameters[1] * Units::RadiansToDegrees);
    case BendType::MM3_inplane:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} : MM3 p_0/k_B={:g} [mdyne A/rad^2], p_1={:g} [degrees]]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKCalPerMol,
                         parameters[1] * Units::RadiansToDegrees);
    default:
      std::unreachable();
  }
}

double BendPotential::generateBendAngle(RandomNumber &random, double beta) const
{
  double theta, sin_theta, energy;
  double temp, temp2;

  switch (type)
  {
    case BendType::Fixed:
      return parameters[0];
    case BendType::Rigid:
      return parameters[0];
    case BendType::Harmonic:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = 0.5 * parameters[0] * (theta - parameters[1]) * (theta - parameters[1]);
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::CoreShell:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = 0.5 * parameters[0] * (theta - parameters[1]) * (theta - parameters[1]);
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::Quartic:
      // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
      // ======================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = 0.5 * parameters[0] * std::pow(theta - parameters[1], 2) +
                 (1.0 / 3.0) * parameters[2] * std::pow(theta - parameters[1], 3) +
                 0.25 * parameters[3] * std::pow(theta - parameters[1], 4);
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::CFF_Quartic:
      // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
      // =====================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = parameters[0] * std::pow(theta - parameters[1], 2) +
                 parameters[2] * std::pow(theta - parameters[1], 3) +
                 parameters[3] * std::pow(theta - parameters[1], 4);
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::HarmonicCosine:
      // (1/2)*p_0*(cos(theta)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = 0.5 * parameters[0] * std::pow(std::cos(theta) - std::cos(parameters[1]), 2);
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::Cosine:
      // p_0*(1+cos(p_1*theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = parameters[0] * (1.0 + std::cos(parameters[1] * theta - parameters[2]));
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::Tafipolsky:
      // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
      // ===============================================
      // p_0/k_B [K]
      do
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        energy = 0.5 * parameters[0] * (1.0 + std::cos(theta)) * (1.0 + std::cos(2.0 * theta));
      } while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    case BendType::MM3:
    case BendType::MM3_inplane:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      {
        theta = std::numbers::pi * random.uniform();
        sin_theta = std::sin(theta);
        temp = (theta - parameters[1]) * Units::RadiansToDegrees;
        temp2 = temp * temp;
        energy = parameters[0] * temp2 *
                 (1.0 - 0.014 * temp + 5.6e-5 * temp2 - 7.0e-7 * temp * temp2 + 2.2e-8 * temp2 * temp2);
      }
      while (random.uniform() > (sin_theta * sin_theta) * std::exp(-beta * energy));
      return theta;
    default:
      std::unreachable();
  }
}

double BendPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                      const std::optional<const double3> &posD) const
{
  double cos_theta, theta;
  double temp, temp2;

  // For four atoms, compute the in-plane angles
  if (posD.has_value())
  {
    double3 dr_ad = posA - posD.value();
    double r_ad = std::sqrt(double3::dot(dr_ad, dr_ad));
    dr_ad /= r_ad;

    double3 dr_bd = posB - posD.value();
    double r_bd = std::sqrt(double3::dot(dr_bd, dr_bd));
    dr_bd /= r_bd;

    double3 dr_cd = posC - posD.value();
    double r_cd = std::sqrt(double3::dot(dr_cd, dr_cd));
    dr_cd /= r_cd;

    double3 t = double3::cross(r_ad, r_cd);
    double rt2 = double3::dot(t, t);
    double delta = -double3::dot(t, r_bd) / rt2;

    double3 ip = posB + delta * t;

    double3 ap = posA - ip;
    double rap2 = double3::dot(ap, ap);
    ap /= rap2;

    double3 cp = posC - ip;
    double rcp2 = double3::dot(cp, cp);
    cp /= rcp2;

    cos_theta = double3::dot(ap, cp);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    theta = std::acos(cos_theta);
  }
  else
  {
    double3 dr_ab = posA - posB;
    double r_ab = std::sqrt(double3::dot(dr_ab, dr_ab));
    dr_ab /= r_ab;

    double3 dr_cb = posC - posB;
    double r_cb = std::sqrt(double3::dot(dr_cb, dr_cb));
    dr_cb /= r_cb;

    cos_theta = double3::dot(dr_ab, dr_cb);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    theta = std::acos(cos_theta);
  }

  switch (type)
  {
    case BendType::Fixed:
      return 0.0;
    case BendType::Rigid:
      return 0.0;
    case BendType::Harmonic:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      return 0.5 * parameters[0] * (theta - parameters[1]) * (theta - parameters[1]);
    case BendType::CoreShell:
      // (1/2)p_0*(theta-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      return 0.5 * parameters[0] * (theta - parameters[1]) * (theta - parameters[1]);
    case BendType::Quartic:
      // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
      // ======================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      temp = theta - parameters[1];
      temp2 = temp * temp;
      return 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
             0.25 * parameters[3] * temp2 * temp2;
    case BendType::CFF_Quartic:
      // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
      // =====================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2/k_B [K/rad^3]
      // p_3/k_B [K/rad^4]
      temp = theta - parameters[1];
      temp2 = temp * temp;
      return parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
    case BendType::HarmonicCosine:
      // (1/2)*p_0*(cos(theta)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      temp = cos_theta - parameters[1];
      temp2 = temp * temp;
      return 0.5 * parameters[0] * temp2;
    case BendType::Cosine:
      // p_0*(1+cos(p_1*theta-p_2))
      // ===============================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      temp = parameters[1] * theta - parameters[2];
      return parameters[0] * (1.0 + std::cos(temp));
    case BendType::Tafipolsky:
      // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
      // ===============================================
      // p_0/k_B [K]
      return 0.5 * parameters[0] * (1.0 + std::cos(theta)) * (1.0 + std::cos(2.0 * theta));
    case BendType::MM3:
    case BendType::MM3_inplane:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      temp = (theta - parameters[1]) * Units::RadiansToDegrees;
      temp2 = temp * temp;
      return parameters[0] * temp2 *
             (1.0 - 0.014 * temp + 5.6e-5 * temp2 - 7.0e-7 * temp * temp2 + 2.2e-8 * temp2 * temp2);
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
