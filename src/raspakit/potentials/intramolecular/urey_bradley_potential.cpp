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

module urey_bradley_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

UreyBradleyPotential::UreyBradleyPotential(std::array<std::size_t, 2> identifiers, UreyBradleyType type,
                                           std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfUreyBradleyParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case UreyBradleyType::Fixed:
      break;
    case UreyBradleyType::Harmonic:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::CoreShellSpring:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::Morse:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::LJ_12_6:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::LennardJones:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::Buckingham:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::RestrainedHarmonic:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::Quartic:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::CFF_Quartic:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case UreyBradleyType::MM3:
      parameters[0] *= 71.94 * Units::KCalPerMolToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string UreyBradleyPotential::print() const
{
  switch (type)
  {
    case UreyBradleyType::Fixed:
      return std::format("{} - {} : FIXED\n", identifiers[0], identifiers[1]);
    case UreyBradleyType::Harmonic:
      return std::format("{} - {} : HARMONIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case UreyBradleyType::CoreShellSpring:
      return std::format("{} - {} : CORE_SHELL_SPRING p_0/k_B={:g} [K/Å^2]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin);
    case UreyBradleyType::Morse:
      return std::format("{} - {} : MORSE p_0/k_B={:g} [K/Å^2], p_1={:g} [Å^-1], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1], parameters[2]);
    case UreyBradleyType::LJ_12_6:
      return std::format("{} - {} : LJ_12_6 p_0/k_B={:g} [K Å^12], p_1={:g} [K Å^6]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case UreyBradleyType::LennardJones:
      return std::format("{} - {} : LENNARD_JONES p_0/k_B={:g} [K], p_1={:g} [Å]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case UreyBradleyType::Buckingham:
      return std::format("{} - {} : BUCKINGHAM p_0/k_B={:g} [K], p_1={:g} [Å^-1], p_2/k_B={:g} [K Å^6]\n",
                         identifiers[0], identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2] * Units::EnergyToKelvin);
    case UreyBradleyType::RestrainedHarmonic:
      return std::format("{} - {} : BUCKINGHAM p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1], parameters[2]);
    case UreyBradleyType::Quartic:
      return std::format("{} - {} : QUARTIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [K/Å^3], p_3={:g} [K/Å^4]\n",
                         identifiers[0], identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2] * Units::EnergyToKelvin, parameters[3] * Units::EnergyToKelvin);
    case UreyBradleyType::CFF_Quartic:
      return std::format(
          "{} - {} : CFF_QUARTIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [K/Å^3], p_3={:g} [K/Å^4]\n",
          identifiers[0], identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
          parameters[2] * Units::EnergyToKelvin, parameters[3] * Units::EnergyToKelvin);
    case UreyBradleyType::MM3:
      return std::format("{} - {} : MM3 p_0/k_B={:g} [mdyne/Å molecule], p_1={:g} [Å]\n", identifiers[0],
                         identifiers[1], parameters[0] / (71.94 * Units::KCalPerMolToEnergy), parameters[1]);
    default:
      std::unreachable();
  }
}

double UreyBradleyPotential::generateUreyBradleyLength(RandomNumber &random, double beta) const
{
  double ureyBradley_length, ran1, ran2;
  double temp, temp2, exp_term;
  double r1, energy;

  switch (type)
  {
    case UreyBradleyType::Fixed:
      return parameters[0];
    case UreyBradleyType::Harmonic:
      // 0.5 * p0 * SQR(r - p1);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      // p_1     [Å]       reference ureyBradley distance
      ran1 = 1.0 / std::sqrt(beta * parameters[0]);
      ran2 = 1.0 / (parameters[1] + 3.0 * ran1) * (parameters[1] + 3.0 * ran1);
      do ureyBradley_length = parameters[1] + ran1 * random.Gaussian();
      while ((random.uniform() > ureyBradley_length * ureyBradley_length * ran2) || (ureyBradley_length <= 0.0));
      return ureyBradley_length;
    case UreyBradleyType::CoreShellSpring:
      // 0.5 * p0 * SQR(r);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      ran1 = 1.0 / std::sqrt(beta * parameters[0]);
      ran2 = 1.0 / (parameters[1] + 3.0 * ran1) * (parameters[1] + 3.0 * ran1);
      do ureyBradley_length = ran1 * random.Gaussian();
      while ((random.uniform() > ureyBradley_length * ureyBradley_length * ran2) || (ureyBradley_length <= 0.0));
      return ureyBradley_length;
    case UreyBradleyType::Morse:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [Å^-1]    parameter
      // p_2     [Å]       reference ureyBradley distance
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        energy =
            parameters[0] * (std::pow(1.0 - std::exp(-parameters[1] * (ureyBradley_length - parameters[2])), 2) - 1.0);
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::LJ_12_6:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K Å^12]
      // p_1/k_B [K Å^6]
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        temp = std::pow(1.0 / (ureyBradley_length * ureyBradley_length), 3);
        energy = parameters[0] * temp * temp - parameters[1] * temp;
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        temp = std::pow(parameters[1] / (ureyBradley_length * ureyBradley_length), 3);
        energy = 4.0 * parameters[0] * (temp * (temp - 1.0));
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      do
      {
        ureyBradley_length = 3.0 * random.uniform() + 0.8;
        temp = parameters[2] * std::pow(1.0 / (ureyBradley_length * ureyBradley_length), 3);
        exp_term = parameters[0] * std::exp(-parameters[1] * ureyBradley_length);
        energy = exp_term - temp;
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::RestrainedHarmonic:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2     [Å]
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        r1 = ureyBradley_length - parameters[1];
        energy = 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
                 parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::Quartic:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        temp = ureyBradley_length - parameters[1];
        temp2 = temp * temp;
        energy = 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
                 0.25 * parameters[3] * temp2 * temp2;
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::CFF_Quartic:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        temp = ureyBradley_length - parameters[1];
        temp2 = temp * temp;
        energy = parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    case UreyBradleyType::MM3:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/Å molecule]
      // p_1     [Å]
      do
      {
        ureyBradley_length = 3.0 * random.uniform();
        temp = ureyBradley_length - parameters[1];
        temp2 = std::pow(ureyBradley_length - parameters[1], 2);
        energy = parameters[0] * temp2 * (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
      } while (random.uniform() > (ureyBradley_length * ureyBradley_length) * std::exp(-beta * energy));
      return ureyBradley_length;
    default:
      std::unreachable();
  }
}

double UreyBradleyPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  double temp, temp2;
  double r1, rri;

  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  switch (type)
  {
    case UreyBradleyType::Fixed:
      return 0.0;
    case UreyBradleyType::Harmonic:
      // 0.5 * p0 * SQR(r - p1);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      // p_1     [Å]       reference ureyBradley distance
      return 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]);
    case UreyBradleyType::CoreShellSpring:
      // 0.5 * p0 * SQR(r);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      return 0.5 * parameters[0] * r * r;
    case UreyBradleyType::Morse:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [Å^-1]    parameter
      // p_2     [Å]       reference ureyBradley distance
      temp = std::exp(parameters[1] * (parameters[2] - r));
      return parameters[0] * ((1.0 - temp) * (1.0 - temp) - 1.0);
    case UreyBradleyType::LJ_12_6:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K Å^12]
      // p_1/k_B [K Å^6]
      rri = (1.0 / rr);
      temp = rri * rri * rri;
      return parameters[0] * temp * temp - parameters[1] * temp;
    case UreyBradleyType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      rri = (parameters[1] / rr);
      temp = rri * rri * rri;
      return 4.0 * parameters[0] * (temp * (temp - 1.0));
    case UreyBradleyType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      rri = (parameters[1] / rr);
      temp = rri * rri * rri;
      return parameters[0] * std::exp(-parameters[1] * r) - parameters[2] * temp;
    case UreyBradleyType::RestrainedHarmonic:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2     [Å]
      r1 = r - parameters[1];
      return 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
             parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
    case UreyBradleyType::Quartic:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      temp = r - parameters[1];
      temp2 = temp * temp;
      return 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
             0.25 * parameters[3] * temp2 * temp2;
    case UreyBradleyType::CFF_Quartic:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      temp = r - parameters[1];
      temp2 = temp * temp;
      return parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
    case UreyBradleyType::MM3:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/Å molecule]
      // p_1     [Å]
      temp = r - parameters[1];
      temp2 = temp * temp;
      return parameters[0] * temp2 * (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const UreyBradleyPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, UreyBradleyPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'UreyBradleyPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("UreyBradleyPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
