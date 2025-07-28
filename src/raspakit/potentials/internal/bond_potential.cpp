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

module bond_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

BondPotential::BondPotential(std::array<std::size_t, 2> identifiers, BondType type,
                             std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBondParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondType::Fixed:
      break;
    case BondType::Harmonic:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case BondType::CoreShellSpring:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case BondType::Morse:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case BondType::LJ_12_6:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      break;
    case BondType::LennardJones:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case BondType::Buckingham:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      break;
    case BondType::RestrainedHarmonic:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case BondType::Quartic:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case BondType::CFF_Quartic:
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case BondType::MM3:
      parameters[0] *= 71.94 * Units::KCalPerMolToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string BondPotential::print() const
{
  switch (type)
  {
    case BondType::Fixed:
      return std::format("{} - {} : FIXED\n", identifiers[0], identifiers[1]);
    case BondType::Harmonic:
      return std::format("{} - {} : HARMONIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case BondType::CoreShellSpring:
      return std::format("{} - {} : CORE_SHELL_SPRING p_0/k_B={:g} [K/Å^2]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin);
    case BondType::Morse:
      return std::format("{} - {} : MORSE p_0/k_B={:g} [K/Å^2], p_1={:g} [Å^-1], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1], parameters[2]);
    case BondType::LJ_12_6:
      return std::format("{} - {} : LJ_12_6 p_0/k_B={:g} [K Å^12], p_1={:g} [K Å^6]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case BondType::LennardJones:
      return std::format("{} - {} : LENNARD_JONES p_0/k_B={:g} [K], p_1={:g} [Å]\n", identifiers[0], identifiers[1],
                         parameters[0] * Units::EnergyToKelvin, parameters[1]);
    case BondType::Buckingham:
      return std::format("{} - {} : BUCKINGHAM p_0/k_B={:g} [K], p_1={:g} [Å^-1], p_2/k_B={:g} [K Å^6]\n",
                         identifiers[0], identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2] * Units::EnergyToKelvin);
    case BondType::RestrainedHarmonic:
      return std::format("{} - {} : BUCKINGHAM p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1], parameters[2]);
    case BondType::Quartic:
      return std::format("{} - {} : QUARTIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [K/Å^3], p_3={:g} [K/Å^4]\n",
                         identifiers[0], identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2] * Units::EnergyToKelvin, parameters[3] * Units::EnergyToKelvin);
    case BondType::CFF_Quartic:
      return std::format(
          "{} - {} : CFF_QUARTIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [K/Å^3], p_3={:g} [K/Å^4]\n",
          identifiers[0], identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1],
          parameters[2] * Units::EnergyToKelvin, parameters[3] * Units::EnergyToKelvin);
    case BondType::MM3:
      return std::format("{} - {} : MM3 p_0/k_B={:g} [mdyne/Å molecule], p_1={:g} [Å]\n", identifiers[0],
                         identifiers[1], parameters[0] / (71.94 * Units::KCalPerMolToEnergy), parameters[1]);
    default:
      std::unreachable();
  }
}

double BondPotential::generateBondLength(RandomNumber &random, double beta) const
{
  double bond_length, ran1, ran2;
  double temp, temp2, exp_term;
  double r1, energy;

  switch (type)
  {
    case BondType::Fixed:
      return parameters[0];
    case BondType::Harmonic:
      // 0.5 * p0 * SQR(r - p1);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      // p_1     [Å]       reference bond distance
      ran1 = 1.0 / std::sqrt(beta * parameters[0]);
      ran2 = 1.0 / (parameters[1] + 3.0 * ran1) * (parameters[1] + 3.0 * ran1);
      do bond_length = parameters[1] + ran1 * random.Gaussian();
      while ((random.uniform() > bond_length * bond_length * ran2) || (bond_length <= 0.0));
      return bond_length;
    case BondType::CoreShellSpring:
      // 0.5 * p0 * SQR(r);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      ran1 = 1.0 / std::sqrt(beta * parameters[0]);
      ran2 = 1.0 / (parameters[1] + 3.0 * ran1) * (parameters[1] + 3.0 * ran1);
      do bond_length = ran1 * random.Gaussian();
      while ((random.uniform() > bond_length * bond_length * ran2) || (bond_length <= 0.0));
      return bond_length;
    case BondType::Morse:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [Å^-1]    parameter
      // p_2     [Å]       reference bond distance
      do
      {
        bond_length = 3.0 * random.uniform();
        energy = parameters[0] * (std::pow(1.0 - std::exp(-parameters[1] * (bond_length - parameters[2])), 2) - 1.0);
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::LJ_12_6:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K Å^12]
      // p_1/k_B [K Å^6]
      do
      {
        bond_length = 3.0 * random.uniform();
        temp = std::pow(1.0 / (bond_length * bond_length), 3);
        energy = parameters[0] * temp * temp - parameters[1] * temp;
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      do
      {
        bond_length = 3.0 * random.uniform();
        temp = std::pow(parameters[1] / (bond_length * bond_length), 3);
        energy = 4.0 * parameters[0] * (temp * (temp - 1.0));
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      do
      {
        bond_length = 3.0 * random.uniform() + 0.8;
        temp = parameters[2] * std::pow(1.0 / (bond_length * bond_length), 3);
        exp_term = parameters[0] * std::exp(-parameters[1] * bond_length);
        energy = exp_term - temp;
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::RestrainedHarmonic:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2     [Å]
      do
      {
        bond_length = 3.0 * random.uniform();
        r1 = bond_length - parameters[1];
        energy = 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
                 parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::Quartic:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      do
      {
        bond_length = 3.0 * random.uniform();
        temp = bond_length - parameters[1];
        temp2 = temp * temp;
        energy = 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
                 0.25 * parameters[3] * temp2 * temp2;
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::CFF_Quartic:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      do
      {
        bond_length = 3.0 * random.uniform();
        temp = bond_length - parameters[1];
        temp2 = temp * temp;
        energy = parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    case BondType::MM3:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/Å molecule]
      // p_1     [Å]
      do
      {
        bond_length = 3.0 * random.uniform();
        temp = bond_length - parameters[1];
        temp2 = std::pow(bond_length - parameters[1], 2);
        energy = parameters[0] * temp2 * (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
      } while (random.uniform() > (bond_length * bond_length) * std::exp(-beta * energy));
      return bond_length;
    default:
      std::unreachable();
  }
}

double BondPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  double temp, temp2;
  double r1, rri;

  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  switch (type)
  {
    case BondType::Fixed:
      return 0.0;
    case BondType::Harmonic:
      // 0.5 * p0 * SQR(r - p1);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      // p_1     [Å]       reference bond distance
      return 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]);
    case BondType::CoreShellSpring:
      // 0.5 * p0 * SQR(r);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      return 0.5 * parameters[0] * r * r;
    case BondType::Morse:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [Å^-1]    parameter
      // p_2     [Å]       reference bond distance
      temp = std::exp(parameters[1] * (parameters[2] - r));
      return parameters[0] * ((1.0 - temp) * (1.0 - temp) - 1.0);
    case BondType::LJ_12_6:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K Å^12]
      // p_1/k_B [K Å^6]
      rri = (1.0 / rr);
      temp = rri * rri * rri;
      return parameters[0] * temp * temp - parameters[1] * temp;
    case BondType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      rri = (parameters[1] / rr);
      temp = rri * rri * rri;
      return 4.0 * parameters[0] * (temp * (temp - 1.0));
    case BondType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      rri = (parameters[1] / rr);
      temp = rri * rri * rri;
      return parameters[0] * std::exp(-parameters[1] * r) - parameters[2] * temp;
    case BondType::RestrainedHarmonic:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2     [Å]
      r1 = r - parameters[1];
      return 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
             parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
    case BondType::Quartic:
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
    case BondType::CFF_Quartic:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      temp = r - parameters[1];
      temp2 = temp * temp;
      return parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
    case BondType::MM3:
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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BondPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
