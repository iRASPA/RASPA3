module;

module urey_bradley_potential;

import std;

import archive;
import randomnumbers;
import rejection_envelope;
import double3;
import double3x3;
import bond_potential;
import distance_potential_gradient_strain;

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

void UreyBradleyPotential::scaleEnergy(double factor)
{
  switch (type)
  {
    case UreyBradleyType::Fixed:
      break;
    case UreyBradleyType::Harmonic:
    case UreyBradleyType::CoreShellSpring:
    case UreyBradleyType::Morse:
    case UreyBradleyType::LennardJones:
    case UreyBradleyType::RestrainedHarmonic:
    case UreyBradleyType::MM3:
      parameters[0] *= factor;
      break;
    case UreyBradleyType::LJ_12_6:
      parameters[0] *= factor;
      parameters[1] *= factor;
      break;
    case UreyBradleyType::Buckingham:
      parameters[0] *= factor;
      parameters[2] *= factor;
      break;
    case UreyBradleyType::Quartic:
    case UreyBradleyType::CFF_Quartic:
      parameters[0] *= factor;
      parameters[2] *= factor;
      parameters[3] *= factor;
      break;
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
      return std::format("{} - {} : RESTRAINED_HARMONIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
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

double UreyBradleyPotential::generateUreyBradleyLength(RandomNumber& random, double beta) const
{
  const Potentials::PotentialTag tag(1, static_cast<int>(type), parameters);
  // The same samplers as BondPotential::generateBondLength: the density r^2 exp(-beta u(r)) is sampled by
  // rejection against a normalized envelope.
  double ureyBradley_length, sigma, envelope;

  switch (type)
  {
    case UreyBradleyType::Fixed:
      return parameters[0];
    case UreyBradleyType::Harmonic:
      // 0.5 * p0 * SQR(r - p1);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      // p_1     [Å]       reference ureyBradley distance
      //
      // Stiff terms: Gaussian proposals against the r^2 Jacobian with the envelope (p_1 + 6 sigma)^2; soft or
      // zero-strength terms use the grid-envelope sampler on a bounded interval instead.
      sigma = std::sqrt(1.0 / (beta * parameters[0]));
      if (!(sigma < 0.5 * std::max(parameters[1], 1.0)))
      {
        return sampleBondLengthRejection(tag, random, 0.0, parameters[1] + 3.0, beta, [&](double r)
                                         { return 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]); });
      }
      envelope = (parameters[1] + 6.0 * sigma) * (parameters[1] + 6.0 * sigma);
      do ureyBradley_length = random.Gaussian(parameters[1], sigma);
      while ((ureyBradley_length <= 0.0) || (random.uniform() > ureyBradley_length * ureyBradley_length / envelope));
      return ureyBradley_length;
    case UreyBradleyType::CoreShellSpring:
      // 0.5 * p0 * SQR(r);
      // ===============================================
      // p_0/k_B [K/Å^2]   force constant
      sigma = std::sqrt(1.0 / (beta * parameters[0]));
      return sampleBondLengthRejection(tag, random, 0.0, 6.0 * sigma, beta,
                                       [&](double r) { return 0.5 * parameters[0] * r * r; });
    case UreyBradleyType::Morse:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [Å^-1]    parameter
      // p_2     [Å]       reference ureyBradley distance
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta, [&](double r)
          { return parameters[0] * (std::pow(1.0 - std::exp(-parameters[1] * (r - parameters[2])), 2) - 1.0); });
    case UreyBradleyType::LJ_12_6:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K Å^12]
      // p_1/k_B [K Å^6]
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = std::pow(1.0 / (r * r), 3);
                                         return parameters[0] * temp * temp - parameters[1] * temp;
                                       });
    case UreyBradleyType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = std::pow((parameters[1] * parameters[1]) / (r * r), 3);
                                         return 4.0 * parameters[0] * (temp * (temp - 1.0));
                                       });
    case UreyBradleyType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      return sampleBondLengthRejection(tag, random, 0.8, 3.8, beta, [&](double r)
          { return parameters[0] * std::exp(-parameters[1] * r) - parameters[2] * std::pow(1.0 / (r * r), 3); });
    case UreyBradleyType::RestrainedHarmonic:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2     [Å]
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta,
          [&](double r)
          {
            double r1 = r - parameters[1];
            return 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
                   parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
          });
    case UreyBradleyType::Quartic:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = r - parameters[1];
                                         double temp2 = temp * temp;
                                         return 0.5 * parameters[0] * temp2 +
                                                (1.0 / 3.0) * parameters[2] * temp * temp2 +
                                                0.25 * parameters[3] * temp2 * temp2;
                                       });
    case UreyBradleyType::CFF_Quartic:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/Å^2]
      // p_1     [Å]
      // p_2/k_B [K/Å^3]
      // p_3/k_B [K/Å^4]
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = r - parameters[1];
                                         double temp2 = temp * temp;
                                         return parameters[0] * temp2 + parameters[2] * temp * temp2 +
                                                parameters[3] * temp2 * temp2;
                                       });
    case UreyBradleyType::MM3:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/Å molecule]
      // p_1     [Å]
      return sampleBondLengthRejection(tag, random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = r - parameters[1];
                                         double temp2 = temp * temp;
                                         return parameters[0] * temp2 *
                                                (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
                                       });
    default:
      std::unreachable();
  }
}

double UreyBradleyPotential::calculateEnergy(const double3& posA, const double3& posB) const
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
      rri = (parameters[1] * parameters[1]) / rr;
      temp = rri * rri * rri;
      return 4.0 * parameters[0] * (temp * (temp - 1.0));
    case UreyBradleyType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      rri = 1.0 / rr;
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

std::tuple<double, std::array<double3, 2>, double3x3> UreyBradleyPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB) const
{
  // 'BondType' has an extra leading 'None' enumerator, so the Urey-Bradley enumerators are shifted by one
  // (UreyBradleyType::Harmonic == 1, BondType::Harmonic == 2).
  const auto bond_type = static_cast<BondType>(std::to_underlying(type) + 1);
  return Potentials::Internal::distancePotentialEnergyGradientStrain(bond_type, parameters, posA, posB, false);
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
