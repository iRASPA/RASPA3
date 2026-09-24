module;

module bond_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import distance_potential_gradient_strain;
import distance_potential_gradient_hessian_strain;

BondPotential::BondPotential(std::array<std::size_t, 2> identifiers, BondType type,
                             std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfBondParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondType::None:
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
    case BondType::None:
      return std::format("{} - {} : NONE\n", identifiers[0], identifiers[1]);
    case BondType::Fixed:
      return std::format("{} - {} : FIXED  p_0={:g} [Å]\n", identifiers[0], identifiers[1], parameters[0]);
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
      return std::format("{} - {} : RESTRAINED_HARMONIC p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
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


// Exact rejection sampling of the bond-length density p(r) ~ r^2 exp(-beta u(r)) on [lo, hi]:
// uniform proposals against a grid-estimated envelope (the 10% headroom covers the discretization
// error of the grid maximum for these smooth one-dimensional densities). The former samplers
// accepted with r^2 exp(-beta u) directly, which exceeds one for r > 1 or u < 0 and silently
// truncates -- not a valid rejection scheme for the anharmonic bond types.
template <typename EnergyFunction>
static double sampleBondLengthRejection(RandomNumber &random, double lo, double hi, double beta,
                                        EnergyFunction energy)
{
  constexpr std::size_t numberOfGridPoints = 1024;
  double envelope = 0.0;
  for (std::size_t i = 0; i != numberOfGridPoints; ++i)
  {
    double r = lo + (hi - lo) * (static_cast<double>(i) + 0.5) / static_cast<double>(numberOfGridPoints);
    envelope = std::max(envelope, r * r * std::exp(-beta * energy(r)));
  }
  envelope *= 1.1;

  double bond_length, density;
  do
  {
    bond_length = lo + (hi - lo) * random.uniform();
    density = bond_length * bond_length * std::exp(-beta * energy(bond_length));
  } while (random.uniform() > density / envelope);
  return bond_length;
}

double BondPotential::generateBondLength(RandomNumber &random, double beta) const
{
  double bond_length, sigma, envelope;

  switch (type)
  {
    case BondType::None:
      return 1.0;
    case BondType::Fixed:
      return parameters[0];
    case BondType::Harmonic:
      // 0.5 * p0 * SQR(r - p1);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      // p_1     [A]       reference bond distance
      //
      // Stiff bonds: Gaussian proposals (the exact Boltzmann factor of the harmonic bond) against
      // the r^2 Jacobian; the envelope (p_1 + 6 sigma)^2 is valid up to a ~1e-9 truncation
      // probability. Soft or zero-strength bonds (where Gaussian proposals rarely land at positive
      // r) use the generic grid-envelope sampler on a bounded interval instead.
      sigma = std::sqrt(1.0 / (beta * parameters[0]));
      if (!(sigma < 0.5 * std::max(parameters[1], 1.0)))
      {
        return sampleBondLengthRejection(
            random, 0.0, parameters[1] + 3.0, beta,
            [&](double r) { return 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]); });
      }
      envelope = (parameters[1] + 6.0 * sigma) * (parameters[1] + 6.0 * sigma);
      do bond_length = random.Gaussian(parameters[1], sigma);
      while ((bond_length <= 0.0) || (random.uniform() > bond_length * bond_length / envelope));
      return bond_length;
    case BondType::CoreShellSpring:
      // 0.5 * p0 * SQR(r);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      sigma = std::sqrt(1.0 / (beta * parameters[0]));
      return sampleBondLengthRejection(random, 0.0, 6.0 * sigma, beta,
                                       [&](double r) { return 0.5 * parameters[0] * r * r; });
    case BondType::Morse:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference bond distance
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         return parameters[0] *
                                                (std::pow(1.0 - std::exp(-parameters[1] * (r - parameters[2])), 2) -
                                                 1.0);
                                       });
    case BondType::LJ_12_6:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = std::pow(1.0 / (r * r), 3);
                                         return parameters[0] * temp * temp - parameters[1] * temp;
                                       });
    case BondType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A]
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = std::pow((parameters[1] * parameters[1]) / (r * r), 3);
                                         return 4.0 * parameters[0] * (temp * (temp - 1.0));
                                       });
    case BondType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      return sampleBondLengthRejection(random, 0.8, 3.8, beta,
                                       [&](double r)
                                       {
                                         return parameters[0] * std::exp(-parameters[1] * r) -
                                                parameters[2] * std::pow(1.0 / (r * r), 3);
                                       });
    case BondType::RestrainedHarmonic:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double r1 = r - parameters[1];
                                         return 0.5 * parameters[0] *
                                                    std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
                                                parameters[0] * parameters[2] *
                                                    std::max(std::fabs(r1) - parameters[2], 0.0);
                                       });
    case BondType::Quartic:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = r - parameters[1];
                                         double temp2 = temp * temp;
                                         return 0.5 * parameters[0] * temp2 +
                                                (1.0 / 3.0) * parameters[2] * temp * temp2 +
                                                0.25 * parameters[3] * temp2 * temp2;
                                       });
    case BondType::CFF_Quartic:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
                                       [&](double r)
                                       {
                                         double temp = r - parameters[1];
                                         double temp2 = temp * temp;
                                         return parameters[0] * temp2 + parameters[2] * temp * temp2 +
                                                parameters[3] * temp2 * temp2;
                                       });
    case BondType::MM3:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/A molecule]
      // p_1     [A]
      return sampleBondLengthRejection(random, 0.0, 3.0, beta,
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

double BondPotential::logBoltzmannVolumeNormalization(double beta) const
{
  // Delta-distributed bond lengths: the density is delta(r - r0) on the radial measure r^2 dr,
  // contributing the fixed Jacobian r0^2.
  switch (type)
  {
    case BondType::None:
      return 0.0;  // r0 = 1.0
    case BondType::Fixed:
      return 2.0 * std::log(parameters[0]);
    default:
      break;
  }

  using ParameterArray = std::remove_cvref_t<decltype(parameters)>;
  using CacheKey = std::tuple<BondType, ParameterArray, double>;
  thread_local std::map<CacheKey, double> cache{};
  const CacheKey key{type, parameters, beta};
  if (auto it = cache.find(key); it != cache.end()) return it->second;

  // The integration interval mirrors the support of the density sampled by generateBondLength.
  double lo = 0.0;
  double hi = 3.0;
  switch (type)
  {
    case BondType::Harmonic:
    {
      const double sigma = std::sqrt(1.0 / (beta * parameters[0]));
      if (!(sigma < 0.5 * std::max(parameters[1], 1.0)))
      {
        hi = parameters[1] + 3.0;
      }
      else
      {
        lo = std::max(0.0, parameters[1] - 10.0 * sigma);
        hi = parameters[1] + 10.0 * sigma;
      }
      break;
    }
    case BondType::CoreShellSpring:
      hi = 6.0 * std::sqrt(1.0 / (beta * parameters[0]));
      break;
    case BondType::Buckingham:
      lo = 0.8;
      hi = 3.8;
      break;
    default:
      break;
  }

  constexpr std::size_t numberOfGridPoints = 4096;
  const double h = (hi - lo) / static_cast<double>(numberOfGridPoints);
  double integral = 0.0;
  for (std::size_t i = 0; i != numberOfGridPoints; ++i)
  {
    const double r = lo + (static_cast<double>(i) + 0.5) * h;
    integral += r * r * std::exp(-beta * calculateEnergy(double3{0.0, 0.0, 0.0}, double3{r, 0.0, 0.0}));
  }
  const double result = std::log(integral * h);

  cache[key] = result;
  return result;
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
    case BondType::None:
    case BondType::Fixed:
      return std::abs(r - parameters[0]) < 1e-10 ? 0.0 : 0.0;
      //return std::abs(r - parameters[0]) < 1e-10 ? 0.0 : std::numeric_limits<double>::max();
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
      rri = (parameters[1] * parameters[1]) / rr;
      temp = rri * rri * rri;
      return 4.0 * parameters[0] * (temp * (temp - 1.0));
    case BondType::Buckingham:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å^-1]
      // p_2/k_B [K Å^6]
      rri = 1.0 / rr;
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

std::tuple<double, std::array<double3, 2>, double3x3> BondPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB) const
{
  return Potentials::Internal::distancePotentialEnergyGradientStrain(type, parameters, posA, posB, false);
}

std::tuple<double, std::array<double3, 2>, double3x3, double, double>
BondPotential::potentialEnergyGradientHessianStrain(const double3 &posA, const double3 &posB) const
{
  return Potentials::Internal::distancePotentialEnergyGradientHessianStrain(type, parameters, posA, posB, false);
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
