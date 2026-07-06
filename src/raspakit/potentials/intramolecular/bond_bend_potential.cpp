module;

module bond_bend_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import units;

BondBendPotential::BondBendPotential(std::array<std::size_t, 4> identifiers, BondBendType type,
                                     std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfBondBendParameters); ++i)
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
      parameters[0] *= 2.51118 * Units::KCalPerMolToEnergy;
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
          "{} - {} - {} : CFF p_0={:g} [degrees], p_1/k_B={:g} [K/A/rad], p_2={:g} [A], p_3/k_B={:g} [K/A/rad], "
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
                         parameters[0] / (2.51118 * Units::KCalPerMolToEnergy), parameters[1], parameters[2],
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
    default:
      std::unreachable();
  }
}

double BondBendPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                          [[maybe_unused]] const double3 &posD) const
{
  double temp, temp2;

  double3 dr_ab = posA - posB;
  double r_ab = std::sqrt(double3::dot(dr_ab, dr_ab));
  dr_ab /= r_ab;

  double3 dr_cb = posC - posB;
  double r_cb = std::sqrt(double3::dot(dr_cb, dr_cb));
  dr_cb /= r_cb;

  double cos_theta = std::clamp(double3::dot(dr_ab, dr_cb), -1.0, 1.0);
  double theta = std::acos(cos_theta);

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
      return (theta - parameters[0]) * (parameters[1] * (r_ab - parameters[2]) + parameters[3] * (r_cb - parameters[4]));
    case BondBendType::MM3:
      // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
      // =====================================
      // p_0     [mdyne/rad]
      // p_1     [A]
      // p_2     [A]
      // p_3     [degrees]
      return parameters[0] * ((r_ab - parameters[1]) + (r_cb - parameters[2])) * Units::RadiansToDegrees *
             (theta - parameters[3]);
    case BondBendType::TruncatedHarmonic:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
      // ================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      temp = theta - parameters[1];
      return 0.5 * parameters[0] * temp * temp *
             std::exp(-(std::pow(r_ab, 8) + std::pow(r_cb, 8)) / std::pow(parameters[2], 8));
    case BondBendType::ScreenedHarmonic:
      // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      temp = theta - parameters[1];
      return 0.5 * parameters[0] * temp * temp * std::exp(-(r_ab / parameters[2] + r_cb / parameters[3]));
    case BondBendType::ScreenedVessal:
      // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
      // ============================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      // p_2     [A]
      // p_3     [A]
      temp = (parameters[1] - std::numbers::pi) * (parameters[1] - std::numbers::pi) -
             (theta - std::numbers::pi) * (theta - std::numbers::pi);
      return (parameters[0] / (8.0 * (theta - std::numbers::pi) * (theta - std::numbers::pi))) * temp * temp *
             std::exp(-(r_ab / parameters[2] + r_cb / parameters[3]));
    case BondBendType::TruncatedVessal:
      // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
      //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
      // ============================================================================
      // p_0/k_B [K/rad^(4+p_2)]
      // p_1     [degrees]
      // p_2     [-]
      // p_3     [A]
      temp = theta - parameters[1];
      temp2 = theta + parameters[1] - 2.0 * std::numbers::pi;
      return parameters[0] *
                 (std::pow(theta, parameters[2]) * temp * temp * temp2 * temp2 -
                  0.5 * parameters[2] * std::pow(std::numbers::pi, parameters[2] - 1.0) * temp * temp *
                      std::pow(std::numbers::pi - parameters[1], 3.0)) *
             std::exp(-(std::pow(r_ab, 8) + std::pow(r_cb, 8)) / std::pow(parameters[3], 8));
    default:
      std::unreachable();
  }
}

std::tuple<double, std::array<double3, 3>, double3x3> BondBendPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB, const double3 &posC) const
{
  double temp, temp2, exp_term, pow_rab8, pow_rbc8;
  double gamma, gamsa, gamsc, pterm;
  double3x3 strain_derivative{};

  double3 Rab_vec = posA - posB;
  double rab = std::sqrt(double3::dot(Rab_vec, Rab_vec));
  double3 Rab = Rab_vec / rab;

  double3 Rbc_vec = posC - posB;
  double rbc = std::sqrt(double3::dot(Rbc_vec, Rbc_vec));
  double3 Rbc = Rbc_vec / rbc;

  double cos_theta = std::clamp(double3::dot(Rab, Rbc), -1.0, 1.0);
  double theta = std::acos(cos_theta);
  double sin_theta = std::max(1.0e-8, std::sqrt(1.0 - cos_theta * cos_theta));

  switch (type)
  {
    case BondBendType::CVFF:
    case BondBendType::CFF:
      pterm = (theta - parameters[0]) * (parameters[1] * (rab - parameters[2]) + parameters[3] * (rbc - parameters[4]));
      gamma = (parameters[1] * (rab - parameters[2]) + parameters[3] * (rbc - parameters[4])) / sin_theta;
      gamsa = -parameters[1] * (theta - parameters[0]);
      gamsc = -parameters[3] * (theta - parameters[0]);
      break;
    case BondBendType::MM3:
      pterm = parameters[0] * ((rab - parameters[1]) + (rbc - parameters[2])) * Units::RadiansToDegrees *
              (theta - parameters[3]);
      gamma = parameters[0] * Units::RadiansToDegrees * ((rab - parameters[1]) + (rbc - parameters[2])) / sin_theta;
      gamsa = -parameters[0] * Units::RadiansToDegrees * (theta - parameters[3]);
      gamsc = gamsa;
      break;
    case BondBendType::TruncatedHarmonic:
      temp = theta - parameters[1];
      exp_term = std::exp(-(std::pow(rab, 8) + std::pow(rbc, 8)) / std::pow(parameters[2], 8));
      pterm = 0.5 * parameters[0] * temp * temp * exp_term;
      gamma = (parameters[0] * (theta - parameters[1]) * exp_term) / sin_theta;
      gamsa = (8.0 * pterm / std::pow(parameters[2], 8)) * std::pow(rab, 7);
      gamsc = (8.0 * pterm / std::pow(parameters[2], 8)) * std::pow(rbc, 7);
      break;
    case BondBendType::ScreenedHarmonic:
      temp = theta - parameters[1];
      exp_term = std::exp(-(rab / parameters[2] + rbc / parameters[3]));
      pterm = 0.5 * parameters[0] * temp * temp * exp_term;
      gamma = (parameters[0] * (theta - parameters[1]) * exp_term) / sin_theta;
      gamsa = pterm / parameters[2];
      gamsc = pterm / parameters[3];
      break;
    case BondBendType::ScreenedVessal:
      temp = (parameters[1] - std::numbers::pi) * (parameters[1] - std::numbers::pi) -
             (theta - std::numbers::pi) * (theta - std::numbers::pi);
      exp_term = std::exp(-(rab / parameters[2] + rbc / parameters[3]));
      pterm = (parameters[0] / (8.0 * (theta - std::numbers::pi) * (theta - std::numbers::pi))) * temp * temp *
              exp_term;
      gamma = (parameters[0] / (2.0 * (theta - std::numbers::pi) * (theta - std::numbers::pi)) *
               ((parameters[1] - std::numbers::pi) * (parameters[1] - std::numbers::pi) -
                (theta - std::numbers::pi) * (theta - std::numbers::pi)) *
               (theta - std::numbers::pi) * exp_term) /
              sin_theta;
      gamsa = pterm / parameters[2];
      gamsc = pterm / parameters[3];
      break;
    case BondBendType::TruncatedVessal:
      temp = theta - parameters[1];
      temp2 = theta + parameters[1] - 2.0 * std::numbers::pi;
      pow_rab8 = std::pow(rab, 8);
      pow_rbc8 = std::pow(rbc, 8);
      exp_term = std::exp(-(pow_rab8 + pow_rbc8) / std::pow(parameters[3], 8));
      pterm = parameters[0] *
              (std::pow(theta, parameters[2]) * temp * temp * temp2 * temp2 -
               0.5 * parameters[2] * std::pow(std::numbers::pi, parameters[2] - 1.0) * temp * temp *
                   std::pow(std::numbers::pi - parameters[1], 3.0)) *
              exp_term;
      gamma = (parameters[0] *
               (std::pow(theta, parameters[2] - 1.0) * (theta - parameters[1]) * (theta + parameters[1] - 2.0 * std::numbers::pi) *
                    ((parameters[2] + 4.0) * (theta - parameters[1]) * (theta - parameters[1]) -
                     2.0 * std::numbers::pi * (parameters[2] + 2.0) * theta + parameters[2] * parameters[1] *
                                                                                      (2.0 * std::numbers::pi - parameters[1])) -
                parameters[2] * std::pow(std::numbers::pi, parameters[2] - 1.0) * (theta - parameters[1]) *
                    std::pow(std::numbers::pi - parameters[1], 3.0)) *
               exp_term) /
              sin_theta;
      gamsa = (8.0 * pterm / std::pow(parameters[3], 8)) * std::pow(rab, 7);
      gamsc = (8.0 * pterm / std::pow(parameters[3], 8)) * std::pow(rbc, 7);
      break;
    default:
      std::unreachable();
  }

  double3 du_da = {
      -(gamma * (Rbc.x - Rab.x * cos_theta) / rab + gamsa * Rab.x),
      -(gamma * (Rbc.y - Rab.y * cos_theta) / rab + gamsa * Rab.y),
      -(gamma * (Rbc.z - Rab.z * cos_theta) / rab + gamsa * Rab.z)};
  double3 du_dc = {
      -(gamma * (Rab.x - Rbc.x * cos_theta) / rbc + gamsc * Rbc.x),
      -(gamma * (Rab.y - Rbc.y * cos_theta) / rbc + gamsc * Rbc.y),
      -(gamma * (Rab.z - Rbc.z * cos_theta) / rbc + gamsc * Rbc.z)};
  double3 du_db = -(du_da + du_dc);

  strain_derivative.ax = rab * Rab_vec.x * du_da.x + rbc * Rbc_vec.x * du_dc.x;
  strain_derivative.bx = rab * Rab_vec.y * du_da.x + rbc * Rbc_vec.y * du_dc.x;
  strain_derivative.cx = rab * Rab_vec.z * du_da.x + rbc * Rbc_vec.z * du_dc.x;
  strain_derivative.ay = rab * Rab_vec.x * du_da.y + rbc * Rbc_vec.x * du_dc.y;
  strain_derivative.by = rab * Rab_vec.y * du_da.y + rbc * Rbc_vec.y * du_dc.y;
  strain_derivative.cy = rab * Rab_vec.z * du_da.y + rbc * Rbc_vec.z * du_dc.y;
  strain_derivative.az = rab * Rab_vec.x * du_da.z + rbc * Rbc_vec.x * du_dc.z;
  strain_derivative.bz = rab * Rab_vec.y * du_da.z + rbc * Rbc_vec.y * du_dc.z;
  strain_derivative.cz = rab * Rab_vec.z * du_da.z + rbc * Rbc_vec.z * du_dc.z;

  return {pterm, {du_da, du_db, du_dc}, strain_derivative};
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
