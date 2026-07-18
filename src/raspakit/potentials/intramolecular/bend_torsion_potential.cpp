module;

module bend_torsion_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import units;

BendTorsionPotential::BendTorsionPotential(std::array<std::size_t, 4> identifiers, BendTorsionType type,
                                           std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfBendTorsionParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BendTorsionType::Smoothed:
      // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
      // ============================================
      // p_0/k_B [K/rad^2]
      // p_1     [-]
      // p_2     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[2] *= Units::DegreesToRadians;
      break;
    case BendTorsionType::SmoothedThreeCosine:
    case BendTorsionType::Nicholas:
    case BendTorsionType::SmoothedCFF:
    case BendTorsionType::SmoothedCFF2:
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      break;
    case BendTorsionType::CVFF:
    case BendTorsionType::CFF:
    case BendTorsionType::SmoothedCFF3:
      // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)  (optionally smoothed)
      // =============================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::DegreesToRadians;
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
      return std::format(
          "{} - {} - {} - {} : SMOOTHED p_0/k_B={:g} [K/rad^2], p_1={:g} [-], p_2={:g} [degrees]\n", identifiers[0],
          identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin, parameters[1],
          parameters[2] * Units::RadiansToDegrees);
    case BendTorsionType::SmoothedThreeCosine:
      return std::format(
          "{} - {} - {} - {} : SMOOTHED_THREE_COSINE p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin);
    case BendTorsionType::Nicholas:
      return std::format("{} - {} - {} - {} : NICHOLAS p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::EnergyToKelvin,
                         parameters[2] * Units::EnergyToKelvin);
    case BendTorsionType::CFF:
      return std::format(
          "{} - {} - {} - {} : CFF p_0/k_B={:g} [K/rad^3], p_1={:g} [degrees], p_2={:g} [degrees]\n", identifiers[0],
          identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2] * Units::RadiansToDegrees);
    case BendTorsionType::SmoothedCFF:
      return std::format("{} - {} - {} - {} : SMOOTHED_CFF p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::EnergyToKelvin,
                         parameters[2] * Units::EnergyToKelvin);
    case BendTorsionType::SmoothedCFF2:
      return std::format("{} - {} - {} - {} : SMOOTHED_CFF2 p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::EnergyToKelvin,
                         parameters[2] * Units::EnergyToKelvin);
    case BendTorsionType::SmoothedCFF3:
      return std::format(
          "{} - {} - {} - {} : SMOOTHED_CFF3 p_0/k_B={:g} [K/rad^3], p_1={:g} [degrees], p_2={:g} [degrees]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2] * Units::RadiansToDegrees);
    case BendTorsionType::CVFF:
      return std::format(
          "{} - {} - {} - {} : CVFF p_0/k_B={:g} [K/rad^3], p_1={:g} [degrees], p_2={:g} [degrees]\n", identifiers[0],
          identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::RadiansToDegrees, parameters[2] * Units::RadiansToDegrees);
    default:
      std::unreachable();
  }
}

// Smoothing function S(theta): 1 below 170 degrees, tapering to 0 at 180 degrees
static double smoothing(double theta)
{
  const double on = 170.0 * (std::numbers::pi / 180.0);
  const double off = 180.0 * (std::numbers::pi / 180.0);

  if (theta < on) return 1.0;
  return (off - theta) * (off - theta) * (off + 2.0 * theta - 3.0 * on) / ((off - on) * (off - on) * (off - on));
}

double BendTorsionPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                             const double3 &posD) const
{
  double phi, sign;

  double3 Dab = posA - posB;
  double r_ab = std::sqrt(double3::dot(Dab, Dab));

  double3 Dbc = posC - posB;
  double r_bc = std::sqrt(double3::dot(Dbc, Dbc));
  Dbc /= r_bc;

  double3 Dcd = posD - posC;
  double r_cd = std::sqrt(double3::dot(Dcd, Dcd));

  double dot_ab = double3::dot(Dab, Dbc);
  double cos_theta1 = std::clamp(dot_ab / r_ab, -1.0, 1.0);
  double theta1 = std::acos(cos_theta1);

  double dot_cd = double3::dot(Dcd, Dbc);
  double cos_theta2 = std::clamp(-dot_cd / r_cd, -1.0, 1.0);
  double theta2 = std::acos(cos_theta2);

  double3 dr = (Dab - dot_ab * Dbc).normalized();
  double3 ds = (Dcd - dot_cd * Dbc).normalized();

  // compute Cos(Phi)
  // Phi is defined in protein convention Phi(trans)=Pi
  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  double cos_phi2 = cos_phi * cos_phi;

  switch (type)
  {
    case BendTorsionType::CVFF:
    case BendTorsionType::CFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
      // =====================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      return parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]) * cos_phi;
    case BendTorsionType::Smoothed:
      // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K/rad^2]
      // p_1     [-]
      // p_2     [degrees]
      sign = double3::dot(Dbc, double3::cross(double3::cross(Dab, Dbc), double3::cross(Dbc, Dcd)));
      phi = std::copysign(std::acos(cos_phi), sign);
      return parameters[0] * (1.0 + std::cos(parameters[1] * phi - parameters[2])) * smoothing(theta1) *
             smoothing(theta2);
    case BendTorsionType::SmoothedThreeCosine:
      // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return (0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
              0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2)) *
             smoothing(theta1) * smoothing(theta2);
    case BendTorsionType::Nicholas:
      // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return (0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
              0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2)) *
             smoothing(theta1);
    case BendTorsionType::SmoothedCFF:
      // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return (parameters[0] * (1.0 - cos_phi) + 2.0 * parameters[1] * (1.0 - cos_phi2) +
              parameters[2] * (1.0 + 3.0 * cos_phi - 4.0 * cos_phi * cos_phi2)) *
             smoothing(theta1) * smoothing(theta2);
    case BendTorsionType::SmoothedCFF2:
      // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return (parameters[0] * (1.0 + cos_phi) + parameters[2] +
              cos_phi * (-3.0 * parameters[2] + 2.0 * cos_phi * (parameters[1] + 2.0 * parameters[2] * cos_phi))) *
             smoothing(theta1) * smoothing(theta2);
    case BendTorsionType::SmoothedCFF3:
      // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
      // ======================================================================================
      // p_0/k_B [K/rad^3]
      // p_1     [degrees]
      // p_2     [degrees]
      return parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]) * cos_phi * smoothing(theta1) *
             smoothing(theta2);
    default:
      std::unreachable();
  }
}

static double smoothingDerivative(double theta)
{
  const double on = 170.0 * (std::numbers::pi / 180.0);
  const double off = 180.0 * (std::numbers::pi / 180.0);

  if (theta < on) return 0.0;
  return 6.0 * (off - theta) * (on - theta) / std::pow(off - on, 3);
}

std::tuple<double, std::array<double3, 4>, double3x3> BendTorsionPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD) const
{
  double phi, sign;
  double3x3 strain_derivative{};
  double3 du_da{}, du_db{}, du_dc{}, du_dd{};

  double3 Dab_vec = posA - posB;
  double rab = std::sqrt(double3::dot(Dab_vec, Dab_vec));

  double3 Dbc_vec = posC - posB;
  double rbc = std::sqrt(double3::dot(Dbc_vec, Dbc_vec));
  double3 Dbc = Dbc_vec / rbc;

  double3 Dcd_vec = posD - posC;
  double rcd = std::sqrt(double3::dot(Dcd_vec, Dcd_vec));

  double dot_ab = double3::dot(Dab_vec, Dbc);
  double cos_theta1 = std::clamp(dot_ab / rab, -1.0, 1.0);
  double theta1 = std::acos(cos_theta1);
  double sin_theta1 = std::max(1.0e-8, std::sqrt(1.0 - cos_theta1 * cos_theta1));

  double dot_cd = double3::dot(Dcd_vec, Dbc);
  double cos_theta2 = std::clamp(-dot_cd / rcd, -1.0, 1.0);
  double theta2 = std::acos(cos_theta2);
  double sin_theta2 = std::max(1.0e-8, std::sqrt(1.0 - cos_theta2 * cos_theta2));

  double3 dr = (Dab_vec - dot_ab * Dbc).normalized();
  double3 ds = (Dcd_vec - dot_cd * Dbc).normalized();

  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  double cos_phi2 = cos_phi * cos_phi;

  double U{}, DCos{}, DTheta1{}, DTheta2{};
  const double s1 = smoothing(theta1);
  const double s2 = smoothing(theta2);
  const double ds1 = smoothingDerivative(theta1);
  const double ds2 = smoothingDerivative(theta2);

  switch (type)
  {
    case BendTorsionType::CVFF:
    case BendTorsionType::CFF:
      U = parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]) * cos_phi;
      DCos = parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]);
      DTheta1 = parameters[0] * (theta2 - parameters[2]) * cos_phi / sin_theta1;
      DTheta2 = parameters[0] * (theta1 - parameters[1]) * cos_phi / sin_theta2;
      break;
    case BendTorsionType::Smoothed:
      sign = double3::dot(Dbc, double3::cross(double3::cross(Dab_vec, Dbc), double3::cross(Dbc, Dcd_vec)));
      phi = std::copysign(std::acos(cos_phi), sign);
      U = parameters[0] * (1.0 + std::cos(parameters[1] * phi - parameters[2])) * s1 * s2;
      DCos = (parameters[0] * parameters[1] * std::sin(parameters[1] * phi - parameters[2]) * s1 * s2) /
             std::copysign(std::max(1.0e-8, std::fabs(std::sin(phi))), std::sin(phi));
      DTheta1 = parameters[0] * (1.0 + std::cos(parameters[1] * phi - parameters[2])) * ds1 * s2 / sin_theta1;
      DTheta2 = parameters[0] * (1.0 + std::cos(parameters[1] * phi - parameters[2])) * s1 * ds2 / sin_theta2;
      break;
    case BendTorsionType::SmoothedThreeCosine:
    {
      const double torsion_part = 0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
                                  0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
      U = torsion_part * s1 * s2;
      DCos = (0.5 * parameters[0] - 2.0 * parameters[1] * cos_phi + 1.5 * parameters[2] * (4.0 * cos_phi2 - 1.0)) * s1 *
             s2;
      DTheta1 = torsion_part * ds1 * s2 / sin_theta1;
      DTheta2 = torsion_part * s2 * ds2 * s1 / sin_theta2;
      break;
    }
    case BendTorsionType::Nicholas:
    {
      const double torsion_part = 0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
                                  0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
      U = torsion_part * s1;
      DCos = (0.5 * parameters[0] - 2.0 * parameters[1] * cos_phi + 1.5 * parameters[2] * (4.0 * cos_phi2 - 1.0)) * s1;
      DTheta1 = torsion_part * ds1 / sin_theta1;
      DTheta2 = 0.0;
      break;
    }
    case BendTorsionType::SmoothedCFF:
    {
      const double torsion_part = parameters[0] * (1.0 - cos_phi) + 2.0 * parameters[1] * (1.0 - cos_phi2) +
                                  parameters[2] * (1.0 + 3.0 * cos_phi - 4.0 * cos_phi * cos_phi2);
      U = torsion_part * s1 * s2;
      DCos = (-parameters[0] - 4.0 * parameters[1] * cos_phi + 3.0 * parameters[2] * (1.0 - 4.0 * cos_phi2)) * s1 * s2;
      DTheta1 = torsion_part * ds1 * s2 / sin_theta1;
      DTheta2 = torsion_part * s1 * ds2 / sin_theta2;
      break;
    }
    case BendTorsionType::SmoothedCFF2:
    {
      const double torsion_part =
          parameters[0] * (1.0 + cos_phi) + parameters[2] +
          cos_phi * (-3.0 * parameters[2] + 2.0 * cos_phi * (parameters[1] + 2.0 * parameters[2] * cos_phi));
      U = torsion_part * s1 * s2;
      DCos = (parameters[0] - 3.0 * parameters[2] + 4.0 * cos_phi * (parameters[1] + 3.0 * parameters[2] * cos_phi)) * s1 *
             s2;
      DTheta1 = torsion_part * ds1 * s2 / sin_theta1;
      DTheta2 = torsion_part * s1 * ds2 / sin_theta2;
      break;
    }
    case BendTorsionType::SmoothedCFF3:
      U = parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]) * cos_phi * s1 * s2;
      DCos = parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]) * s1 * s2;
      DTheta1 = cos_phi * parameters[0] * (theta2 - parameters[2]) * s2 *
                (s1 + (theta1 - parameters[1]) * ds1) / sin_theta1;
      DTheta2 = cos_phi * parameters[0] * (theta1 - parameters[1]) * s1 *
                (s2 + (theta2 - parameters[2]) * ds2) / sin_theta2;
      break;
    default:
      std::unreachable();
  }

  const double r = std::sqrt(double3::dot(Dab_vec - dot_ab * Dbc, Dab_vec - dot_ab * Dbc));
  const double s = std::sqrt(double3::dot(Dcd_vec - dot_cd * Dbc, Dcd_vec - dot_cd * Dbc));
  const double d = dot_ab / rbc;
  const double e = dot_cd / rbc;

  double3 dtA = (ds - cos_phi * dr) / r;
  double3 dtD = (dr - cos_phi * ds) / s;
  double3 dtB = dtA * (d - 1.0) + e * dtD;
  double3 dtC = -dtD * (e + 1.0) - d * dtA;

  du_da += DCos * dtA;
  du_db += DCos * dtB;
  du_dc += DCos * dtC;
  du_dd += DCos * dtD;

  double3 Dab = Dab_vec / rab;
  double3 Dcd = Dcd_vec / rcd;

  // Bend theta1 (A-B-C, vertex B), unit arms Dab (towards A) and Dbc (towards C);
  // DTheta already carries the (dU/dtheta)/sin(theta) factor, so dU/dcos(theta) = -DTheta.
  if (DTheta1 != 0.0)
  {
    const double3 fa = (-DTheta1 / rab) * (Dbc - cos_theta1 * Dab);
    const double3 fc = (-DTheta1 / rbc) * (Dab - cos_theta1 * Dbc);
    du_da += fa;
    du_db -= fa + fc;
    du_dc += fc;
  }

  // Bend theta2 (B-C-D, vertex C), unit arms -Dbc (towards B) and Dcd (towards D).
  if (DTheta2 != 0.0)
  {
    const double3 fb = (-DTheta2 / rbc) * (Dcd + cos_theta2 * Dbc);
    const double3 fd = (-DTheta2 / rcd) * (-1.0 * Dbc - cos_theta2 * Dcd);
    du_db += fb;
    du_dc -= fb + fd;
    du_dd += fd;
  }

  // Molecular virial: the total gradient is translation invariant, so the strain derivative is
  // sum_i (x_i - x_B) (x) g_i with arms A-B, C-B, and D-B = (C-B) + (D-C).
  const double3 armA = Dab_vec;
  const double3 armC = Dbc_vec;
  const double3 armD = Dbc_vec + Dcd_vec;
  strain_derivative.ax = armA.x * du_da.x + armC.x * du_dc.x + armD.x * du_dd.x;
  strain_derivative.bx = armA.y * du_da.x + armC.y * du_dc.x + armD.y * du_dd.x;
  strain_derivative.cx = armA.z * du_da.x + armC.z * du_dc.x + armD.z * du_dd.x;
  strain_derivative.ay = armA.x * du_da.y + armC.x * du_dc.y + armD.x * du_dd.y;
  strain_derivative.by = armA.y * du_da.y + armC.y * du_dc.y + armD.y * du_dd.y;
  strain_derivative.cy = armA.z * du_da.y + armC.z * du_dc.y + armD.z * du_dd.y;
  strain_derivative.az = armA.x * du_da.z + armC.x * du_dc.z + armD.x * du_dd.z;
  strain_derivative.bz = armA.y * du_da.z + armC.y * du_dc.z + armD.y * du_dd.z;
  strain_derivative.cz = armA.z * du_da.z + armC.z * du_dc.z + armD.z * du_dd.z;

  return {U, {du_da, du_db, du_dc, du_dd}, strain_derivative};
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
