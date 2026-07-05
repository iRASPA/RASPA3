module;

module bend_torsion_potential;

import std;

import archive;
import randomnumbers;
import double3;

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
