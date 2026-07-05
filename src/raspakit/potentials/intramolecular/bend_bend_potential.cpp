module;

module bend_bend_potential;

import std;

import archive;
import randomnumbers;
import double3;

BendBendPotential::BendBendPotential(std::array<std::size_t, 4> identifiers, BendBendType type,
                                     std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfBendBendParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BendBendType::CVFF:
    case BendBendType::CFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::DegreesToRadians;
      break;
    case BendBendType::MM3:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      parameters[0] *= 0.02191418 * Units::KCalPerMolToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      parameters[2] *= Units::DegreesToRadians;
      break;
    default:
      std::unreachable();
  }
}

std::string BendBendPotential::print() const
{
  switch (type)
  {
    case BendBendType::CVFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : CVFF p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::RadiansToDegrees,
                         parameters[2] * Units::RadiansToDegrees);
    case BendBendType::CFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : CFF p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::RadiansToDegrees,
                         parameters[2] * Units::RadiansToDegrees);
    case BendBendType::MM3:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : MM3 p_0={:g} [mdyne A/rad^2], p_1={:g} [degrees], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKCalPerMol / 0.02191418,
                         parameters[1] * Units::RadiansToDegrees, parameters[2] * Units::RadiansToDegrees);
    default:
      std::unreachable();
  }
}

double BendBendPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                          const double3 &posD) const
{
  // Theta1 is the angle A-B-C, Theta2 is the angle A-B-D (atom B is the central atom)
  double3 dr_ab = posA - posB;
  double r_ab = std::sqrt(double3::dot(dr_ab, dr_ab));
  dr_ab /= r_ab;

  double3 dr_bc = posC - posB;
  double r_bc = std::sqrt(double3::dot(dr_bc, dr_bc));
  dr_bc /= r_bc;

  double3 dr_bd = posD - posB;
  double r_bd = std::sqrt(double3::dot(dr_bd, dr_bd));
  dr_bd /= r_bd;

  double cos_theta1 = std::clamp(double3::dot(dr_ab, dr_bc), -1.0, 1.0);
  double theta1 = std::acos(cos_theta1);

  double cos_theta2 = std::clamp(double3::dot(dr_ab, dr_bd), -1.0, 1.0);
  double theta2 = std::acos(cos_theta2);

  switch (type)
  {
    case BendBendType::CVFF:
    case BendBendType::CFF:
      // p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0/k_B [K/rad^2)]
      // p_1     [degrees]
      // p_2     [degrees]
      return parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]);
    case BendBendType::MM3:
      // -p_0*(Theta1-p_1)*(Theta2-p_2)
      // ===================================
      // p_0     [mdyne A/rad^2]
      // p_1     [degrees]
      // p_2     [degrees]
      return -parameters[0] * (Units::RadiansToDegrees * Units::RadiansToDegrees) * (theta1 - parameters[1]) *
             (theta2 - parameters[2]);
    default:
      std::unreachable();
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendBendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendBendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BendBendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BendBendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
