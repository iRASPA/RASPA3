module;

module bend_bend_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import units;

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

std::tuple<double, std::array<double3, 4>, double3x3> BendBendPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD) const
{
  double3x3 strain_derivative{};
  double3 du_da{}, du_db{}, du_dc{}, du_dd{};

  double3 Dab_vec = posA - posB;
  double rab = std::sqrt(double3::dot(Dab_vec, Dab_vec));
  double3 Dab = Dab_vec / rab;

  double3 Dbc_vec = posC - posB;
  double rbc = std::sqrt(double3::dot(Dbc_vec, Dbc_vec));
  double3 Dbc = Dbc_vec / rbc;

  double3 Dbd_vec = posD - posB;
  double rbd = std::sqrt(double3::dot(Dbd_vec, Dbd_vec));
  double3 Dbd = Dbd_vec / rbd;

  double cos_theta1 = std::clamp(double3::dot(Dab, Dbc), -1.0, 1.0);
  double theta1 = std::acos(cos_theta1);
  double sin_theta1 = std::max(1.0e-8, std::sqrt(1.0 - cos_theta1 * cos_theta1));

  double cos_theta2 = std::clamp(double3::dot(Dab, Dbd), -1.0, 1.0);
  double theta2 = std::acos(cos_theta2);
  double sin_theta2 = std::max(1.0e-8, std::sqrt(1.0 - cos_theta2 * cos_theta2));

  double U{};
  double d_theta1{}, d_theta2{};

  switch (type)
  {
    case BendBendType::CVFF:
    case BendBendType::CFF:
      U = parameters[0] * (theta1 - parameters[1]) * (theta2 - parameters[2]);
      d_theta1 = parameters[0] * (theta2 - parameters[2]) / sin_theta1;
      d_theta2 = parameters[0] * (theta1 - parameters[1]) / sin_theta2;
      break;
    case BendBendType::MM3:
      U = -parameters[0] * (Units::RadiansToDegrees * Units::RadiansToDegrees) * (theta1 - parameters[1]) *
          (theta2 - parameters[2]);
      d_theta1 = parameters[0] * (Units::RadiansToDegrees * Units::RadiansToDegrees) * (theta2 - parameters[2]) /
                 sin_theta1;
      d_theta2 = parameters[0] * (Units::RadiansToDegrees * Units::RadiansToDegrees) * (theta1 - parameters[1]) /
                 sin_theta2;
      break;
    default:
      std::unreachable();
  }

  auto add_bend_contribution = [&](double d_theta, double cos_theta, const double3 &dir_ab, const double3 &dir_other,
                                 double r_ab, double r_other, const double3 &vec_ab, const double3 &vec_other,
                                 double3 &grad_a, double3 &grad_b, double3 &grad_other)
  {
    double3 fa = {
        -d_theta * (dir_other.x - cos_theta * dir_ab.x) / r_ab,
        -d_theta * (dir_other.y - cos_theta * dir_ab.y) / r_ab,
        -d_theta * (dir_other.z - cos_theta * dir_ab.z) / r_ab};
    double3 fo = {
        -d_theta * (dir_ab.x - cos_theta * dir_other.x) / r_other,
        -d_theta * (dir_ab.y - cos_theta * dir_other.y) / r_other,
        -d_theta * (dir_ab.z - cos_theta * dir_other.z) / r_other};
    double3 fb = -(fa + fo);
    grad_a += fa;
    grad_b += fb;
    grad_other += fo;

    strain_derivative.ax += r_ab * vec_ab.x * fa.x + r_other * vec_other.x * fo.x;
    strain_derivative.bx += r_ab * vec_ab.y * fa.x + r_other * vec_other.y * fo.x;
    strain_derivative.cx += r_ab * vec_ab.z * fa.x + r_other * vec_other.z * fo.x;
    strain_derivative.ay += r_ab * vec_ab.x * fa.y + r_other * vec_other.x * fo.y;
    strain_derivative.by += r_ab * vec_ab.y * fa.y + r_other * vec_other.y * fo.y;
    strain_derivative.cy += r_ab * vec_ab.z * fa.y + r_other * vec_other.z * fo.y;
    strain_derivative.az += r_ab * vec_ab.x * fa.z + r_other * vec_other.x * fo.z;
    strain_derivative.bz += r_ab * vec_ab.y * fa.z + r_other * vec_other.y * fo.z;
    strain_derivative.cz += r_ab * vec_ab.z * fa.z + r_other * vec_other.z * fo.z;
  };

  add_bend_contribution(d_theta1, cos_theta1, Dab, Dbc, rab, rbc, Dab_vec, Dbc_vec, du_da, du_db, du_dc);
  add_bend_contribution(d_theta2, cos_theta2, Dab, Dbd, rab, rbd, Dab_vec, Dbd_vec, du_da, du_db, du_dd);

  return {U, {du_da, du_db, du_dc, du_dd}, strain_derivative};
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
