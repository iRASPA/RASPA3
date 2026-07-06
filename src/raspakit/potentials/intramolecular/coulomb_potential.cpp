module;

module coulomb_potential;

import std;

import archive;
import randomnumbers;
import units;
import double3;
import double3x3;

CoulombPotential::CoulombPotential(std::array<std::size_t, 2> identifiers, CoulombType type,
                                   double chargeA, double chargeB, double scaling)
    : identifiers(identifiers), type(type), chargeA(chargeA), chargeB(chargeB), scaling(scaling)
{
  switch (type)
  {
    case CoulombType::Coulomb:
      break;
    default:
      std::unreachable();
  }
}

std::string CoulombPotential::print() const
{
  switch (type)
  {
    case CoulombType::Coulomb:
      return std::format("{} - {} : COULOMB p_0={:g} [{}], p_1={:g} [{}] scaling: {} [-]\n", identifiers[0], identifiers[1],
                         chargeA, "e", chargeB, "e", scaling);
    default:
      std::unreachable();
  }
}

double CoulombPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  switch (type)
  {
    case CoulombType::Coulomb:
      return scaling * Units::CoulombicConversionFactor * chargeA * chargeB / r;
    default:
      std::unreachable();
  }
}

std::tuple<double, std::array<double3, 2>, double3x3> CoulombPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB) const
{
  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  double U{};
  double DF{};

  switch (type)
  {
    case CoulombType::Coulomb:
      U = scaling * Units::CoulombicConversionFactor * chargeA * chargeB / r;
      if (r > 0.0)
      {
        DF = scaling * Units::CoulombicConversionFactor * chargeA * chargeB / (rr * r);
      }
      break;
    default:
      std::unreachable();
  }

  double3 du_dr = DF * dr;
  double3 du_da = du_dr;
  double3 du_db = -du_dr;

  double3x3 strain_derivative{};
  strain_derivative.ax = dr.x * du_dr.x;
  strain_derivative.bx = dr.y * du_dr.x;
  strain_derivative.cx = dr.z * du_dr.x;
  strain_derivative.ay = dr.x * du_dr.y;
  strain_derivative.by = dr.y * du_dr.y;
  strain_derivative.cy = dr.z * du_dr.y;
  strain_derivative.az = dr.x * du_dr.z;
  strain_derivative.bz = dr.y * du_dr.z;
  strain_derivative.cz = dr.z * du_dr.z;

  return {U, {du_da, du_db}, strain_derivative};
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CoulombPotential &b)
{
  archive << b.versionNumber;

  archive << b.type;
  archive << b.identifiers;
  archive << b.chargeA;
  archive << b.chargeB;
  archive << b.scaling;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CoulombPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'CoulombPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.type;
  archive >> b.identifiers;
  archive >> b.chargeA;
  archive >> b.chargeB;
  archive >> b.scaling;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("CoulombPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
