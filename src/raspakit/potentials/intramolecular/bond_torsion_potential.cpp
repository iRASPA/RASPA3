module;

module bond_torsion_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import units;

BondTorsionPotential::BondTorsionPotential(std::array<std::size_t, 4> identifiers, BondTorsionType type,
                                           std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfBondTorsionParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondTorsionType::MM3:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      parameters[0] *= Units::KCalPerMolToEnergy;
      parameters[1] *= Units::KCalPerMolToEnergy;
      parameters[2] *= Units::KCalPerMolToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string BondTorsionPotential::print() const
{
  switch (type)
  {
    case BondTorsionType::MM3:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      return std::format(
          "{} - {} - {} - {} : MM3 p_0={:g} [kcal/A mole], p_1={:g} [kcal/A mole], p_2={:g} [kcal/A mole], "
          "p_3/k_B={:g} [A]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKCalPerMol,
          parameters[1] * Units::EnergyToKCalPerMol, parameters[2] * Units::EnergyToKCalPerMol, parameters[3]);
      break;
    default:
      std::unreachable();
  }
}

double BondTorsionPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                             const double3 &posD) const
{
  double temp;

  double3 Dab = posA - posB;

  double3 Dcb = posC - posB;
  double r_bc = std::sqrt(double3::dot(Dcb, Dcb));
  Dcb /= r_bc;

  double3 Ddc = posD - posC;

  double dot_ab = double3::dot(Dab, Dcb);
  double dot_dc = double3::dot(Ddc, Dcb);

  double3 dr = (Dab - dot_ab * Dcb).normalized();
  double3 ds = (Ddc - dot_dc * Dcb).normalized();

  // compute Cos(Phi)
  // Phi is defined in protein convention Phi(trans)=Pi
  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  double cos_phi2 = cos_phi * cos_phi;

  switch (type)
  {
    case BondTorsionType::MM3:
      // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
      // =====================================================================================
      // p_0     [kcal/A mole]
      // p_1     [kcal/A mole]
      // p_2     [kcal/A mole]
      // p_3     [A]
      temp = (r_bc - parameters[3]);
      return parameters[0] * temp * cos_phi + parameters[1] * temp * (2.0 * cos_phi2 - 1.0) +
             parameters[2] * temp * (4.0 * cos_phi2 * cos_phi - 3.0 * cos_phi);
    default:
      std::unreachable();
  }
}

std::tuple<double, std::array<double3, 4>, double3x3> BondTorsionPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD) const
{
  double3x3 strain_derivative{};
  double U{}, DCos{}, gamsa{}, gamsb{}, gamsc{};

  double3 Dab = posA - posB;
  double rab = std::sqrt(double3::dot(Dab, Dab));

  double3 Dcb = posC - posB;
  double rbc = std::sqrt(double3::dot(Dcb, Dcb));
  double3 Dcb_unit = Dcb / rbc;

  double3 Ddc = posD - posC;
  double rcd = std::sqrt(double3::dot(Ddc, Ddc));

  double dot_ab = double3::dot(Dab, Dcb_unit);
  double dot_cd = double3::dot(Ddc, Dcb_unit);

  double3 dr = Dab - dot_ab * Dcb_unit;
  double r = std::sqrt(double3::dot(dr, dr));
  dr /= r;

  double3 ds = Ddc - dot_cd * Dcb_unit;
  double s = std::sqrt(double3::dot(ds, ds));
  ds /= s;

  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  double cos_phi2 = cos_phi * cos_phi;

  switch (type)
  {
    case BondTorsionType::MM3:
    {
      const double temp = (rbc - parameters[3]);
      U = parameters[0] * temp * cos_phi + parameters[1] * temp * (2.0 * cos_phi2 - 1.0) +
          parameters[2] * temp * (4.0 * cos_phi2 * cos_phi - 3.0 * cos_phi);
      DCos = parameters[0] * temp + 4.0 * parameters[1] * temp * cos_phi +
             parameters[2] * temp * (12.0 * cos_phi2 - 3.0);
      gamsa = 0.0;
      gamsb = -(parameters[0] * cos_phi + parameters[1] * (2.0 * cos_phi2 - 1.0) +
                parameters[2] * (4.0 * cos_phi2 * cos_phi - 3.0 * cos_phi));
      gamsc = 0.0;
      break;
    }
    default:
      std::unreachable();
  }

  const double d = dot_ab / rbc;
  const double e = dot_cd / rbc;

  double3 dtA = (ds - cos_phi * dr) / r;
  double3 dtD = (dr - cos_phi * ds) / s;
  double3 dtB = dtA * (d - 1.0) + e * dtD;
  double3 dtC = -dtD * (e + 1.0) - d * dtA;

  double3 Dab_unit = Dab / rab;
  double3 Ddc_unit = Ddc / rcd;

  double3 du_da = DCos * dtA - gamsa * Dab_unit;
  double3 du_db = DCos * dtB + gamsb * Dcb_unit + gamsa * Dab_unit;
  double3 du_dc = DCos * dtC - gamsb * Dcb_unit + gamsc * Ddc_unit;
  double3 du_dd = DCos * dtD - gamsc * Ddc_unit;

  strain_derivative.ax += rbc * Dcb.x * (du_dc.x + du_dd.x) + rab * Dab.x * du_da.x + rcd * Ddc.x * du_dd.x;
  strain_derivative.bx += rbc * Dcb.y * (du_dc.x + du_dd.x) + rab * Dab.y * du_da.x + rcd * Ddc.y * du_dd.x;
  strain_derivative.cx += rbc * Dcb.z * (du_dc.x + du_dd.x) + rab * Dab.z * du_da.x + rcd * Ddc.z * du_dd.x;
  strain_derivative.ay += rbc * Dcb.x * (du_dc.y + du_dd.y) + rab * Dab.x * du_da.y + rcd * Ddc.x * du_dd.y;
  strain_derivative.by += rbc * Dcb.y * (du_dc.y + du_dd.y) + rab * Dab.y * du_da.y + rcd * Ddc.y * du_dd.y;
  strain_derivative.cy += rbc * Dcb.z * (du_dc.y + du_dd.y) + rab * Dab.z * du_da.y + rcd * Ddc.z * du_dd.y;
  strain_derivative.az += rbc * Dcb.x * (du_dc.z + du_dd.z) + rab * Dab.x * du_da.z + rcd * Ddc.x * du_dd.z;
  strain_derivative.bz += rbc * Dcb.y * (du_dc.z + du_dd.z) + rab * Dab.y * du_da.z + rcd * Ddc.y * du_dd.z;
  strain_derivative.cz += rbc * Dcb.z * (du_dc.z + du_dd.z) + rab * Dab.z * du_da.z + rcd * Ddc.z * du_dd.z;

  return {U, {du_da, du_db, du_dc, du_dd}, strain_derivative};
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondTorsionPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondTorsionPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondTorsionPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BondTorsionPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
