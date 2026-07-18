module;

module inversion_bend_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import units;

InversionBendPotential::InversionBendPotential(std::array<std::size_t, 4> identifiers, InversionBendType type,
                                               std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfInversionBendParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::Harmonic2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case InversionBendType::HarmonicCosine:
    case InversionBendType::HarmonicCosine2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] = std::cos(parameters[1] * Units::DegreesToRadians);
      break;
    case InversionBendType::Planar:
    case InversionBendType::Planar2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      break;
    case InversionBendType::MM3:
      // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      parameters[0] *= 0.02191418 * Units::KCalPerMolToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    default:
      std::unreachable();
  }
}

std::string InversionBendPotential::print() const
{
  switch (type)
  {
    case InversionBendType::Harmonic:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case InversionBendType::HarmonicCosine:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC_COSINE p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, std::acos(parameters[1]) * Units::RadiansToDegrees);
    case InversionBendType::Planar:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      return std::format("{} - {} - {} - {} : PLANAR p_0/k_B={:g} [K]\n", identifiers[0], identifiers[1],
                         identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin);
    case InversionBendType::Harmonic2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC2 p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case InversionBendType::HarmonicCosine2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC_COSINE2 p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, std::acos(parameters[1]) * Units::RadiansToDegrees);
    case InversionBendType::Planar2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      return std::format("{} - {} - {} - {} : PLANAR2 p_0/k_B={:g} [K]\n", identifiers[0], identifiers[1],
                         identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin);
    case InversionBendType::MM3:
      // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : MM3 p_0/k_B={:g} [mdyne A/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    default:
      std::unreachable();
  }
}

double InversionBendPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                               const double3 &posD) const
{
  double c, chi;
  double temp, temp2;

  double3 Rab = posA - posB;
  double rab2 = double3::dot(Rab, Rab);
  double rrab = std::sqrt(rab2);

  double3 Rbc = posC - posB;
  double3 Rbd = posD - posB;
  double3 Rcd = posD - posC;
  double3 Rad = posD - posA;

  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::HarmonicCosine:
    case InversionBendType::Planar:
      // w is a vector perpendicular to the B-C-D plane
      // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
      temp = double3::dot(Rbc, Rbd);
      c = double3::dot(Rbc, Rbc) * double3::dot(Rbd, Rbd) - temp * temp;
      break;
    case InversionBendType::Harmonic2:
    case InversionBendType::HarmonicCosine2:
    case InversionBendType::Planar2:
    case InversionBendType::MM3:
      // w is a vector perpendicular to the A-C-D plane
      // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
      temp = double3::dot(Rad, Rcd);
      c = double3::dot(Rcd, Rcd) * double3::dot(Rad, Rad) - temp * temp;
      break;
    default:
      std::unreachable();
  }

  double e = double3::dot(Rab, double3::cross(Rbd, Rbc));
  double cos_chi = std::sqrt(double3::dot(Rab, Rab) - e * e / c) / rrab;
  cos_chi = std::clamp(cos_chi, -1.0, 1.0);

  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::Harmonic2:
      // (1/2)*p_0*(chi-p_1)^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      chi = std::acos(cos_chi);
      temp = chi - parameters[1];
      return 0.5 * parameters[0] * temp * temp;
    case InversionBendType::HarmonicCosine:
    case InversionBendType::HarmonicCosine2:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      temp = cos_chi - parameters[1];
      return 0.5 * parameters[0] * temp * temp;
      break;
    case InversionBendType::Planar:
    case InversionBendType::Planar2:
      // (1/2)*p_0*(1-cos(phi))
      // ===============================================
      // p_0/k_B [K]
      return parameters[0] * (1.0 - cos_chi);
    case InversionBendType::MM3:
      // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
      // =================================================================================================
      // p_0/k_B [mdyne A/rad^2]
      // p_1     [degrees]
      chi = std::acos(cos_chi);
      temp = (chi - parameters[1]) * Units::RadiansToDegrees;
      temp2 = temp * temp;
      return parameters[0] * temp2 *
             (1.0 - 0.014 * temp + 5.6e-5 * temp2 - 7.0e-7 * temp * temp2 + 2.2e-8 * temp2 * temp2);
    default:
      std::unreachable();
  }
}

std::tuple<double, std::array<double3, 4>, double3x3> InversionBendPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD) const
{
  double c, e, dot, chi, temp, temp2, dedcos, term;
  double3 dccd_a, dccd_c, dccd_d, deed_a, deed_c, deed_d;
  double3 du_da, du_db, du_dc, du_dd;
  double3x3 strain_derivative{};

  double3 Rab = posA - posB;
  double rab2 = double3::dot(Rab, Rab);
  double rrab = std::sqrt(rab2);

  double3 Rbc = posC - posB;
  double rbc2 = double3::dot(Rbc, Rbc);

  double3 Rbd = posD - posB;
  double rbd2 = double3::dot(Rbd, Rbd);

  double3 Rac = posC - posA;
  double rac2 = double3::dot(Rac, Rac);

  double3 Rad = posD - posA;
  double rad2 = double3::dot(Rad, Rad);

  const bool use_acd_plane = type == InversionBendType::Harmonic2 || type == InversionBendType::HarmonicCosine2 ||
                             type == InversionBendType::Planar2 || type == InversionBendType::MM3;

  if (use_acd_plane)
  {
    dot = double3::dot(Rad, Rac);
    c = rac2 * rad2 - dot * dot;
  }
  else
  {
    dot = double3::dot(Rbc, Rbd);
    c = rbc2 * rbd2 - dot * dot;
  }

  e = double3::dot(Rab, double3::cross(Rbd, Rbc));
  double cos_chi = std::sqrt(std::max(0.0, rab2 - e * e / c)) / rrab;
  cos_chi = std::clamp(cos_chi, -1.0, 1.0);

  double U{};
  switch (type)
  {
    case InversionBendType::Harmonic:
    case InversionBendType::Harmonic2:
      chi = std::acos(cos_chi);
      temp = chi - parameters[1];
      U = 0.5 * parameters[0] * temp * temp;
      dedcos = -std::copysign(1.0, e) * (parameters[0] * temp / std::sqrt(c * (rab2 - e * e / c)));
      break;
    case InversionBendType::HarmonicCosine:
    case InversionBendType::HarmonicCosine2:
      chi = std::acos(cos_chi);
      temp = cos_chi - parameters[1];
      U = 0.5 * parameters[0] * temp * temp;
      dedcos = std::copysign(1.0, e) * parameters[0] * temp * std::sin(chi) / std::sqrt(c * (rab2 - e * e / c));
      break;
    case InversionBendType::Planar:
    case InversionBendType::Planar2:
      chi = std::acos(cos_chi);
      U = parameters[0] * (1.0 - cos_chi);
      dedcos = -std::copysign(1.0, e) * parameters[0] * std::sin(chi) / std::sqrt(c * (rab2 - e * e / c));
      break;
    case InversionBendType::MM3:
      chi = std::acos(cos_chi);
      temp = (chi - parameters[1]) * Units::RadiansToDegrees;
      temp2 = temp * temp;
      U = parameters[0] * temp2 *
          (1.0 - 0.014 * temp + 5.6e-5 * temp2 - 7.0e-7 * temp * temp2 + 2.2e-8 * temp2 * temp2);
      dedcos = -std::copysign(1.0, e) * parameters[0] * temp * Units::RadiansToDegrees *
               (2.0 - 3.0 * 0.014 * temp + 4.0 * 5.6e-5 * temp2 - 5.0 * 7.0e-7 * temp * temp2 +
                6.0 * 2.2e-8 * temp2 * temp2) /
               std::sqrt(c * (rab2 - e * e / c));
      break;
    default:
      std::unreachable();
  }

  if (use_acd_plane)
  {
    term = e / c;
    dccd_c = {(Rac.x * rad2 - Rad.x * dot) * term, (Rac.y * rad2 - Rad.y * dot) * term,
              (Rac.z * rad2 - Rad.z * dot) * term};
    dccd_d = {(Rad.x * rac2 - Rac.x * dot) * term, (Rad.y * rac2 - Rac.y * dot) * term,
              (Rad.z * rac2 - Rac.z * dot) * term};
    dccd_a = -(dccd_c + dccd_d);
  }
  else
  {
    term = e / c;
    dccd_a = {0.0, 0.0, 0.0};
    dccd_c = {(Rbc.x * rbd2 - Rbd.x * dot) * term, (Rbc.y * rbd2 - Rbd.y * dot) * term,
              (Rbc.z * rbd2 - Rbd.z * dot) * term};
    dccd_d = {(Rbd.x * rbc2 - Rbc.x * dot) * term, (Rbd.y * rbc2 - Rbc.y * dot) * term,
              (Rbd.z * rbc2 - Rbc.z * dot) * term};
  }

  // Combined e- and rab-dependent parts of dcos(chi)/datom (validated against finite differences
  // of the energy; the cross-product operand order fixes a sign error in the original RASPA2
  // transcription that produced gradients inconsistent with calculateEnergy).
  term = e / rab2;
  deed_a = {Rbc.y * Rbd.z - Rbc.z * Rbd.y + Rab.x * term, Rbc.z * Rbd.x - Rbc.x * Rbd.z + Rab.y * term,
            Rbc.x * Rbd.y - Rbc.y * Rbd.x + Rab.z * term};
  deed_c = {Rbd.y * Rab.z - Rbd.z * Rab.y, Rbd.z * Rab.x - Rbd.x * Rab.z, Rbd.x * Rab.y - Rbd.y * Rab.x};
  deed_d = {Rab.y * Rbc.z - Rab.z * Rbc.y, Rab.z * Rbc.x - Rab.x * Rbc.z, Rab.x * Rbc.y - Rab.y * Rbc.x};

  du_da = dedcos * (dccd_a + deed_a);
  du_dc = dedcos * (dccd_c + deed_c);
  du_dd = dedcos * (dccd_d + deed_d);
  du_db = -(du_da + du_dc + du_dd);

  strain_derivative.ax = Rab.x * du_da.x + Rbc.x * du_dc.x + Rbd.x * du_dd.x;
  strain_derivative.bx = Rab.y * du_da.x + Rbc.y * du_dc.x + Rbd.y * du_dd.x;
  strain_derivative.cx = Rab.z * du_da.x + Rbc.z * du_dc.x + Rbd.z * du_dd.x;
  strain_derivative.ay = Rab.x * du_da.y + Rbc.x * du_dc.y + Rbd.x * du_dd.y;
  strain_derivative.by = Rab.y * du_da.y + Rbc.y * du_dc.y + Rbd.y * du_dd.y;
  strain_derivative.cy = Rab.z * du_da.y + Rbc.z * du_dc.y + Rbd.z * du_dd.y;
  strain_derivative.az = Rab.x * du_da.z + Rbc.x * du_dc.z + Rbd.x * du_dd.z;
  strain_derivative.bz = Rab.y * du_da.z + Rbc.y * du_dc.z + Rbd.y * du_dd.z;
  strain_derivative.cz = Rab.z * du_da.z + Rbc.z * du_dc.z + Rbd.z * du_dd.z;

  return {U, {du_da, du_db, du_dc, du_dd}, strain_derivative};
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const InversionBendPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, InversionBendPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'InversionBendPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("InversionBendPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
