module;

module bond_bond_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;

BondBondPotential::BondBondPotential(std::array<std::size_t, 3> identifiers, BondBondType type,
                                     std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(parameters.size(), maximumNumberOfBondBondParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case BondBondType::CVFF:
    case BondBondType::CFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      parameters[0] *= Units::KelvinToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string BondBondPotential::print() const
{
  switch (type)
  {
    case BondBondType::CVFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return std::format("{} - {} - {} : CVFF p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2]);
    case BondBondType::CFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return std::format("{} - {} - {} : CFF p_0/k_B={:g} [K/Å^2], p_1={:g} [Å], p_2={:g} [Å]\n", identifiers[0],
                         identifiers[1], identifiers[2], parameters[0] * Units::EnergyToKelvin, parameters[1],
                         parameters[2]);
    default:
      std::unreachable();
  }
}

double BondBondPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC) const
{
  double3 dr_ab = posA - posB;
  double rr_ab = double3::dot(dr_ab, dr_ab);
  double r_ab = std::sqrt(rr_ab);

  double3 dr_cb = posC - posB;
  double rr_cb = double3::dot(dr_cb, dr_cb);
  double r_cb = std::sqrt(rr_cb);

  switch (type)
  {
    case BondBondType::CVFF:
    case BondBondType::CFF:
      // p_0*(rab-p_1)*(rbc-p_2)
      // =======================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return parameters[0] * (r_ab - parameters[1]) * (r_cb - parameters[2]);
    default:
      std::unreachable();
  }
}

std::tuple<double, std::array<double3, 3>, double3x3> BondBondPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB, const double3 &posC) const
{
  double3x3 strain_derivative{};
  double3 Rab = posA - posB;
  double rab = std::sqrt(double3::dot(Rab, Rab));
  double3 Rbc = posC - posB;
  double rbc = std::sqrt(double3::dot(Rbc, Rbc));

  double U{};
  double3 du_da, du_db, du_dc;

  switch (type)
  {
    case BondBondType::CVFF:
    case BondBondType::CFF:
      U = parameters[0] * (rab - parameters[1]) * (rbc - parameters[2]);
      {
        const double gamma = -parameters[0] * (rbc - parameters[2]) / rab;
        const double gammc = -parameters[0] * (rab - parameters[1]) / rbc;
        du_da = {-gamma * Rab.x, -gamma * Rab.y, -gamma * Rab.z};
        du_dc = {-gammc * Rbc.x, -gammc * Rbc.y, -gammc * Rbc.z};
        du_db = -(du_da + du_dc);
      }
      break;
    default:
      std::unreachable();
  }

  strain_derivative.ax = Rab.x * du_da.x + Rbc.x * du_dc.x;
  strain_derivative.bx = Rab.y * du_da.x + Rbc.y * du_dc.x;
  strain_derivative.cx = Rab.z * du_da.x + Rbc.z * du_dc.x;
  strain_derivative.ay = Rab.x * du_da.y + Rbc.x * du_dc.y;
  strain_derivative.by = Rab.y * du_da.y + Rbc.y * du_dc.y;
  strain_derivative.cy = Rab.z * du_da.y + Rbc.z * du_dc.y;
  strain_derivative.az = Rab.x * du_da.z + Rbc.x * du_dc.z;
  strain_derivative.bz = Rab.y * du_da.z + Rbc.y * du_dc.z;
  strain_derivative.cz = Rab.z * du_da.z + Rbc.z * du_dc.z;

  return {U, {du_da, du_db, du_dc}, strain_derivative};
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondBondPotential &b)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondBondPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BondBondPotential' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("BondBondPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
