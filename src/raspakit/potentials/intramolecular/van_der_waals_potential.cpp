module;

module van_der_waals_potential;

import std;

import archive;
import randomnumbers;
import double3;
import double3x3;
import units;

VanDerWaalsPotential::VanDerWaalsPotential(std::array<std::size_t, 2> identifiers, VanDerWaalsType type,
                                           std::vector<double> vector_parameters, double scaling)
    : identifiers(identifiers), type(type), scaling(scaling)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfVanDerWaalsParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      parameters[0] *= Units::KelvinToEnergy;
      break;
    default:
      std::unreachable();
  }
}

std::string VanDerWaalsPotential::print() const
{
  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      return std::format("{} - {} : LENNARD_JONES p_0/k_B={:g} [K], p_1={:g} [Å], scaling/k_B={:g}\n", identifiers[0],
                         identifiers[1], parameters[0] * Units::EnergyToKelvin, parameters[1], scaling);
    default:
      std::unreachable();
  }
}

double VanDerWaalsPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  double temp;
  double rri;

  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);

  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [Å]
      rri = (parameters[1] * parameters[1]) / rr;
      temp = rri * rri * rri;
      return scaling * 4.0 * parameters[0] * (temp * (temp - 1.0));
    default:
      std::unreachable();
  }
}

std::tuple<double, std::array<double3, 2>, double3x3> VanDerWaalsPotential::potentialEnergyGradientStrain(
    const double3 &posA, const double3 &posB) const
{
  double temp;
  double U{}, DF{};
  double3 du_dr, du_da, du_db;
  double3x3 strain_derivative{};

  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  switch (type)
  {
    case VanDerWaalsType::LennardJones:
      temp = (parameters[1] / rr) * (parameters[1] / rr) * (parameters[1] / rr);
      U = scaling * 4.0 * parameters[0] * (temp * (temp - 1.0));
      if (r > 0.0)
      {
        DF = scaling * 24.0 * parameters[0] * (temp * (1.0 - 2.0 * temp)) / rr;
      }
      break;
    default:
      std::unreachable();
  }

  du_dr = DF * dr;
  du_da = du_dr;
  du_db = -du_dr;

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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VanDerWaalsPotential &b)
{
  archive << b.versionNumber;

  archive << b.type;
  archive << b.identifiers;
  archive << b.scaling;
  archive << b.parameters;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VanDerWaalsPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'VanDerWaalsPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.type;
  archive >> b.identifiers;
  archive >> b.scaling;
  archive >> b.parameters;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("VanDerWaalsPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
