module;

module bond_torsion_potential;

import std;

import archive;
import randomnumbers;
import double3;

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
