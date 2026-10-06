module;

module van_der_waals_potential;

import std;

import archive;
import double3;
import units;
import vdwparameters;

std::string VanDerWaalsPotential::print() const
{
  const VDWParameters::ParameterMetadata metadata = VDWParameters::parameterMetadata(parameters.type);
  std::string parameterString{};
  for (std::size_t k = 0; k < metadata.count; ++k)
  {
    double value = k < 4 ? parameters.parameters[k] : parameters.parameters2[k - 4];
    if (metadata.isEnergy[k]) value *= Units::EnergyToKelvin;
    std::format_to(std::back_inserter(parameterString), "p_{}{}={:g}{}", k, metadata.isEnergy[k] ? "/k_B" : "",
                   value, k + 1 < metadata.count ? ", " : "");
  }
  return std::format("{} - {} : {} {}, scaling={:g}\n", identifiers[0], identifiers[1],
                     VDWParameters::nameOfType(parameters.type), parameterString, scaling);
}

double VanDerWaalsPotential::calculateEnergy(const double3 &posA, const double3 &posB) const
{
  const double3 dr = posA - posB;
  const double rr = double3::dot(dr, dr);
  return scaling * parameters.potentialEnergyAtFullCoupling(rr);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VanDerWaalsPotential &b)
{
  archive << b.versionNumber;

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
