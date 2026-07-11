module;

module partial_molar_properties_data;

import std;

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PartialMolarPropertiesData &p)
{
  archive << p.size;
  archive << p.partialMolarEnergy;
  archive << p.partialMolarVolume;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PartialMolarPropertiesData &p)
{
  archive >> p.size;
  archive >> p.partialMolarEnergy;
  archive >> p.partialMolarVolume;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PartialMolarPropertiesData: Error in binary restart\n"));
  }
#endif

  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PartialMolarPropertiesTerms &p)
{
  archive << p.size;
  archive << p.swappableComponents;
  archive << p.totalEnergyTimesNumberOfMolecules;
  archive << p.volumeTimesNumberOfMolecules;
  archive << p.numberOfMoleculesSquared;
  archive << p.numberOfMolecules;
  archive << p.totalEnergy;
  archive << p.volume;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PartialMolarPropertiesTerms &p)
{
  archive >> p.size;
  archive >> p.swappableComponents;
  archive >> p.totalEnergyTimesNumberOfMolecules;
  archive >> p.volumeTimesNumberOfMolecules;
  archive >> p.numberOfMoleculesSquared;
  archive >> p.numberOfMolecules;
  archive >> p.totalEnergy;
  archive >> p.volume;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PartialMolarPropertiesTerms: Error in binary restart\n"));
  }
#endif

  return archive;
}
