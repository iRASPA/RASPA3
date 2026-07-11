module;

module property_partial_molar_properties;

import std;

import archive;
import stringutils;
import component;
import units;
import averages;
import json;

std::string PropertyPartialMolarProperties::writeAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                                                    std::vector<Component> &components) const
{
  std::ostringstream stream;

  // Conversion from simulation volume units (Angstrom^3 per molecule) to cm^3/mol.
  const double volumeToCm3PerMol = Units::VolumeConversionFactor * Units::AvogadroConstant * 1.0e6;

  std::print(stream, "Partial molar properties\n");
  std::print(stream, "===============================================================================\n\n");

  if (swappableComponents.empty())
  {
    std::print(stream, "No fluctuating components present.\n\n");
  }
  else
  {
    std::pair<PartialMolarPropertiesData, PartialMolarPropertiesData> properties = average();
    for (std::size_t k = 0; k < swappableComponents.size(); k++)
    {
      std::size_t index = swappableComponents[k];
      std::print(stream, "Component {} [{}]\n", index, components[index].name);
      std::print(stream, "-------------------------------------------------------------------------------\n");

      std::print(stream, "  Partial molar internal energy:\n");
      for (std::size_t i = 0; i < numberOfBlocks; ++i)
      {
        PartialMolarPropertiesData average = averaged(i);
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, Units::EnergyToKelvin * average.partialMolarEnergy[k]);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Partial molar energy:   {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * properties.first.partialMolarEnergy[k],
                 Units::EnergyToKelvin * properties.second.partialMolarEnergy[k]);
      std::print(stream, "                            {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * properties.first.partialMolarEnergy[k],
                 Units::EnergyToKJPerMol * properties.second.partialMolarEnergy[k]);

      std::print(stream, "  Partial molar volume:\n");
      for (std::size_t i = 0; i < numberOfBlocks; ++i)
      {
        PartialMolarPropertiesData average = averaged(i);
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, average.partialMolarVolume[k]);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Partial molar volume:   {: .6e} +/- {: .6e} [A^3]\n",
                 properties.first.partialMolarVolume[k], properties.second.partialMolarVolume[k]);
      std::print(stream, "                            {: .6e} +/- {: .6e} [cm^3/mol]\n",
                 volumeToCm3PerMol * properties.first.partialMolarVolume[k],
                 volumeToCm3PerMol * properties.second.partialMolarVolume[k]);
      std::print(stream, "\n");
    }
  }

  std::print(stream, "\n");

  return stream.str();
}

nlohmann::json PropertyPartialMolarProperties::jsonAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                                                      std::vector<Component> &components) const
{
  nlohmann::json status;

  const double volumeToCm3PerMol = Units::VolumeConversionFactor * Units::AvogadroConstant * 1.0e6;

  if (!swappableComponents.empty())
  {
    std::pair<PartialMolarPropertiesData, PartialMolarPropertiesData> properties = average();
    for (std::size_t k = 0; k < swappableComponents.size(); k++)
    {
      std::size_t index = swappableComponents[k];

      std::vector<double> blockEnergy(numberOfBlocks);
      std::vector<double> blockVolume(numberOfBlocks);
      for (std::size_t i = 0; i < numberOfBlocks; ++i)
      {
        PartialMolarPropertiesData average = averaged(i);
        blockEnergy[i] = Units::EnergyToKelvin * average.partialMolarEnergy[k];
        blockVolume[i] = average.partialMolarVolume[k];
      }

      status[components[index].name]["partialMolarEnergy"]["block"] = blockEnergy;
      status[components[index].name]["partialMolarEnergy"]["mean"]["[K]"] =
          Units::EnergyToKelvin * properties.first.partialMolarEnergy[k];
      status[components[index].name]["partialMolarEnergy"]["confidence"]["[K]"] =
          Units::EnergyToKelvin * properties.second.partialMolarEnergy[k];
      status[components[index].name]["partialMolarEnergy"]["mean"]["[kJ/mol]"] =
          Units::EnergyToKJPerMol * properties.first.partialMolarEnergy[k];
      status[components[index].name]["partialMolarEnergy"]["confidence"]["[kJ/mol]"] =
          Units::EnergyToKJPerMol * properties.second.partialMolarEnergy[k];

      status[components[index].name]["partialMolarVolume"]["block"] = blockVolume;
      status[components[index].name]["partialMolarVolume"]["mean"]["[A^3]"] = properties.first.partialMolarVolume[k];
      status[components[index].name]["partialMolarVolume"]["confidence"]["[A^3]"] =
          properties.second.partialMolarVolume[k];
      status[components[index].name]["partialMolarVolume"]["mean"]["[cm^3/mol]"] =
          volumeToCm3PerMol * properties.first.partialMolarVolume[k];
      status[components[index].name]["partialMolarVolume"]["confidence"]["[cm^3/mol]"] =
          volumeToCm3PerMol * properties.second.partialMolarVolume[k];
    }
  }

  return status;
}

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
