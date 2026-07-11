module;

module property_enthalpy;

import std;

import archive;
import stringutils;
import component;
import units;
import averages;
import json;
import averages;

std::string PropertyEnthalpy::writeAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                                      std::vector<Component> &components) const
{
  std::ostringstream stream;

  std::print(stream, "Enthalpy of adsorption\n");
  std::print(stream, "===============================================================================\n\n");

  if (swappableComponents.empty())
  {
    std::print(stream, "No fluctuating components present.\n\n");
  }
  else
  {
    std::pair<EnthalpyOfAdsorptionData, EnthalpyOfAdsorptionData> enthalpy = average();
    for (std::size_t k = 0; k < swappableComponents.size(); k++)
    {
      std::size_t index = swappableComponents[k];
      double idealGasTerm = components[index].idealGasEnergy.value_or(0.0);
      std::print(stream, "Component {} [{}]\n", index, components[index].name);
      std::print(stream, "-------------------------------------------------------------------------------\n");
      for (std::size_t i = 0; i < numberOfBlocks; ++i)
      {
        EnthalpyOfAdsorptionData average = averaged(i);
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i,
                   Units::EnergyToKelvin * (average.values[k] - idealGasTerm));
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Enthalpy of adsorption: {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * (enthalpy.first.values[k] - idealGasTerm),
                 Units::EnergyToKelvin * enthalpy.second.values[k]);
      std::print(stream, "                            {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * (enthalpy.first.values[k] - idealGasTerm),
                 Units::EnergyToKJPerMol * enthalpy.second.values[k]);
      if (!components[index].idealGasEnergy)
      {
        std::print(stream, "    Note: need to subtract the ideal-gas energy.\n");
      }
      std::print(stream, "\n");
    }
    if (swappableComponents.size() > 1)
    {
      std::print(stream, "Total enthalpy of adsorption\n");
      std::print(stream, "-------------------------------------------------------------------------------\n");

      std::vector<double> totalEnthalpyBlocks(numberOfBlocks);
      for (std::size_t i = 0; i < numberOfBlocks; ++i)
      {
        double totalEnthalpyOfAdsorption = 0.0;
        EnthalpyOfAdsorptionData average = averaged(i);
        for (std::size_t k = 0; k < swappableComponents.size(); k++)
        {
          std::size_t index = swappableComponents[k];
          double idealGasTerm = components[index].idealGasEnergy.value_or(0.0);
          totalEnthalpyOfAdsorption += components[index].molFraction * (average.values[k] - idealGasTerm);
        }
        totalEnthalpyBlocks[i] = totalEnthalpyOfAdsorption;
        std::print(stream, "    Block[ {:2d}] {}\n", i, Units::EnergyToKelvin * totalEnthalpyOfAdsorption);
      }
      std::pair<double, double> totalEnthalpy = meanConfidence(totalEnthalpyBlocks);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Enthalpy of adsorption: {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * totalEnthalpy.first, Units::EnergyToKelvin * totalEnthalpy.second);
      std::print(stream, "                            {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * totalEnthalpy.first, Units::EnergyToKJPerMol * totalEnthalpy.second);
      for (std::size_t k = 0; k < swappableComponents.size(); k++)
      {
        std::size_t index = swappableComponents[k];
        if (!components[index].idealGasEnergy)
        {
          std::print(stream, "    Note: Recompute value of the total enthalpy of adsorption by hand.\n");
          std::print(stream, "          Need to subtract the ideal-gas energy of component {}.\n", index);
        }
      }
      std::print(stream, "\n");
    }
  }

  std::print(stream, "\n");

  return stream.str();
}

nlohmann::json PropertyEnthalpy::jsonAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                                        std::vector<Component> &components) const
{
  nlohmann::json status;

  if (!swappableComponents.empty())
  {
    std::pair<EnthalpyOfAdsorptionData, EnthalpyOfAdsorptionData> enthalpy = average();
    for (std::size_t k = 0; k < swappableComponents.size(); k++)
    {
      std::size_t index = swappableComponents[k];
      double idealGasTerm = components[index].idealGasEnergy.value_or(0.0);

      std::vector<double> blockEnthalpy = blockEnthalpies(k, idealGasTerm);
      status[components[index].name]["block"] = blockEnthalpy;
      status[components[index].name]["mean"]["[K]"] = Units::EnergyToKelvin * (enthalpy.first.values[k] - idealGasTerm);
      status[components[index].name]["confidence"]["[K]"] = Units::EnergyToKelvin * enthalpy.second.values[k];
      status[components[index].name]["mean"]["[kJ/mol]"] =
          Units::EnergyToKJPerMol * (enthalpy.first.values[k] - idealGasTerm);
      status[components[index].name]["confidence"]["[kJ/mol]"] = Units::EnergyToKJPerMol * enthalpy.second.values[k];

      if (!components[index].idealGasEnergy)
      {
        status[components[index].name]["warning"] =
            "Note: Recompute value of the total enthalpy of adsorption by hand. Need to subtract the ideal-gas "
            "energy";
      }
    }
    if (swappableComponents.size() > 1)
    {
      std::vector<double> blockTotalEnthalpy(numberOfBlocks);

      for (std::size_t k = 0; k < swappableComponents.size(); k++)
      {
        std::size_t index = swappableComponents[k];
        double idealGasTerm = components[index].idealGasEnergy.value_or(0.0);
        std::vector<double> blockEnthalpy = blockEnthalpies(k, idealGasTerm);
        for (std::size_t i = 0; i < numberOfBlocks; i++)
        {
          blockTotalEnthalpy[i] += components[index].molFraction * blockEnthalpy[i];
        }
      }

      std::pair<double, double> totalEnthalpy = meanConfidence(blockTotalEnthalpy);
      status["total"]["block"] = blockTotalEnthalpy;
      status["total"]["mean"]["[K]"] = totalEnthalpy.first;
      status["total"]["confidence"]["[K]"] = totalEnthalpy.second;
      status["total"]["mean"]["[kJ/mol]"] = totalEnthalpy.first;
      status["total"]["confidence"]["[kJ/mol]"] = totalEnthalpy.second;
    }
  }

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnthalpyOfAdsorptionData &p)
{
  archive << p.size;
  archive << p.values;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnthalpyOfAdsorptionData &p)
{
  archive >> p.size;
  archive >> p.values;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyEnthalpy: Error in binary restart\n"));
  }
#endif

  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnthalpyOfAdsorptionTerms &p)
{
  archive << p.size;
  archive << p.swappableComponents;
  archive << p.totalEnergyTimesNumberOfMolecules;
  archive << p.numberOfMoleculesSquared;
  archive << p.numberOfMolecules;
  archive << p.temperature;
  archive << p.totalEnergy;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnthalpyOfAdsorptionTerms &p)
{
  archive >> p.size;
  archive >> p.swappableComponents;
  archive >> p.totalEnergyTimesNumberOfMolecules;
  archive >> p.numberOfMoleculesSquared;
  archive >> p.numberOfMolecules;
  archive >> p.temperature;
  archive >> p.totalEnergy;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("EnthalpyOfAdsorptionTerms: Error in binary restart\n"));
  }
#endif

  return archive;
}
