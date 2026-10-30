module;

module property_widom;

import archive;
import int3;
import units;
import averages;
import stringutils;
import json;

std::string PropertyWidom::writeAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                                   std::optional<double> frameworkMass,
                                                                   std::optional<int3> number_of_unit_cells) const
{
  std::ostringstream stream;

  std::pair<double, double> average_rosenbluth_weight = result();

  std::print(stream, "    Widom insertion Rosenbluth weight statistics:\n");
  std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    double blockAverage = averagedRosenbluthWeight(blockIndex);
    std::print(stream, "        Block[ {:2d}] {: .6e}\n", blockIndex, blockAverage);
  }
  std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
  std::print(stream, "    Average Rosenbluth weight:   {: .6e} +/- {: .6e} [-]\n", 
      average_rosenbluth_weight.first, average_rosenbluth_weight.second);
  std::print(stream, "\n\n");

  if (frameworkMass.has_value())
  {
    double frameworkDensity =
        1e-3 * frameworkMass.value() /
        (volume * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit * Units::AvogadroConstant);
    double conversion_factor_mol_per_kg = 1.0 / (Units::MolarGasConstant * temperature * frameworkDensity);

    std::print(stream, "    Henry coefficient based on Rosenbluth weight:\n");
    std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
    for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      double blockAverage = conversion_factor_mol_per_kg * averagedRosenbluthWeight(blockIndex);
      std::print(stream, "        Block[ {:2d}] {: .6e}\n", blockIndex, blockAverage);
    }
    std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
    std::print(stream, "    Average Henry coefficient:   {: .6e} +/- {: .6e} [mol/kg/Pa]\n",
               average_rosenbluth_weight.first * conversion_factor_mol_per_kg,
               average_rosenbluth_weight.second * conversion_factor_mol_per_kg);
    if (number_of_unit_cells.has_value())
    {
      double conversion_factor_molecules_per_uc =
          Units::AvogadroConstant * volume * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit /
          (Units::MolarGasConstant * temperature *
           static_cast<double>(number_of_unit_cells->x * number_of_unit_cells->y * number_of_unit_cells->z));

      std::print(stream, "    Average Henry coefficient:   {: .6e} +/- {: .6e} [molec./uc/Pa]\n",
                 average_rosenbluth_weight.first * conversion_factor_molecules_per_uc,
                 average_rosenbluth_weight.second * conversion_factor_molecules_per_uc);
    }
    std::print(stream, "\n\n");
  }

  return stream.str();
}

std::string PropertyWidom::writeAveragesChemicalPotentialStatistics(double beta,
                                                                    std::optional<double> imposedChemicalPotential,
                                                                    std::optional<double> imposedFugacity) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;

  std::pair<WidomData, WidomData> average_chemical_potential = chemicalPotentialResult(beta);
  std::pair<double, double> average_fugacity = fugacityResult(beta);

  switch (Units::unitSystem)
  {
    case Units::System::RASPA:
    {
      std::print(stream, "    Widom insertion chemical potential statistics:\n");
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedChemicalPotential(blockIndex, beta).excess;
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, conv * blockAverage);
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");

      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_chemical_potential.first.excess,
                 Units::EnergyToKelvin * average_chemical_potential.second.excess);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_chemical_potential.first.idealGas,
                 Units::EnergyToKelvin * average_chemical_potential.second.idealGas);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_chemical_potential.first.total,
                 Units::EnergyToKelvin * average_chemical_potential.second.total);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:         {: .6e} [K]\n",
                   Units::EnergyToKelvin * imposedChemicalPotential.value());
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_chemical_potential.first.excess,
                 Units::EnergyToKJPerMol * average_chemical_potential.second.excess);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_chemical_potential.first.idealGas,
                 Units::EnergyToKJPerMol * average_chemical_potential.second.idealGas);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_chemical_potential.first.total,
                 Units::EnergyToKJPerMol * average_chemical_potential.second.total);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:         {: .6e} [kJ/mol]\n",
                   Units::EnergyToKJPerMol * imposedChemicalPotential.value());
      }
      std::print(stream, "\n");

      std::print(stream, "    Widom insertion fugacity statistics:\n");
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedFugacity(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, Units::PressureConversionFactor * blockAverage);
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      std::print(stream, "    Measured fugacity:          {: .6e} +/- {: .6e} [Pa]\n",
                 Units::PressureConversionFactor * average_fugacity.first,
                 Units::PressureConversionFactor * average_fugacity.second);
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:         {: .6e} [K]\n",
                   Units::PressureConversionFactor * imposedFugacity.value());
      }
    }
    break;
    case Units::System::ReducedUnits:
    {
      std::print(stream, "    Widom insertion chemical potential statistics:\n");
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedChemicalPotential(blockIndex, beta).excess;
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, beta * blockAverage);
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      std::print(stream, "    Beta * Excess chemical potential:          {: .6e} +/- {: .6e} [-]\n",
                 beta * average_chemical_potential.first.excess,
                 beta * average_chemical_potential.second.excess);
      std::print(stream, "    Beta * Ideal chemical potential:           {: .6e} +/- {: .6e} [-]\n",
                 beta * average_chemical_potential.first.idealGas,
                 beta * average_chemical_potential.second.idealGas);
      std::print(stream, "    Beta * Total chemical potential:           {: .6e} +/- {: .6e} [-]\n",
                 beta * average_chemical_potential.first.total,
                 beta * average_chemical_potential.second.total);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Beta * Imposed chemical potential:  {: .6e} [-]\n",
                   beta * imposedChemicalPotential.value());
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [{}]\n",
                 average_chemical_potential.first.excess, average_chemical_potential.second.excess,
                 Units::unitOfEnergyString);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [{}]\n",
                 average_chemical_potential.first.idealGas, average_chemical_potential.second.idealGas,
                 Units::unitOfEnergyString);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [{}]\n",
                 average_chemical_potential.first.total, average_chemical_potential.second.total,
                 Units::unitOfEnergyString);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:  {: .6e} [{}]\n", imposedChemicalPotential.value(),
                   Units::unitOfEnergyString);
      }
      std::print(stream, "\n");

      std::print(stream, "    Widom insertion fugacity statistics:\n");
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedFugacity(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, blockAverage);
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
      std::print(stream, "    Total fugacity:           {: .6e} +/- {: .6e} [{}]\n",
                 Units::PressureConversionFactor * average_fugacity.first, 
                 Units::PressureConversionFactor * average_fugacity.second,
                 Units::unitOfPressureString);
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:  {: .6e} [{}]\n",
                   Units::PressureConversionFactor * imposedFugacity.value(), Units::unitOfPressureString);
      }
      break;
    }
  }

  std::print(stream, "\n\n");

  return stream.str();
}

nlohmann::json PropertyWidom::jsonAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                                     std::optional<double> frameworkMass,
                                                                     std::optional<int3> number_of_unit_cells) const
{
  nlohmann::json status;

  std::pair<double, double> average_rosenbluth_weight = result();

  std::vector<double> blocks(numberOfBlocks);
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    blocks[blockIndex] = averagedRosenbluthWeight(blockIndex);
  }
  status["rosenbluthWeight"]["blocks"]["[-]"] = blocks;
  status["rosenbluthWeight"]["mean"]["[-]"] = average_rosenbluth_weight.first;
  status["rosenbluthWeight"]["confidence"]["[-]"] = average_rosenbluth_weight.second;

  if (frameworkMass.has_value())
  {
    double frameworkDensity =
        1e-3 * frameworkMass.value() /
        (volume * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit * Units::AvogadroConstant);
    double conversion_factor_mol_per_kg = 1.0 / (Units::MolarGasConstant * temperature * frameworkDensity);

    std::vector<double> henryBlocks(numberOfBlocks);
    for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      henryBlocks[blockIndex] = conversion_factor_mol_per_kg * averagedRosenbluthWeight(blockIndex);
    }
    status["henryCoefficient"]["blocks"]["[mol/kg/Pa]"] = henryBlocks;
    status["henryCoefficient"]["mean"]["[mol/kg/Pa]"] = average_rosenbluth_weight.first * conversion_factor_mol_per_kg;
    status["henryCoefficient"]["confidence"]["[mol/kg/Pa]"] =
        average_rosenbluth_weight.second * conversion_factor_mol_per_kg;

    if (number_of_unit_cells.has_value())
    {
      double conversion_factor_molecules_per_uc =
          Units::AvogadroConstant * volume * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit /
          (Units::MolarGasConstant * temperature *
           static_cast<double>(number_of_unit_cells->x * number_of_unit_cells->y * number_of_unit_cells->z));

      status["henryCoefficient"]["mean"]["[molec./uc/Pa]"] =
          average_rosenbluth_weight.first * conversion_factor_molecules_per_uc;
      status["henryCoefficient"]["confidence"]["[molec./uc/Pa]"] =
          average_rosenbluth_weight.second * conversion_factor_molecules_per_uc;
    }
  }

  return status;
}

nlohmann::json PropertyWidom::jsonAveragesChemicalPotentialStatistics(double beta,
                                                                      std::optional<double> imposedChemicalPotential,
                                                                      std::optional<double> imposedFugacity) const
{
  nlohmann::json status;

  std::pair<WidomData, WidomData> average_chemical_potential = chemicalPotentialResult(beta);
  std::pair<double, double> average_fugacity = fugacityResult(beta);

  std::vector<double> excessBlocksKelvin(numberOfBlocks);
  std::vector<double> fugacityBlocks(numberOfBlocks);
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    excessBlocksKelvin[blockIndex] = Units::EnergyToKelvin * averagedChemicalPotential(blockIndex, beta).excess;
    fugacityBlocks[blockIndex] = Units::PressureConversionFactor * averagedFugacity(blockIndex, beta);
  }

  status["excessChemicalPotential"]["blocks"]["[K]"] = excessBlocksKelvin;
  status["excessChemicalPotential"]["mean"]["[K]"] = Units::EnergyToKelvin * average_chemical_potential.first.excess;
  status["excessChemicalPotential"]["confidence"]["[K]"] =
      Units::EnergyToKelvin * average_chemical_potential.second.excess;
  status["idealGasChemicalPotential"]["mean"]["[K]"] =
      Units::EnergyToKelvin * average_chemical_potential.first.idealGas;
  status["idealGasChemicalPotential"]["confidence"]["[K]"] =
      Units::EnergyToKelvin * average_chemical_potential.second.idealGas;
  status["totalChemicalPotential"]["mean"]["[K]"] = Units::EnergyToKelvin * average_chemical_potential.first.total;
  status["totalChemicalPotential"]["confidence"]["[K]"] =
      Units::EnergyToKelvin * average_chemical_potential.second.total;

  status["excessChemicalPotential"]["mean"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * average_chemical_potential.first.excess;
  status["excessChemicalPotential"]["confidence"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * average_chemical_potential.second.excess;
  status["idealGasChemicalPotential"]["mean"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * average_chemical_potential.first.idealGas;
  status["idealGasChemicalPotential"]["confidence"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * average_chemical_potential.second.idealGas;
  status["totalChemicalPotential"]["mean"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * average_chemical_potential.first.total;
  status["totalChemicalPotential"]["confidence"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * average_chemical_potential.second.total;

  if (imposedChemicalPotential)
  {
    status["imposedChemicalPotential"]["[K]"] = Units::EnergyToKelvin * imposedChemicalPotential.value();
    status["imposedChemicalPotential"]["[kJ/mol]"] = Units::EnergyToKJPerMol * imposedChemicalPotential.value();
  }

  status["measuredFugacity"]["blocks"] = fugacityBlocks;
  status["measuredFugacity"]["mean"] = Units::PressureConversionFactor * average_fugacity.first;
  status["measuredFugacity"]["confidence"] = Units::PressureConversionFactor * average_fugacity.second;
  if (imposedFugacity)
  {
    status["imposedFugacity"] = Units::PressureConversionFactor * imposedFugacity.value();
  }

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w)
{
  archive << w.versionNumber;

  archive << w.numberOfBlocks;
  archive << w.rosenbluthWeight;
  archive << w.chemicalPotentialTerms;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > w.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyWidom' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> w.numberOfBlocks;
  archive >> w.rosenbluthWeight;
  archive >> w.chemicalPotentialTerms;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyWidom: Error in binary restart\n"));
  }
#endif

  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const WidomData &l)
{
  archive << l.versionNumber;

  archive << l.total;
  archive << l.excess;
  archive << l.idealGas;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, WidomData &l)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > l.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'WidomData' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> l.total;
  archive >> l.excess;
  archive >> l.idealGas;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("WidomData: Error in binary restart\n"));
  }
#endif

  return archive;
}
