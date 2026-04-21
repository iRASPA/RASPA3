module;

module loading_data;

import std;

import archive;
import int3;
import stringutils;
import component;
import units;

std::string LoadingData::printStatus(std::size_t componentId, const std::string &componentName, double componentTotalMass,
                                  double componentAmountOfExcessMolecules, std::optional<double> frameworkMass,
                                  std::optional<int3> numberOfUnitCells) const
{
  std::ostringstream stream;

  if (frameworkMass.has_value())
  {
    std::print(stream, "Component {} ({})\n", componentId, componentName);

    const double toMolePerKg = 1000.0 / frameworkMass.value();
    const double toMgPerG = 1000.0 * componentTotalMass / frameworkMass.value();

    int3 number_of_unit_cells = numberOfUnitCells.value_or(int3{1, 1, 1});
    double to_molecules_per_unit_cell =
        1.0 / (static_cast<double>(number_of_unit_cells.x * number_of_unit_cells.y * number_of_unit_cells.z));

    double loading = numberOfMolecules[componentId];
    double excess_loading = numberOfMolecules[componentId] - componentAmountOfExcessMolecules;
    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
        std::print(stream, "    absolute adsorption: {: .6e} molecules\n", loading);
        std::print(stream, "                         {: .6e} molecules/uc\n", loading * to_molecules_per_unit_cell);
        std::print(stream, "                         {: .6e} mol/kg-framework\n", loading * toMolePerKg);
        std::print(stream, "                         {: .6e} mg/g-framework\n", loading * toMgPerG);

        std::print(stream, "    excess adsorption:   {: .6e} molecules\n", excess_loading);
        std::print(stream, "                         {: .6e} molecules/uc\n",
                   excess_loading * to_molecules_per_unit_cell);
        std::print(stream, "                         {: .6e} mol/kg-framework\n", excess_loading * toMolePerKg);
        std::print(stream, "                         {: .6e} mg/g-framework\n", excess_loading * toMgPerG);
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    absolute adsorption: {: .6e} molecules\n", loading);
        std::print(stream, "    excess adsorption:   {: .6e} molecules\n", excess_loading);
        break;
    }
  }
  else
  {
    const double densityConversionFactor =
        1.0 / (1000.0 * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::print(stream, "Component {} ({})\n", componentId, componentName);
    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
        std::print(stream, "    molecules:        {: .6e} molecules\n", numberOfMolecules[componentId]);
        std::print(stream, "    number density:   {: .6e} molec/A^3\n", numberDensities[componentId]);
        std::print(stream, "    density:          {: .6e} kg/m^3\n",
                   densityConversionFactor * componentTotalMass * numberDensities[componentId]);
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    molecules:        {: .6e} molecules\n", numberOfMolecules[componentId]);
        std::print(stream, "    number density:   {: .6e} molec./{}^3\n", numberDensities[componentId],
                   Units::displayedUnitOfLengthString);
        break;
    }
  }

  return stream.str();
}

std::string LoadingData::printStatus(std::size_t componentId, const std::string &componentName, double componentTotalMass,
                                  double componentAmountOfExcessMolecules, const LoadingData &average, const LoadingData &error,
                                  std::optional<double> frameworkMass, std::optional<int3> numberOfUnitCells) const
{
  std::ostringstream stream;

  if (frameworkMass)
  {
    std::print(stream, "Component {} ({})\n", componentId, componentName);

    const double toMolePerKg = 1000.0 / frameworkMass.value();
    const double toMgPerKg = 1000.0 * componentTotalMass / frameworkMass.value();

    int3 number_of_unit_cells = numberOfUnitCells.value_or(int3{1, 1, 1});
    double to_molecules_per_unit_cell =
        1.0 / (static_cast<double>(number_of_unit_cells.x * number_of_unit_cells.y * number_of_unit_cells.z));

    double loading = numberOfMolecules[componentId];
    double loading_avg = average.numberOfMolecules[componentId];
    double loading_error = error.numberOfMolecules[componentId];

    double excess_loading = numberOfMolecules[componentId] - componentAmountOfExcessMolecules;
    double excess_loading_avg = average.numberOfMolecules[componentId] - componentAmountOfExcessMolecules;
    double excess_loading_error = error.numberOfMolecules[componentId];

    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
        std::print(stream, "    absolute adsorption: {:.6e} molecules ({:.6e} +/- {:.6e})\n", loading, loading_avg,
                   loading_error);
        std::print(stream, "                         {:.6e} molec./uc ({:.6e} +/- {:.6e})\n",
                   loading * to_molecules_per_unit_cell, loading_avg * to_molecules_per_unit_cell,
                   loading_error * to_molecules_per_unit_cell);
        std::print(stream, "                         {:.6e} mol/kg    ({:.6e} +/- {:.6e})\n", loading * toMolePerKg,
                   loading_avg * toMolePerKg, loading_error * toMolePerKg);
        std::print(stream, "                         {:.6e} mg/g      ({:.6e} +/- {:.6e})\n", loading * toMgPerKg,
                   loading_avg * toMgPerKg, loading_error * toMgPerKg);

        std::print(stream, "    excess adsorption:   {:.6e} molecules ({:.6e} +/- {:.6e})\n", excess_loading,
                   excess_loading_avg, excess_loading_error);
        std::print(stream, "                         {:.6e} molec./uc ({:.6e} +/- {:.6e})\n",
                   excess_loading * to_molecules_per_unit_cell, excess_loading_avg * to_molecules_per_unit_cell,
                   excess_loading_error * to_molecules_per_unit_cell);
        std::print(stream, "                         {:.6e} mol/kg    ({:.6e} +/- {:.6e})\n",
                   excess_loading * toMolePerKg, excess_loading_avg * toMolePerKg, excess_loading_error * toMolePerKg);
        std::print(stream, "                         {:.6e} mg/g      ({:.6e} +/- {:.6e})\n",
                   excess_loading * toMgPerKg, excess_loading_avg * toMgPerKg, excess_loading_error * toMgPerKg);
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    absolute adsorption: {:.6e} molecules ({:.6e} +/- {:.6e})\n", loading, loading_avg,
                   loading_error);
        std::print(stream, "    excess adsorption:   {:.6e} molecules ({:.6e} +/- {:.6e})\n", excess_loading,
                   excess_loading_avg, excess_loading_error);
        break;
    }
  }
  else
  {
    const double densityConversionFactor =
        1.0 / (1000.0 * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::print(stream, "Component {} ({})\n", componentId, componentName);
    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
        std::print(stream, "    molecules:      {:.6e} molecules  ({:.6e} +/- {:.6e})\n",
                   numberOfMolecules[componentId], average.numberOfMolecules[componentId],
                   error.numberOfMolecules[componentId]);
        std::print(stream, "    number density: {:.6e} molec./A^3 ({:.6e} +/- {:.6e})\n",
                   numberDensities[componentId], average.numberDensities[componentId],
                   error.numberDensities[componentId]);
        std::print(stream, "    density:        {:.6e} kg/m^3     ({:.6e} +/- {:.6e})\n",
                   densityConversionFactor * componentTotalMass * numberDensities[componentId],
                   densityConversionFactor * componentTotalMass * average.numberDensities[componentId],
                   densityConversionFactor * componentTotalMass * error.numberDensities[componentId]);
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    molecules:      {:.6e} molecules  ({:.6e} +/- {:.6e})\n",
                   numberOfMolecules[componentId], average.numberOfMolecules[componentId],
                   error.numberOfMolecules[componentId]);
        std::print(stream, "    number density: {:.6e} molec./{}^3 ({:.6e} +/- {:.6e})\n",
                   numberDensities[componentId], Units::displayedUnitOfLengthString,
                   average.numberDensities[componentId], error.numberDensities[componentId]);
        break;
    }
  }

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const LoadingData &l)
{
  archive << l.versionNumber;

  archive << l.size;
  archive << l.totalNumberOfMolecules;
  archive << l.totalDensity;
  archive << l.numberOfMolecules;
  archive << l.numberDensities;
  archive << l.inverseNumberDensities;

#if DEBUG_ARCHIVE
  archive << static_caststd::<int64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, LoadingData &l)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > l.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'LoadingData' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> l.size;
  archive >> l.totalNumberOfMolecules;
  archive >> l.totalDensity;
  archive >> l.numberOfMolecules;
  archive >> l.numberDensities;
  archive >> l.inverseNumberDensities;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("LoadingData: Error in binary restart\n"));
  }
#endif

  return archive;
}
