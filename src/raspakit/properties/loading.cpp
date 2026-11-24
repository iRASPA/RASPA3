module;

module property_loading;

import std;

import archive;
import int3;
import stringutils;
import units;
import component;
import json;

std::string PropertyLoading::writeAveragesStatistics(std::vector<Component> components,
                                                     std::optional<double> frameworkMass,
                                                     std::optional<int3> numberOfUnitCells) const
{
  std::ostringstream stream;

  if (frameworkMass.has_value())
  {
    const double toMolePerKg = 1000.0 / frameworkMass.value();

    std::pair<LoadingData, LoadingData> loadingAverage = average();

    int3 number_of_unit_cells = numberOfUnitCells.value_or(int3{1, 1, 1});
    double to_molecules_per_unit_cell =
        1.0 / (static_cast<double>(number_of_unit_cells.x * number_of_unit_cells.y * number_of_unit_cells.z));

    std::print(stream, "LoadingData\n");
    std::print(stream, "========================================================================================================================\n\n");

    for (std::size_t i = 0; i < components.size(); ++i)
    {
      const double toMgPerG = 1000.0 * components[i].totalMass / frameworkMass.value();

      std::print(stream, "Component {} ({})\n", i, components[i].name);

      for (std::size_t j = 0; j < numberOfBlocks; ++j)
      {
        LoadingData blockAverage = averaged(j);
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", j, blockAverage.numberOfMolecules[i]);
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");

      switch (Units::unitSystem)
      {
        case Units::System::RASPA:
          std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                     loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [molecules/uc]\n",
                     to_molecules_per_unit_cell * loadingAverage.first.numberOfMolecules[i],
                     to_molecules_per_unit_cell * loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [mol/kg-framework]\n",
                     toMolePerKg * loadingAverage.first.numberOfMolecules[i],
                     toMolePerKg * loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [mg/g-framework]\n",
                     toMgPerG * loadingAverage.first.numberOfMolecules[i],
                     toMgPerG * loadingAverage.second.numberOfMolecules[i]);
          break;
        case Units::System::ReducedUnits:
          std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                     loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
          break;
      }

      std::print(stream, "\n");

      for (std::size_t j = 0; j < numberOfBlocks; ++j)
      {
        LoadingData blockAverage = averaged(j);
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", j,
                   blockAverage.numberOfMolecules[i] - components[i].amountOfExcessMolecules);
      }
      std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");

      switch (Units::unitSystem)
      {
        case Units::System::RASPA:
          std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                     loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules,
                     loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [molecules/uc]\n",
                     to_molecules_per_unit_cell *
                         (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules),
                     to_molecules_per_unit_cell * loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [mol/kg-framework]\n",
                     toMolePerKg * (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules),
                     toMolePerKg * loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [mg/g-framework]\n",
                     toMgPerG * (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules),
                     toMgPerG * loadingAverage.second.numberOfMolecules[i]);
          break;
        case Units::System::ReducedUnits:
          std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                     loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules,
                     loadingAverage.second.numberOfMolecules[i]);
          break;
      }

      std::print(stream, "\n\n");
    }
  }
  else
  {
    const double densityConversionFactor =
        1.0 / (1000.0 * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::pair<LoadingData, LoadingData> loadingAverage = average();

    std::print(stream, "Densities\n");
    std::print(stream, "========================================================================================================================\n\n");

    for (std::size_t i = 0; i < components.size(); ++i)
    {
      std::print(stream, "Component {} ({})\n", i, components[i].name);

      for (std::size_t j = 0; j < numberOfBlocks; ++j)
      {
        LoadingData blockAverage = averaged(j);
        std::print(stream, "    Block[ {:2d}] {}\n", j, blockAverage.numberOfMolecules[i]);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");

      switch (Units::unitSystem)
      {
        case Units::System::RASPA:
          std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molecules]\n",
                     loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molec/A^3]\n",
                     loadingAverage.first.numberDensities[i], loadingAverage.second.numberDensities[i]);
          std::print(stream, "    Density average  {: .6e} +/- {: .6e} [{}]\n",
                     densityConversionFactor * components[i].totalMass * loadingAverage.first.numberDensities[i],
                     densityConversionFactor * components[i].totalMass * loadingAverage.second.numberDensities[i],
                     Units::unitOfDensityString);
          if (loadingAverage.first.totalNumberOfMolecules > 0.0)
          {
            // (the error of the ratio over the blocks: the error of N_i over the mean of N gave a single
            // component, whose mole fraction is one, an error)
            const std::pair<double, double> moleFraction = averageMoleFraction(i);
            std::print(stream, "    Mol-fraction average {: .6e} +/- {: .6e} [-]\n", moleFraction.first,
                       moleFraction.second);
          }
          break;
        case Units::System::ReducedUnits:
          std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molecules]\n",
                     loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
          std::print(stream, "    Density average  {: .6e} +/- {: .6e} [{}]\n", loadingAverage.first.numberDensities[i],
                     loadingAverage.second.numberDensities[i], Units::unitOfDensityString);
          if (loadingAverage.first.totalNumberOfMolecules > 0.0)
          {
            // (the error of the ratio over the blocks: the error of N_i over the mean of N gave a single
            // component, whose mole fraction is one, an error)
            const std::pair<double, double> moleFraction = averageMoleFraction(i);
            std::print(stream, "    Mol-fraction average {: .6e} +/- {: .6e} [-]\n", moleFraction.first,
                       moleFraction.second);
          }
          break;
      }
    }
  }

  std::print(stream, "\n");

  return stream.str();
}

std::pair<double, double> PropertyLoading::averageMoleFraction(std::size_t comp) const
{
  return statistics(
      [comp](const LoadingData &data)
      { return data.totalNumberOfMolecules > 0.0 ? data.numberOfMolecules[comp] / data.totalNumberOfMolecules : 0.0; });
}

// The loadings (with a framework) or the densities (without one) of the text statistics, per component.
nlohmann::json PropertyLoading::jsonAveragesStatistics(const std::vector<Component> &components,
                                                       std::optional<double> frameworkMass,
                                                       std::optional<int3> numberOfUnitCells) const
{
  nlohmann::json status;
  const std::pair<LoadingData, LoadingData> loadingAverage = average();
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;

  for (std::size_t i = 0; i < components.size(); ++i)
  {
    nlohmann::json component;
    std::vector<double> blocks(numberOfBlocks);
    for (std::size_t j = 0; j < numberOfBlocks; ++j) blocks[j] = averaged(j).numberOfMolecules[i];
    const double mean = loadingAverage.first.numberOfMolecules[i];
    const double error = loadingAverage.second.numberOfMolecules[i];

    component["block"] = blocks;
    component["mean"]["[molecules/cell]"] = mean;
    component["confidence"]["[molecules/cell]"] = error;

    if (frameworkMass.has_value())
    {
      const double excess = components[i].amountOfExcessMolecules;
      component["excess"]["mean"]["[molecules/cell]"] = mean - excess;
      component["excess"]["confidence"]["[molecules/cell]"] = error;
      if (!reduced)
      {
        const int3 cells = numberOfUnitCells.value_or(int3{1, 1, 1});
        const double toMoleculesPerUnitCell = 1.0 / static_cast<double>(cells.x * cells.y * cells.z);
        const double toMolePerKg = 1000.0 / frameworkMass.value();
        const double toMgPerG = 1000.0 * components[i].totalMass / frameworkMass.value();
        const std::array<std::pair<const char *, double>, 3> units{
            std::pair{"[molecules/uc]", toMoleculesPerUnitCell}, std::pair{"[mol/kg-framework]", toMolePerKg},
            std::pair{"[mg/g-framework]", toMgPerG}};
        for (const auto &[unit, factor] : units)
        {
          component["mean"][unit] = factor * mean;
          component["confidence"][unit] = factor * error;
          component["excess"]["mean"][unit] = factor * (mean - excess);
          component["excess"]["confidence"][unit] = factor * error;
        }
      }
    }
    else
    {
      const double densityConversionFactor =
          reduced ? 1.0 : components[i].totalMass / (1000.0 * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                                                    Units::AvogadroConstant);
      const std::string densityUnit = std::format("[{}]", Units::unitOfDensityString);
      component["density"]["mean"][densityUnit] = densityConversionFactor * loadingAverage.first.numberDensities[i];
      component["density"]["confidence"][densityUnit] =
          densityConversionFactor * loadingAverage.second.numberDensities[i];
      if (loadingAverage.first.totalNumberOfMolecules > 0.0)
      {
        const std::pair<double, double> moleFraction = averageMoleFraction(i);
        component["moleFraction"]["mean"] = moleFraction.first;
        component["moleFraction"]["confidence"] = moleFraction.second;
      }
    }
    status[components[i].name] = component;
  }
  return status;
}

std::pair<double, double> PropertyLoading::averageLoadingNumberOfMolecules(std::size_t comp) const
{
  std::pair<LoadingData, LoadingData> loadingAverage = average();
  if (comp >= loadingAverage.first.numberOfMolecules.size())
    throw std::out_of_range(std::format("PropertyLoading: component {} out of range ({})", comp,
                                        loadingAverage.first.numberOfMolecules.size()));
  return {loadingAverage.first.numberOfMolecules[comp], loadingAverage.second.numberOfMolecules[comp]};
}

std::string PropertyLoading::repr() const { return std::string("PropertyLoading test"); }

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
        if (totalNumberOfMolecules > 0.0)
        {
          std::print(stream, "    mol-fraction:     {: .6e} [-]\n",
                     numberOfMolecules[componentId] / totalNumberOfMolecules);
        }
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    molecules:        {: .6e} molecules\n", numberOfMolecules[componentId]);
        std::print(stream, "    number density:   {: .6e} molec./{}^3\n", numberDensities[componentId],
                   Units::displayedUnitOfLengthString);
        if (totalNumberOfMolecules > 0.0)
        {
          std::print(stream, "    mol-fraction:     {: .6e} [-]\n",
                     numberOfMolecules[componentId] / totalNumberOfMolecules);
        }
        break;
    }
  }

  return stream.str();
}

std::string LoadingData::printStatus(std::size_t componentId, const std::string &componentName, double componentTotalMass,
                                  double componentAmountOfExcessMolecules, const LoadingData &average, const LoadingData &error,
                                  std::optional<double> frameworkMass, std::optional<int3> numberOfUnitCells,
                                  std::optional<std::pair<double, double>> moleFraction) const
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
        if (totalNumberOfMolecules > 0.0 && average.totalNumberOfMolecules > 0.0)
        {
          const double molFraction = numberOfMolecules[componentId] / totalNumberOfMolecules;
          // The mole fraction and its error from the mole fractions of the blocks, where the caller has them. The
          // error of the number of molecules over the total says nothing of the fraction: a pure fluid, whose
          // fraction is one in every block, was given the relative error of its loading.
          const double averageMolFraction =
              moleFraction.has_value() ? moleFraction->first
                                       : average.numberOfMolecules[componentId] / average.totalNumberOfMolecules;
          const double errorMolFraction =
              moleFraction.has_value() ? moleFraction->second
                                       : error.numberOfMolecules[componentId] / average.totalNumberOfMolecules;
          std::print(stream, "    mol-fraction:   {:.6e} [-]        ({:.6e} +/- {:.6e})\n", molFraction,
                     averageMolFraction, errorMolFraction);
        }
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    molecules:      {:.6e} molecules  ({:.6e} +/- {:.6e})\n",
                   numberOfMolecules[componentId], average.numberOfMolecules[componentId],
                   error.numberOfMolecules[componentId]);
        std::print(stream, "    number density: {:.6e} molec./{}^3 ({:.6e} +/- {:.6e})\n",
                   numberDensities[componentId], Units::displayedUnitOfLengthString,
                   average.numberDensities[componentId], error.numberDensities[componentId]);
        if (totalNumberOfMolecules > 0.0 && average.totalNumberOfMolecules > 0.0)
        {
          const double molFraction = numberOfMolecules[componentId] / totalNumberOfMolecules;
          // The mole fraction and its error from the mole fractions of the blocks, where the caller has them. The
          // error of the number of molecules over the total says nothing of the fraction: a pure fluid, whose
          // fraction is one in every block, was given the relative error of its loading.
          const double averageMolFraction =
              moleFraction.has_value() ? moleFraction->first
                                       : average.numberOfMolecules[componentId] / average.totalNumberOfMolecules;
          const double errorMolFraction =
              moleFraction.has_value() ? moleFraction->second
                                       : error.numberOfMolecules[componentId] / average.totalNumberOfMolecules;
          std::print(stream, "    mol-fraction:   {:.6e} [-]        ({:.6e} +/- {:.6e})\n", molFraction,
                     averageMolFraction, errorMolFraction);
        }
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
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
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
