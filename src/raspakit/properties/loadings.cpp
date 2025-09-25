module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
#endif

module loadings;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import int3;
import stringutils;
import component;
import units;

std::string Loadings::printStatus(const Component &comp, std::optional<double> frameworkMass,
                                  std::optional<int3> numberOfUnitCells) const
{
  std::ostringstream stream;

  if (frameworkMass.has_value())
  {
    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);

    const double toMolePerKg = 1000.0 / frameworkMass.value();
    const double toMgPerG = 1000.0 * comp.totalMass / frameworkMass.value();

    int3 number_of_unit_cells = numberOfUnitCells.value_or(int3{1, 1, 1});
    double to_molecules_per_unit_cell =
        1.0 / (static_cast<double>(number_of_unit_cells.x * number_of_unit_cells.y * number_of_unit_cells.z));

    double loading = numberOfMolecules[comp.componentId];
    double excess_loading = numberOfMolecules[comp.componentId] - comp.amountOfExcessMolecules;
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

    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);
    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
        std::print(stream, "    molecules:        {: .6e} molecules\n", numberOfMolecules[comp.componentId]);
        std::print(stream, "    number density:   {: .6e} molec/A^3\n", numberDensities[comp.componentId]);
        std::print(stream, "    density:          {: .6e} kg/m^3\n",
                   densityConversionFactor * comp.totalMass * numberDensities[comp.componentId]);
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    molecules:        {: .6e} molecules\n", numberOfMolecules[comp.componentId]);
        std::print(stream, "    number density:   {: .6e} molec./{}^3\n", numberDensities[comp.componentId],
                   Units::displayedUnitOfLengthString);
        break;
    }
  }

  return stream.str();
}

std::string Loadings::printStatus(const Component &comp, const Loadings &average, const Loadings &error,
                                  std::optional<double> frameworkMass, std::optional<int3> numberOfUnitCells) const
{
  std::ostringstream stream;

  if (frameworkMass)
  {
    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);

    const double toMolePerKg = 1000.0 / frameworkMass.value();
    const double toMgPerKg = 1000.0 * comp.totalMass / frameworkMass.value();

    int3 number_of_unit_cells = numberOfUnitCells.value_or(int3{1, 1, 1});
    double to_molecules_per_unit_cell =
        1.0 / (static_cast<double>(number_of_unit_cells.x * number_of_unit_cells.y * number_of_unit_cells.z));

    double loading = numberOfMolecules[comp.componentId];
    double loading_avg = average.numberOfMolecules[comp.componentId];
    double loading_error = error.numberOfMolecules[comp.componentId];

    double excess_loading = numberOfMolecules[comp.componentId] - comp.amountOfExcessMolecules;
    double excess_loading_avg = average.numberOfMolecules[comp.componentId] - comp.amountOfExcessMolecules;
    double excess_loading_error = error.numberOfMolecules[comp.componentId];

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

    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);
    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
        std::print(stream, "    molecules:      {:.6e} molecules  ({:.6e} +/- {:.6e})\n",
                   numberOfMolecules[comp.componentId], average.numberOfMolecules[comp.componentId],
                   error.numberOfMolecules[comp.componentId]);
        std::print(stream, "    number density: {:.6e} molec./A^3 ({:.6e} +/- {:.6e})\n",
                   numberDensities[comp.componentId], average.numberDensities[comp.componentId],
                   error.numberDensities[comp.componentId]);
        std::print(stream, "    density:        {:.6e} kg/m^3     ({:.6e} +/- {:.6e})\n",
                   densityConversionFactor * comp.totalMass * numberDensities[comp.componentId],
                   densityConversionFactor * comp.totalMass * average.numberDensities[comp.componentId],
                   densityConversionFactor * comp.totalMass * error.numberDensities[comp.componentId]);
        break;
      case Units::System::ReducedUnits:
        std::print(stream, "    molecules:      {:.6e} molecules  ({:.6e} +/- {:.6e})\n",
                   numberOfMolecules[comp.componentId], average.numberOfMolecules[comp.componentId],
                   error.numberOfMolecules[comp.componentId]);
        std::print(stream, "    number density: {:.6e} molec./{}^3 ({:.6e} +/- {:.6e})\n",
                   numberDensities[comp.componentId], Units::displayedUnitOfLengthString,
                   average.numberDensities[comp.componentId], error.numberDensities[comp.componentId]);
        break;
    }
  }

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Loadings &l)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Loadings &l)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > l.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Loadings' at line {} in file {}\n", location.line(),
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
    throw std::runtime_error(std::format("Loadings: Error in binary restart\n"));
  }
#endif

  return archive;
}
