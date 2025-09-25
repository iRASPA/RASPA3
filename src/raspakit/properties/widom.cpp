module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
#endif

module property_widom;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import int3;
import units;
import averages;
import stringutils;

PropertyWidom::PropertyWidom() {}

std::string PropertyWidom::writeAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                                   std::optional<double> frameworkMass,
                                                                   std::optional<int3> number_of_unit_cells) const
{
  std::ostringstream stream;

  std::print(stream, "    Widom insertion Rosenbluth weight statistics:\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    double blockAverage = averagedRosenbluthWeight(blockIndex);
    std::print(stream, "        Block[ {:2d}] {: .6e}\n", blockIndex, blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::pair<double, double> averageRosenbluthWeightValue = averageRosenbluthWeight();
  std::print(stream, "    Average Rosenbluth weight:   {: .6e} +/- {: .6e} [-]\n", averageRosenbluthWeightValue.first,
             averageRosenbluthWeightValue.second);
  std::print(stream, "\n\n");

  if (frameworkMass.has_value())
  {
    double frameworkDensity =
        1e-3 * frameworkMass.value() /
        (volume * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit * Units::AvogadroConstant);
    double conversion_factor_mol_per_kg = 1.0 / (Units::MolarGasConstant * temperature * frameworkDensity);

    std::print(stream, "    Henry coefficient based on Rosenbluth weight:\n");
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      double blockAverage = conversion_factor_mol_per_kg * averagedRosenbluthWeight(blockIndex);
      std::print(stream, "        Block[ {:2d}] {: .6e}\n", blockIndex, blockAverage);
    }
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    std::print(stream, "    Average Henry coefficient:   {: .6e} +/- {: .6e} [mol/kg/Pa]\n",
               averageRosenbluthWeightValue.first * conversion_factor_mol_per_kg,
               averageRosenbluthWeightValue.second * conversion_factor_mol_per_kg);
    if (number_of_unit_cells.has_value())
    {
      double conversion_factor_molecules_per_uc =
          Units::AvogadroConstant * volume * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit /
          (Units::MolarGasConstant * temperature *
           static_cast<double>(number_of_unit_cells->x * number_of_unit_cells->y * number_of_unit_cells->z));

      std::print(stream, "    Average Henry coefficient:   {: .6e} +/- {: .6e} [molec./uc/Pa]\n",
                 averageRosenbluthWeightValue.first * conversion_factor_molecules_per_uc,
                 averageRosenbluthWeightValue.second * conversion_factor_molecules_per_uc);
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

  switch (Units::unitSystem)
  {
    case Units::System::RASPA:
    {
      std::print(stream, "    Widom insertion chemical potential  statistics:\n");
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::pair<double, double> average_excess_widom_chemical_potential = averageExcessChemicalPotential(beta);
      std::pair<double, double> average_ideal_gas_widom_chemical_potential = averageIdealGasChemicalPotential(beta);
      std::pair<double, double> average_total_widom_chemical_potential = averageTotalChemicalPotential(beta);
      std::pair<double, double> average_widom_fugacity = averageFugacity(beta);

      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_excess_widom_chemical_potential.first,
                 Units::EnergyToKelvin * average_excess_widom_chemical_potential.second);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_ideal_gas_widom_chemical_potential.first,
                 Units::EnergyToKelvin * average_ideal_gas_widom_chemical_potential.second);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_total_widom_chemical_potential.first,
                 Units::EnergyToKelvin * average_total_widom_chemical_potential.second);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:         {: .6e} [K]\n",
                   Units::EnergyToKelvin * imposedChemicalPotential.value());
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_excess_widom_chemical_potential.first,
                 Units::EnergyToKJPerMol * average_excess_widom_chemical_potential.second);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_ideal_gas_widom_chemical_potential.first,
                 Units::EnergyToKJPerMol * average_ideal_gas_widom_chemical_potential.second);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_total_widom_chemical_potential.first,
                 Units::EnergyToKJPerMol * average_total_widom_chemical_potential.second);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:         {: .6e} [kJ/mol]\n",
                   Units::EnergyToKJPerMol * imposedChemicalPotential.value());
        std::print(stream, "    ---------------------------------------------------------------------------\n");
      }
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:                   {: .6e} [Pa]\n",
                   Units::PressureConversionFactor * imposedFugacity.value());

        std::print(stream, "    Measured fugacity:                  {: .6e} +/- {: .6e} [Pa]\n",
                   Units::PressureConversionFactor * average_widom_fugacity.first,
                   Units::PressureConversionFactor * average_widom_fugacity.second);
      }
    }
    break;
    case Units::System::ReducedUnits:
    {
      std::print(stream, "    Widom insertion chemical potential  statistics:\n");
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, beta * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::pair<double, double> average_excess_widom_chemical_potential = averageExcessChemicalPotential(beta);
      std::pair<double, double> average_ideal_gas_widom_chemical_potential = averageIdealGasChemicalPotential(beta);
      std::pair<double, double> average_total_widom_chemical_potential = averageTotalChemicalPotential(beta);
      std::pair<double, double> average_widom_fugacity = averageFugacity(beta);

      std::print(stream, "    Beta * Excess chemical potential:          {: .6e} +/- {: .6e} [-]\n",
                 beta * average_excess_widom_chemical_potential.first,
                 beta * average_excess_widom_chemical_potential.second);
      std::print(stream, "    Beta * Ideal chemical potential:           {: .6e} +/- {: .6e} [-]\n",
                 beta * average_ideal_gas_widom_chemical_potential.first,
                 beta * average_ideal_gas_widom_chemical_potential.second);
      std::print(stream, "    Beta * Total chemical potential:           {: .6e} +/- {: .6e} [-]\n",
                 beta * average_total_widom_chemical_potential.first,
                 beta * average_total_widom_chemical_potential.second);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Beta * Imposed chemical potential:  {: .6e} [-]\n",
                   beta * imposedChemicalPotential.value());
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [{}]\n",
                 average_excess_widom_chemical_potential.first, average_excess_widom_chemical_potential.second,
                 Units::unitOfEnergyString);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [{}]\n",
                 average_ideal_gas_widom_chemical_potential.first, average_ideal_gas_widom_chemical_potential.second,
                 Units::unitOfEnergyString);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [{}]\n",
                 average_total_widom_chemical_potential.first, average_total_widom_chemical_potential.second,
                 Units::unitOfEnergyString);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:  {: .6e} [{}]\n", imposedChemicalPotential.value(),
                   Units::unitOfEnergyString);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:            {: .6e} [{}]\n",
                   Units::PressureConversionFactor * imposedFugacity.value(), Units::unitOfPressureString);

        std::print(stream, "    Measured fugacity:           {: .6e} +/- {: .6e} [{}]\n",
                   Units::PressureConversionFactor * average_widom_fugacity.first,
                   Units::PressureConversionFactor * average_widom_fugacity.second, Units::unitOfPressureString);
      }
      break;
    }
  }

  std::print(stream, "\n\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w)
{
  archive << w.versionNumber;

  archive << w.numberOfBlocks;
  archive << w.bookKeepingWidom;
  archive << w.bookKeepingDensity;

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
  archive >> w.bookKeepingWidom;
  archive >> w.bookKeepingDensity;

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
