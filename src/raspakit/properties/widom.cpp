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
import <vector>;
import <map>;
import <array>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <cmath>;
import <numbers>;
import <optional>;
import <fstream>;
import <exception>;
import <source_location>;
import <complex>;
import <print>;
#endif

import archive;
import units;
import averages;
import stringutils;

PropertyWidom::PropertyWidom() {}

std::string PropertyWidom::writeAveragesStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                                   std::optional<double> imposedFugacity) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;

  std::print(stream, "    Widom insertion Rosenbluth weight  statistics:\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    double blockAverage = averagedRosenbluthWeight(blockIndex);
    std::print(stream, "        Block[ {:2d}] {: .6e}\n", blockIndex, blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::pair<double, double> averageRosenbluthWeightValue = averageRosenbluthWeight();
  std::print(stream, "    Average Rosenbluth weight:   {: .6e} +/- {: .6e} [-]\n", averageRosenbluthWeightValue.first,
             averageRosenbluthWeightValue.second);
  std::print(stream, "\n\n");

  switch (Units::unitSystem)
  {
    case Units::System::RASPA:
    {
      std::print(stream, "    Widom insertion chemical potential  statistics:\n");
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::pair<double, double> average_excess_widom_chemical_potential = averageExcessChemicalPotential(beta);
      std::pair<double, double> average_chemical_potential_tail_correction = averageChemicalPotentialTailCorrection();
      std::pair<double, double> average_ideal_gas_widom_chemical_potential = averageIdealGasChemicalPotential(beta);
      std::pair<double, double> average_total_widom_chemical_potential = averageTotalChemicalPotential(beta);
      std::pair<double, double> average_widom_fugacity = averageFugacity(beta);

      std::print(stream, "    Excess chemical potential:          {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_excess_widom_chemical_potential.first,
                 Units::EnergyToKelvin * average_excess_widom_chemical_potential.second);
      std::print(stream, "    Tail-correction chemical potential: {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_chemical_potential_tail_correction.first,
                 Units::EnergyToKelvin * average_chemical_potential_tail_correction.second);
      std::print(stream, "    Ideal chemical potential:           {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_ideal_gas_widom_chemical_potential.first,
                 Units::EnergyToKelvin * average_ideal_gas_widom_chemical_potential.second);
      std::print(stream, "    Total chemical potential:           {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * average_total_widom_chemical_potential.first,
                 Units::EnergyToKelvin * average_total_widom_chemical_potential.second);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:  {: .6e} [K]\n",
                   Units::EnergyToKelvin * imposedChemicalPotential.value());
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:   {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_excess_widom_chemical_potential.first,
                 Units::EnergyToKJPerMol * average_excess_widom_chemical_potential.second);
      std::print(stream, "    Tail-correction chemical potential: {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKJPerMol * average_chemical_potential_tail_correction.first,
                 Units::EnergyToKJPerMol * average_chemical_potential_tail_correction.second);
      std::print(stream, "    Ideal chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_ideal_gas_widom_chemical_potential.first,
                 Units::EnergyToKJPerMol * average_ideal_gas_widom_chemical_potential.second);
      std::print(stream, "    Total chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * average_total_widom_chemical_potential.first,
                 Units::EnergyToKJPerMol * average_total_widom_chemical_potential.second);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Imposed chemical potential:  {: .6e} [kJ/mol]\n",
                   Units::EnergyToKJPerMol * imposedChemicalPotential.value());
        std::print(stream, "    ---------------------------------------------------------------------------\n");
      }
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:            {: .6e} [Pa]\n",
                   Units::PressureConversionFactor * imposedFugacity.value());

        std::print(stream, "    Measured fugacity:           {: .6e} +/- {: .6e} [Pa]\n",
                   Units::PressureConversionFactor * average_widom_fugacity.first,
                   Units::PressureConversionFactor * average_widom_fugacity.second);
      }
    }
    break;
    case Units::System::ReducedUnits:
    {
      std::print(stream, "    Widom insertion chemical potential  statistics:\n");
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, beta * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::pair<double, double> average_excess_widom_chemical_potential = averageExcessChemicalPotential(beta);
      std::pair<double, double> average_chemical_potential_tail_correction = averageChemicalPotentialTailCorrection();
      std::pair<double, double> average_ideal_gas_widom_chemical_potential = averageIdealGasChemicalPotential(beta);
      std::pair<double, double> average_total_widom_chemical_potential = averageTotalChemicalPotential(beta);
      std::pair<double, double> average_widom_fugacity = averageFugacity(beta);

      std::print(stream, "    Beta * Excess chemical potential:          {: .6e} +/- {: .6e} [-]\n",
                 beta * average_excess_widom_chemical_potential.first,
                 beta * average_excess_widom_chemical_potential.second);
      std::print(stream, "    Beta * Tail-correction chemical potential: {: .6e} +/- {: .6e} [-]\n",
                 beta * average_chemical_potential_tail_correction.first,
                 beta * average_chemical_potential_tail_correction.second);
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
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w)
{
  uint64_t versionNumber;
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
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyWidom: Error in binary restart\n"));
  }
#endif

  return archive;
}
