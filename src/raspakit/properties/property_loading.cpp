module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <optional>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <exception>
#include <source_location>
#include <complex>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module property_loading;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <iostream>;
import <sstream>;
import <fstream>;
import <optional>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;
import stringutils;
import units;
import loadings;
import component;

std::string PropertyLoading::writeAveragesStatistics(std::vector<Component> components, 
                                                     std::optional<double> frameworkMass) const
{
  std::ostringstream stream;

  if (frameworkMass.has_value())
  {
    const double toMolePerKg = 1000.0 / frameworkMass.value();

    std::pair<Loadings, Loadings> loadingAverage = averageLoading();

    std::print(stream, "Loadings\n");
    std::print(stream, "===============================================================================\n\n");


    for(size_t i = 0; i < components.size(); ++i)
    {
      if(components[i].type != Component::Type::Framework)
      {
        const double toMgPerG = 1000.0 * components[i].totalMass / frameworkMass.value();
        
        std::print(stream, "Component {} ({})\n", components[i].componentId, components[i].name);
        
        for (size_t j = 0; j < bookKeepingLoadings.size(); ++j)
        {
          Loadings blockAverage = averagedLoading(j);
          std::print(stream, "    Block[ {:2d}] {: .6e}\n", j, blockAverage.numberOfMolecules[i]);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                           loadingAverage.first.numberOfMolecules[i], 
                           loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [mol/kg-framework]\n",
                           toMolePerKg * loadingAverage.first.numberOfMolecules[i], 
                           toMolePerKg * loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Abs. loading average  {: .6e} +/- {: .6e} [mg/g-framework]\n",
                           toMgPerG * loadingAverage.first.numberOfMolecules[i], 
                           toMgPerG * loadingAverage.second.numberOfMolecules[i]);
        
        std::print(stream, "\n");

        for (size_t j = 0; j < bookKeepingLoadings.size(); ++j)
        {
          Loadings blockAverage = averagedLoading(j);
          std::print(stream, "    Block[ {:2d}] {: .6e}\n", j, 
                             blockAverage.numberOfMolecules[i] - components[i].amountOfExcessMolecules);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [molecules/cell]\n",
                           loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules, 
                           loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [mol/kg-framework]\n",
                   toMolePerKg * (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules), 
                   toMolePerKg * loadingAverage.second.numberOfMolecules[i]);
        std::print(stream, "    Excess loading average  {: .6e} +/- {: .6e} [mg/g-framework]\n",
                   toMgPerG * (loadingAverage.first.numberOfMolecules[i] - components[i].amountOfExcessMolecules), 
                   toMgPerG * loadingAverage.second.numberOfMolecules[i]);

        std::print(stream, "\n\n");
      }
    }
  }
  else
  {
    const double densityConversionFactor = 
      1.0 / (1000.0 * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::pair<Loadings, Loadings> loadingAverage = averageLoading();

    std::print(stream, "Densities\n");
    std::print(stream, "===============================================================================\n\n");

    for (size_t i = 0; i < components.size(); ++i)
    {
      std::print(stream, "Component {} ({})\n", components[i].componentId, components[i].name);

      for (size_t j = 0; j < bookKeepingLoadings.size(); ++j)
      {
        Loadings blockAverage = averagedLoading(j);
        std::print(stream, "    Block[ {:2d}] {}\n", j, blockAverage.numberOfMolecules[i]);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molecules]\n",
        loadingAverage.first.numberOfMolecules[i], loadingAverage.second.numberOfMolecules[i]);
      std::print(stream, "    Density average  {: .6e} +/- {: .6e} [molec/A^3]\n",
        loadingAverage.first.numberDensities[i], loadingAverage.second.numberDensities[i]);
      std::print(stream, "    Density average  {: .6e} +/- {: .6e} [kg/m^3]\n",
                         densityConversionFactor * components[i].totalMass * loadingAverage.first.numberDensities[i], 
                         densityConversionFactor * components[i].totalMass * loadingAverage.second.numberDensities[i]);
    }
  }

  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyLoading &l)
{
  archive << l.versionNumber;

  archive << l.numberOfBlocks;
  archive << l.numberOfComponents;
  archive << l.bookKeepingLoadings;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyLoading &l)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > l.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyLoading' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> l.numberOfBlocks;
  archive >> l.numberOfComponents;
  archive >> l.bookKeepingLoadings;
  return archive;
}
