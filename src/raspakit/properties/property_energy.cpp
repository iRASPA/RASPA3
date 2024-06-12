module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <random>
#include <sstream>
#include <fstream>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <print>
#endif

module property_energy;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <random>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <print>;
#endif


import archive;
import stringutils;
import framework;
import component;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import averages;
import units;


std::string PropertyEnergy::writeAveragesStatistics(bool externalField, std::vector<Framework>& frameworkComponents, std::vector<Component>& components) const
{
  std::ostringstream stream;

  std::pair<EnergyStatus, EnergyStatus> computedAverage = averageEnergy();

  std::print(stream, "Energy averages and statistics:\n");
  std::print(stream, "===============================================================================\n\n");

  if(externalField)
  {
    std::print(stream, "ExternalField-molecular energy:\n");
    std::print(stream, "-------------------------------------------------------------------------------\n\n");

    for (size_t k = 0; k < 1; k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if (k == l)
        {
          std::print(stream, "    ExternalField-molecule energy {}-{} [{}-{}]:\n", k, l,
            components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    ExternalField-molecule energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
            components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.externalFieldComponentEnergy(k, l).totalInter.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
          prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).totalInter.energy,
          prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).totalInter.energy);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n\n");

    std::print(stream, "ExternalField-molecule energy contributions per energy type:\n");
    std::print(stream, "-------------------------------------------------------------------------------\n\n");

    for (size_t k = 0; k < 1; k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if (k == l)
        {
          std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n", k, l,
            components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
            components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.externalFieldComponentEnergy(k, l).VanDerWaals.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
          prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).VanDerWaals.energy,
          prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).VanDerWaals.energy);
        std::print(stream, "\n");
      }
    }

    for (size_t k = 0; k < 1; k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if (k == l)
        {
          std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n", k, l,
            components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
            components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, 
                             prefactor * blockAverage.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
          prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy,
          prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
        std::print(stream, "\n");
      }
    }

    for (size_t k = 0; k < 1; k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if (k == l)
        {
          std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n", k, l,
            components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
            components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.externalFieldComponentEnergy(k, l).CoulombicReal.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
          prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).CoulombicReal.energy,
          prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).CoulombicReal.energy);
        std::print(stream, "\n");
      }
    }

    for (size_t k = 0; k < 1; k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if (k == l)
        {
          std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n", k, l,
            components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
            components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", 
                             i, prefactor * blockAverage.externalFieldComponentEnergy(k, l).CoulombicFourier.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
          prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).CoulombicFourier.energy,
          prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).CoulombicFourier.energy);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n\n");
  }


  std::print(stream, "Framework-molecule energy:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  for (size_t k = 0; k < frameworkComponents.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Framework-molecule energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Framework-molecule energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.frameworkComponentEnergy(k, l).totalInter.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.frameworkComponentEnergy(k, l).totalInter.energy,
        prefactor * computedAverage.second.frameworkComponentEnergy(k, l).totalInter.energy);
      std::print(stream, "\n");
    }
  }
  std::print(stream, "\n\n");

  std::print(stream, "Framework-molecule energy contributions per energy type:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  for (size_t k = 0; k < frameworkComponents.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.frameworkComponentEnergy(k, l).VanDerWaals.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.frameworkComponentEnergy(k, l).VanDerWaals.energy,
        prefactor * computedAverage.second.frameworkComponentEnergy(k, l).VanDerWaals.energy);
      std::print(stream, "\n");
    }
  }

  for (size_t k = 0; k < frameworkComponents.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, 
                           prefactor * blockAverage.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy,
        prefactor * computedAverage.second.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
      std::print(stream, "\n");
    }
  }

  for (size_t k = 0; k < frameworkComponents.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.frameworkComponentEnergy(k, l).CoulombicReal.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.frameworkComponentEnergy(k, l).CoulombicReal.energy,
        prefactor * computedAverage.second.frameworkComponentEnergy(k, l).CoulombicReal.energy);
      std::print(stream, "\n");
    }
  }

  for (size_t k = 0; k < frameworkComponents.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", 
                           i, prefactor * blockAverage.frameworkComponentEnergy(k, l).CoulombicFourier.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.frameworkComponentEnergy(k, l).CoulombicFourier.energy,
        prefactor * computedAverage.second.frameworkComponentEnergy(k, l).CoulombicFourier.energy);
      std::print(stream, "\n");
    }
  }
  std::print(stream, "\n\n");


  std::print(stream, "Inter-molecular energy:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  for (size_t k = 0; k < components.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Inter-molecular energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Inter-molecular energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.componentEnergy(k, l).totalInter.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.componentEnergy(k, l).totalInter.energy,
        prefactor * computedAverage.second.componentEnergy(k, l).totalInter.energy);
      std::print(stream, "\n");
    }
  }
  std::print(stream, "\n\n");

  std::print(stream, "Inter-molecular energy contributions per energy type:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  // Write Van der Waals intermolecular energy
  for (size_t k = 0; k < components.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.componentEnergy(k, l).VanDerWaals.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.componentEnergy(k, l).VanDerWaals.energy,
        prefactor * computedAverage.second.componentEnergy(k, l).VanDerWaals.energy);
      std::print(stream, "\n");
    }
  }

  for (size_t k = 0; k < components.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, 
                           prefactor * blockAverage.componentEnergy(k, l).VanDerWaalsTailCorrection.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.componentEnergy(k, l).VanDerWaalsTailCorrection.energy,
        prefactor * computedAverage.second.componentEnergy(k, l).VanDerWaalsTailCorrection.energy);
      std::print(stream, "\n");
    }
  }

  for (size_t k = 0; k < components.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.componentEnergy(k, l).CoulombicReal.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.componentEnergy(k, l).CoulombicReal.energy,
        prefactor * computedAverage.second.componentEnergy(k, l).CoulombicReal.energy);
      std::print(stream, "\n");
    }
  }

  for (size_t k = 0; k < components.size(); k++)
  {
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n", k, l,
          components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k,
          components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", 
                           i, prefactor * blockAverage.componentEnergy(k, l).CoulombicFourier.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first.componentEnergy(k, l).CoulombicFourier.energy,
        prefactor * computedAverage.second.componentEnergy(k, l).CoulombicFourier.energy);
      std::print(stream, "\n");
    }
  }
  std::print(stream, "\n");

  std::print(stream, "Total energy:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  double prefactor = Units::EnergyToKelvin;
  for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
  {
    EnergyStatus blockAverage = averagedEnergy(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.totalEnergy.energy);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Average  {: .6e} +/- {: .6e} [K]\n", prefactor * computedAverage.first.totalEnergy.energy, 
                                                               prefactor * computedAverage.second.totalEnergy.energy);

  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergy &e)
{
  archive << e.versionNumber;

  archive << e.numberOfBlocks;
  archive << e.numberOfExternalFields;
  archive << e.numberOfFrameworks;
  archive << e.numberOfComponents;
  archive << e.bookKeepingEnergyStatus;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergy &e)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > e.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyEnergy' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfBlocks;
  archive >> e.numberOfExternalFields;
  archive >> e.numberOfFrameworks;
  archive >> e.numberOfComponents;
  archive >> e.bookKeepingEnergyStatus;

  return archive;
}
