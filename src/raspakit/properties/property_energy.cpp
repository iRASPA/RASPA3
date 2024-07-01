module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <print>
#include <random>
#include <source_location>
#include <sstream>
#include <vector>
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
import json;

std::string PropertyEnergy::writeAveragesStatistics(bool externalField, std::vector<Framework> &frameworkComponents,
                                                    std::vector<Component> &components) const
{
  std::ostringstream stream;

  std::pair<EnergyStatus, EnergyStatus> computedAverage = averageEnergy();

  std::print(stream, "Energy averages and statistics:\n");
  std::print(stream, "===============================================================================\n\n");

  if (externalField)
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
          std::print(stream, "    ExternalField-molecule energy {}-{} [{}-{}]:\n", k, l, components[k].name,
                     components[l].name);
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
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).totalInter.energy);
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
          std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                     components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).VanDerWaals.energy);
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
          std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n", k, l, components[k].name,
                     components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                     components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(
            stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
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
          std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                     components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).CoulombicReal.energy);
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
          std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n", k, l, components[k].name,
                     components[l].name);
        }
        else
        {
          prefactor *= 2.0;
          std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                     components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).CoulombicFourier.energy);
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
        std::print(stream, "    Framework-molecule energy {}-{} [{}-{}]:\n", k, l, components[k].name,
                   components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Framework-molecule energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.frameworkComponentEnergy(k, l).totalInter.energy);
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
        std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.frameworkComponentEnergy(k, l).VanDerWaals.energy);
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
        std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n", k, l, components[k].name,
                   components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
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
        std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.frameworkComponentEnergy(k, l).CoulombicReal.energy);
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
        std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.frameworkComponentEnergy(k, l).CoulombicFourier.energy);
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
        std::print(stream, "    Inter-molecular energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Inter-molecular energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).totalInter.energy);
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
        std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).VanDerWaals.energy);
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
        std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n", k, l, components[k].name,
                   components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
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
        std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).CoulombicReal.energy);
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
        std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n", k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n", k, l, l, k, components[k].name,
                   components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).CoulombicFourier.energy);
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

nlohmann::json PropertyEnergy::jsonAveragesStatistics(bool externalField, std::vector<Framework> &frameworkComponents,
                                                      std::vector<Component> &components) const
{
  nlohmann::json status;

  std::pair<EnergyStatus, EnergyStatus> computedAverage = averageEnergy();
  std::vector<EnergyStatus> blockEnergies = blockEnergy();

  for (size_t k = 0; k < components.size(); k++)
  {
    if (externalField)
    {
      size_t l = 0;
      std::string pair = std::format("{}", components[k].name);
      double prefactor = Units::EnergyToKelvin;

      std::vector<double> tmp(numberOfBlocks);
      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.externalFieldComponentEnergy(l, k).totalInter.energy; });
      status["ExternalField-Molecule"][pair]["total"]["block"] = tmp;
      status["ExternalField-Molecule"][pair]["total"]["mean"] =
          prefactor * computedAverage.first.externalFieldComponentEnergy(l, k).totalInter.energy;
      status["ExternalField-Molecule"][pair]["total"]["confidence"] =
          prefactor * computedAverage.second.externalFieldComponentEnergy(l, k).totalInter.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.externalFieldComponentEnergy(l, k).VanDerWaals.energy; });
      status["ExternalField-Molecule"][pair]["vanDerWaals"]["block"] = tmp;
      status["ExternalField-Molecule"][pair]["vanDerWaals"]["mean"] =
          prefactor * computedAverage.first.externalFieldComponentEnergy(l, k).VanDerWaals.energy;
      status["ExternalField-Molecule"][pair]["vanDerWaals"]["confidence"] =
          prefactor * computedAverage.second.externalFieldComponentEnergy(l, k).VanDerWaals.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.externalFieldComponentEnergy(l, k).VanDerWaalsTailCorrection.energy; });
      status["ExternalField-Molecule"][pair]["tailCorrection"]["block"] = tmp;
      status["ExternalField-Molecule"][pair]["tailCorrection"]["mean"] =
          prefactor * computedAverage.first.externalFieldComponentEnergy(l, k).VanDerWaalsTailCorrection.energy;
      status["ExternalField-Molecule"][pair]["tailCorrection"]["confidence"] =
          prefactor * computedAverage.second.externalFieldComponentEnergy(l, k).VanDerWaalsTailCorrection.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.externalFieldComponentEnergy(l, k).CoulombicReal.energy; });
      status["ExternalField-Molecule"][pair]["coulombReal"]["block"] = tmp;
      status["ExternalField-Molecule"][pair]["coulombReal"]["mean"] =
          prefactor * computedAverage.first.externalFieldComponentEnergy(l, k).CoulombicReal.energy;
      status["ExternalField-Molecule"][pair]["coulombReal"]["confidence"] =
          prefactor * computedAverage.second.externalFieldComponentEnergy(l, k).CoulombicReal.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.externalFieldComponentEnergy(l, k).CoulombicFourier.energy; });
      status["ExternalField-Molecule"][pair]["coulombFourier"]["block"] = tmp;
      status["ExternalField-Molecule"][pair]["coulombFourier"]["mean"] =
          prefactor * computedAverage.first.externalFieldComponentEnergy(l, k).CoulombicFourier.energy;
      status["ExternalField-Molecule"][pair]["coulombFourier"]["confidence"] =
          prefactor * computedAverage.second.externalFieldComponentEnergy(l, k).CoulombicFourier.energy;
    }

    for (size_t l = 0; l < frameworkComponents.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      std::string pair = std::format("{}-{}", components[k].name, frameworkComponents[l].name);
      std::vector<double> tmp(numberOfBlocks);

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.frameworkComponentEnergy(l, k).totalInter.energy; });
      status["Framework-Molecule"][pair]["total"]["block"] = tmp;
      status["Framework-Molecule"][pair]["total"]["mean"] =
          prefactor * computedAverage.first.frameworkComponentEnergy(k, l).totalInter.energy;
      status["Framework-Molecule"][pair]["total"]["confidence"] =
          prefactor * computedAverage.second.frameworkComponentEnergy(k, l).totalInter.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.frameworkComponentEnergy(l, k).VanDerWaals.energy; });
      status["Framework-Molecule"][pair]["vanDerWaals"]["block"] = tmp;
      status["Framework-Molecule"][pair]["vanDerWaals"]["mean"] =
          prefactor * computedAverage.first.frameworkComponentEnergy(k, l).VanDerWaals.energy;
      status["Framework-Molecule"][pair]["vanDerWaals"]["confidence"] =
          prefactor * computedAverage.second.frameworkComponentEnergy(k, l).VanDerWaals.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.frameworkComponentEnergy(l, k).VanDerWaalsTailCorrection.energy; });
      status["Framework-Molecule"][pair]["tailCorrection"]["block"] = tmp;
      status["Framework-Molecule"][pair]["tailCorrection"]["mean"] =
          prefactor * computedAverage.first.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy;
      status["Framework-Molecule"][pair]["tailCorrection"]["confidence"] =
          prefactor * computedAverage.second.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.frameworkComponentEnergy(l, k).CoulombicReal.energy; });
      status["Framework-Molecule"][pair]["coulombReal"]["block"] = tmp;
      status["Framework-Molecule"][pair]["coulombReal"]["mean"] =
          prefactor * computedAverage.first.frameworkComponentEnergy(k, l).CoulombicReal.energy;
      status["Framework-Molecule"][pair]["coulombReal"]["confidence"] =
          prefactor * computedAverage.second.frameworkComponentEnergy(k, l).CoulombicReal.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.frameworkComponentEnergy(l, k).CoulombicFourier.energy; });
      status["Framework-Molecule"][pair]["coulombFourier"]["block"] = tmp;
      status["Framework-Molecule"][pair]["coulombFourier"]["mean"] =
          prefactor * computedAverage.first.frameworkComponentEnergy(k, l).CoulombicFourier.energy;
      status["Framework-Molecule"][pair]["coulombFourier"]["confidence"] =
          prefactor * computedAverage.second.frameworkComponentEnergy(k, l).CoulombicFourier.energy;
    }
    for (size_t l = k; l < components.size(); l++)
    {
      double prefactor = (k == l) ? 2 * Units::EnergyToKelvin : Units::EnergyToKelvin;
      std::string pair = std::format("{}-{}", components[k].name, components[l].name);

      std::vector<double> tmp(numberOfBlocks);
      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.componentEnergy(l, k).totalInter.energy; });
      status["Molecule-Molecule"][pair]["total"]["block"] = tmp;
      status["Molecule-Molecule"][pair]["total"]["mean"] =
          prefactor * computedAverage.first.componentEnergy(k, l).totalInter.energy;
      status["Molecule-Molecule"][pair]["total"]["confidence"] =
          prefactor * computedAverage.second.componentEnergy(k, l).totalInter.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.componentEnergy(l, k).VanDerWaals.energy; });
      status["Molecule-Molecule"][pair]["vanDerWaals"]["block"] = tmp;
      status["Molecule-Molecule"][pair]["vanDerWaals"]["mean"] =
          prefactor * computedAverage.first.componentEnergy(k, l).VanDerWaals.energy;
      status["Molecule-Molecule"][pair]["vanDerWaals"]["confidence"] =
          prefactor * computedAverage.second.componentEnergy(k, l).VanDerWaals.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.componentEnergy(l, k).VanDerWaalsTailCorrection.energy; });
      status["Molecule-Molecule"][pair]["tailCorrection"]["block"] = tmp;
      status["Molecule-Molecule"][pair]["tailCorrection"]["mean"] =
          prefactor * computedAverage.first.componentEnergy(k, l).VanDerWaalsTailCorrection.energy;
      status["Molecule-Molecule"][pair]["tailCorrection"]["confidence"] =
          prefactor * computedAverage.second.componentEnergy(k, l).VanDerWaalsTailCorrection.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.componentEnergy(l, k).CoulombicReal.energy; });
      status["Molecule-Molecule"][pair]["coulombReal"]["block"] = tmp;
      status["Molecule-Molecule"][pair]["coulombReal"]["mean"] =
          prefactor * computedAverage.first.componentEnergy(k, l).CoulombicReal.energy;
      status["Molecule-Molecule"][pair]["coulombReal"]["confidence"] =
          prefactor * computedAverage.second.componentEnergy(k, l).CoulombicReal.energy;

      std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(), [prefactor, k, l](EnergyStatus &block)
                     { return prefactor * block.componentEnergy(l, k).CoulombicFourier.energy; });
      status["Molecule-Molecule"][pair]["coulombFourier"]["block"] = tmp;
      status["Molecule-Molecule"][pair]["coulombFourier"]["mean"] =
          prefactor * computedAverage.first.componentEnergy(k, l).CoulombicFourier.energy;
      status["Molecule-Molecule"][pair]["coulombFourier"]["confidence"] =
          prefactor * computedAverage.second.componentEnergy(k, l).CoulombicFourier.energy;
    }
  }

  double prefactor = Units::EnergyToKelvin;
  std::vector<double> tmp(numberOfBlocks);
  std::transform(blockEnergies.begin(), blockEnergies.end(), tmp.begin(),
                 [prefactor](EnergyStatus &block) { return prefactor * block.totalEnergy.energy; });
  status["totalEnergy"]["block"] = tmp;
  status["totalEnergy"]["mean"] = prefactor * computedAverage.first.totalEnergy.energy;
  status["totalEnergy"]["confidence"] = prefactor * computedAverage.second.totalEnergy.energy;

  return status;
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
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
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
