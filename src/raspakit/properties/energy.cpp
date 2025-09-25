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
#include <print>
#include <random>
#include <source_location>
#include <sstream>
#include <vector>
#endif

module property_energy;

#ifndef USE_LEGACY_HEADERS
import std;
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

std::string PropertyEnergy::writeAveragesStatistics(bool externalField, std::optional<Framework> &framework,
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

    for (std::size_t k = 0; k < 1; k++)
    {
      for (std::size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    ExternalField-molecule energy{} {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[k].name, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).totalInter.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).totalInter.energy,
                   prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).totalInter.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n\n");

    std::print(stream, "ExternalField-molecule energy contributions per energy type:\n");
    std::print(stream, "-------------------------------------------------------------------------------\n\n");

    for (std::size_t k = 0; k < 1; k++)
    {
      for (std::size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Van der Waals energy{} {}-{} [{}-{}]:\n", Units::displayedUnitOfEnergyConversionString,
                   k, l, components[k].name, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).VanDerWaals.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).VanDerWaals.energy,
                   prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).VanDerWaals.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }

    for (std::size_t k = 0; k < 1; k++)
    {
      for (std::size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    VDW Tail-Correction energy{} {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[k].name, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(
            stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
            prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy,
            prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).VanDerWaalsTailCorrection.energy,
            Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }

    for (std::size_t k = 0; k < 1; k++)
    {
      for (std::size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Coulomb Real energy{} {}-{} [{}-{}]:\n", Units::displayedUnitOfEnergyConversionString,
                   k, l, components[k].name, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).CoulombicReal.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).CoulombicReal.energy,
                   prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).CoulombicReal.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }

    for (std::size_t k = 0; k < 1; k++)
    {
      for (std::size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Coulomb Fourier energy{} {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[k].name, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.externalFieldComponentEnergy(k, l).CoulombicFourier.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.externalFieldComponentEnergy(k, l).CoulombicFourier.energy,
                   prefactor * computedAverage.second.externalFieldComponentEnergy(k, l).CoulombicFourier.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n\n");
  }

  if (framework.has_value())
  {
    std::print(stream, "Framework-molecule energy:\n");
    std::print(stream, "-------------------------------------------------------------------------------\n\n");

    for (std::size_t k = 0; k < framework->numberOfComponents; k++)
    {
      for (std::size_t l = 0; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Framework-molecule energy{} {}-{} [{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.frameworkComponentEnergy(k, l).totalInter.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.frameworkComponentEnergy(k, l).totalInter.energy,
                   prefactor * computedAverage.second.frameworkComponentEnergy(k, l).totalInter.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n\n");

    std::print(stream, "Framework-molecule energy contributions per energy type:\n");
    std::print(stream, "-------------------------------------------------------------------------------\n\n");

    for (std::size_t k = 0; k < framework->numberOfComponents; k++)
    {
      for (std::size_t l = 0; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Van der Waals energy{} {}-{} [{}]:\n", Units::displayedUnitOfEnergyConversionString, k,
                   l, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.frameworkComponentEnergy(k, l).VanDerWaals.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.frameworkComponentEnergy(k, l).VanDerWaals.energy,
                   prefactor * computedAverage.second.frameworkComponentEnergy(k, l).VanDerWaals.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }

    for (std::size_t k = 0; k < framework->numberOfComponents; k++)
    {
      for (std::size_t l = 0; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    VDW Tail-Correction energy{} {}-{} [{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy,
                   prefactor * computedAverage.second.frameworkComponentEnergy(k, l).VanDerWaalsTailCorrection.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }

    for (std::size_t k = 0; k < framework->numberOfComponents; k++)
    {
      for (std::size_t l = 0; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Coulomb Real energy{} {}-{} [{}]:\n", Units::displayedUnitOfEnergyConversionString, k,
                   l, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.frameworkComponentEnergy(k, l).CoulombicReal.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.frameworkComponentEnergy(k, l).CoulombicReal.energy,
                   prefactor * computedAverage.second.frameworkComponentEnergy(k, l).CoulombicReal.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }

    for (std::size_t k = 0; k < framework->numberOfComponents; k++)
    {
      for (std::size_t l = 0; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::print(stream, "    Coulomb Fourier energy{} {}-{} [{}]:\n", Units::displayedUnitOfEnergyConversionString,
                   k, l, components[l].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
        {
          EnergyStatus blockAverage = averagedEnergy(i);
          std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                     prefactor * blockAverage.frameworkComponentEnergy(k, l).CoulombicFourier.energy);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                   prefactor * computedAverage.first.frameworkComponentEnergy(k, l).CoulombicFourier.energy,
                   prefactor * computedAverage.second.frameworkComponentEnergy(k, l).CoulombicFourier.energy,
                   Units::displayedUnitOfEnergyString);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n\n");
  }

  std::print(stream, "Inter-molecular energy:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  for (std::size_t k = 0; k < components.size(); k++)
  {
    for (std::size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Inter-molecular energy{} {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Inter-molecular energy{} {}-{} + {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, l, k, components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).totalInter.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.componentEnergy(k, l).totalInter.energy,
                 prefactor * computedAverage.second.componentEnergy(k, l).totalInter.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }
  }
  std::print(stream, "\n\n");

  std::print(stream, "Inter-molecular energy contributions per energy type:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  // Write Van der Waals intermolecular energy
  for (std::size_t k = 0; k < components.size(); k++)
  {
    for (std::size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Van der Waals energy{} {}-{} [{}-{}]:\n", Units::displayedUnitOfEnergyConversionString,
                   k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Van der Waals energy{} {}-{} + {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, l, k, components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).VanDerWaals.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.componentEnergy(k, l).VanDerWaals.energy,
                 prefactor * computedAverage.second.componentEnergy(k, l).VanDerWaals.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }
  }

  for (std::size_t k = 0; k < components.size(); k++)
  {
    for (std::size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    VDW Tail-Correction energy{} {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    VDW Tail-Correction energy{} {}-{} + {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, l, k, components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).VanDerWaalsTailCorrection.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.componentEnergy(k, l).VanDerWaalsTailCorrection.energy,
                 prefactor * computedAverage.second.componentEnergy(k, l).VanDerWaalsTailCorrection.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }
  }

  for (std::size_t k = 0; k < components.size(); k++)
  {
    for (std::size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Coulomb Real energy{} {}-{} [{}-{}]:\n", Units::displayedUnitOfEnergyConversionString,
                   k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Real energy{} {}-{} + {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, l, k, components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).CoulombicReal.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.componentEnergy(k, l).CoulombicReal.energy,
                 prefactor * computedAverage.second.componentEnergy(k, l).CoulombicReal.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }
  }

  for (std::size_t k = 0; k < components.size(); k++)
  {
    for (std::size_t l = k; l < components.size(); l++)
    {
      double prefactor = Units::EnergyToKelvin;
      if (k == l)
      {
        std::print(stream, "    Coulomb Fourier energy{} {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, components[k].name, components[l].name);
      }
      else
      {
        prefactor *= 2.0;
        std::print(stream, "    Coulomb Fourier energy{} {}-{} + {}-{} [{}-{}]:\n",
                   Units::displayedUnitOfEnergyConversionString, k, l, l, k, components[k].name, components[l].name);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.componentEnergy(k, l).CoulombicFourier.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.componentEnergy(k, l).CoulombicFourier.energy,
                 prefactor * computedAverage.second.componentEnergy(k, l).CoulombicFourier.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }
  }

  std::print(stream, "Intra-molecular energies per energy type:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  for (std::size_t k = 0; k < components.size(); k++)
  {
    if (!components[k].intraMolecularPotentials.bonds.empty())
    {
      std::print(stream, "    Bond energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bond);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bond,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bond, Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.ureyBradleys.empty())
    {
      std::print(stream, "    Urey-Bradley energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].ureyBradley);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].ureyBradley,
                 prefactor * computedAverage.second.intraComponentEnergies[k].ureyBradley,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.bends.empty())
    {
      std::print(stream, "    Bend energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bend);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bend,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bend, Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.inversionBends.empty())
    {
      std::print(stream, "    Inversion-bend energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].inversionBend);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].inversionBend,
                 prefactor * computedAverage.second.intraComponentEnergies[k].inversionBend,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.outOfPlaneBends.empty())
    {
      std::print(stream, "    Out-of-plane bend energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].outOfPlaneBend);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].outOfPlaneBend,
                 prefactor * computedAverage.second.intraComponentEnergies[k].outOfPlaneBend,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.torsions.empty())
    {
      std::print(stream, "    Torsion energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].torsion);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].torsion,
                 prefactor * computedAverage.second.intraComponentEnergies[k].torsion,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.improperTorsions.empty())
    {
      std::print(stream, "    Improper torsion energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].improperTorsion);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].improperTorsion,
                 prefactor * computedAverage.second.intraComponentEnergies[k].improperTorsion,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.bondBonds.empty())
    {
      std::print(stream, "    Bond-bond energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bondBond);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bondBond,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bondBond,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.bondBends.empty())
    {
      std::print(stream, "    Bond-bend energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bondBend);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bondBend,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bondBend,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.bondTorsions.empty())
    {
      std::print(stream, "    Bond-torsion energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bondTorsion);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bondTorsion,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bondTorsion,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.bendBends.empty())
    {
      std::print(stream, "    Bend-bend energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bendBend);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bendBend,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bendBend,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.bendTorsions.empty())
    {
      std::print(stream, "    Bend-Torsion energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].bendTorsion);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].bendTorsion,
                 prefactor * computedAverage.second.intraComponentEnergies[k].bendTorsion,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.vanDerWaals.empty())
    {
      std::print(stream, "    Intra Van Der Waals energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].vanDerWaals);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].vanDerWaals,
                 prefactor * computedAverage.second.intraComponentEnergies[k].vanDerWaals,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }

    if (!components[k].intraMolecularPotentials.coulombs.empty())
    {
      std::print(stream, "    Intra Van Der Waals energy{} {} [{}]\n", Units::displayedUnitOfEnergyConversionString, k,
                 components[k].name);
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      double prefactor = Units::EnergyToKelvin;
      for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
      {
        EnergyStatus blockAverage = averagedEnergy(i);
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i,
                   prefactor * blockAverage.intraComponentEnergies[k].coulomb);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
                 prefactor * computedAverage.first.intraComponentEnergies[k].coulomb,
                 prefactor * computedAverage.second.intraComponentEnergies[k].coulomb,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "\n");
    }
  }
  std::print(stream, "\n");

  double prefactor = Units::EnergyToKelvin;

  std::print(stream, "Kinetic Energies:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");
  std::print(stream, "    Translational Kinetic energy{}\n", Units::displayedUnitOfEnergyConversionString);
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
  {
    EnergyStatus blockAverage = averagedEnergy(i);
    std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.translationalKineticEnergy);
  }
  std::print(stream, "        -----------------------------------------------------------------------\n");
  std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
             prefactor * computedAverage.first.translationalKineticEnergy,
             prefactor * computedAverage.second.translationalKineticEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  std::print(stream, "    Rotational Kinetic energy{}\n", Units::displayedUnitOfEnergyConversionString);
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
  {
    EnergyStatus blockAverage = averagedEnergy(i);
    std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.rotationalKineticEnergy);
  }
  std::print(stream, "        -----------------------------------------------------------------------\n");
  std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n",
             prefactor * computedAverage.first.rotationalKineticEnergy,
             prefactor * computedAverage.second.rotationalKineticEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  std::print(stream, "    Nose Hoover energy{}\n", Units::displayedUnitOfEnergyConversionString);
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
  {
    EnergyStatus blockAverage = averagedEnergy(i);
    std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.noseHooverEnergy);
  }
  std::print(stream, "        -----------------------------------------------------------------------\n");
  std::print(stream, "        Average  {: .6e} +/- {: .6e} [{}]\n", prefactor * computedAverage.first.noseHooverEnergy,
             prefactor * computedAverage.second.noseHooverEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  std::print(stream, "Polarization energy{}\n", Units::displayedUnitOfEnergyConversionString);
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
  {
    EnergyStatus blockAverage = averagedEnergy(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.polarizationEnergy.energy);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Average  {: .6e} +/- {: .6e} [{}]\n",
             prefactor * computedAverage.first.polarizationEnergy.energy,
             prefactor * computedAverage.second.polarizationEnergy.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  std::print(stream, "Total energy{}\n", Units::displayedUnitOfEnergyConversionString);
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (std::size_t i = 0; i < bookKeepingEnergyStatus.size(); ++i)
  {
    EnergyStatus blockAverage = averagedEnergy(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage.totalEnergy.energy);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Average  {: .6e} +/- {: .6e} [{}]\n", prefactor * computedAverage.first.totalEnergy.energy,
             prefactor * computedAverage.second.totalEnergy.energy, Units::displayedUnitOfEnergyString);

  std::print(stream, "\n");

  return stream.str();
}

nlohmann::json PropertyEnergy::jsonAveragesStatistics(bool externalField, std::optional<Framework> &framework,
                                                      std::vector<Component> &components) const
{
  nlohmann::json status;

  std::pair<EnergyStatus, EnergyStatus> computedAverage = averageEnergy();
  std::vector<EnergyStatus> blockEnergies = blockEnergy();

  for (std::size_t k = 0; k < components.size(); k++)
  {
    if (externalField)
    {
      std::size_t l = 0;
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

    if (framework.has_value())
    {
      for (std::size_t l = 0; l < framework->numberOfComponents; l++)
      {
        double prefactor = Units::EnergyToKelvin;
        std::string pair = std::format("{}", components[k].name);
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

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergy &e)
{
  std::uint64_t versionNumber;
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

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyEnergy: Error in binary restart\n"));
  }
#endif

  return archive;
}
