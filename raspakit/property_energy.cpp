module;

module property_energy;

import <iostream>;
import <random>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#else
  import print;
#endif


import archive;
import stringutils;
import component;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import averages;
import units;


std::string PropertyEnergy::writeAveragesStatistics(std::vector<Component>& components) const
{
  std::ostringstream stream;

  std::pair<EnergyStatus, EnergyStatus> computedAverage = averageEnergy();

  std::print(stream, "Energy averages and statistics:\n");
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "Inter-molecular energy:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n\n");

  // Write total intermolecular energy
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
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage(k, l).totalInter.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first(k, l).totalInter.energy,
        prefactor * computedAverage.second(k, l).totalInter.energy);
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
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage(k, l).VanDerWaals.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first(k, l).VanDerWaals.energy,
        prefactor * computedAverage.second(k, l).VanDerWaals.energy);
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
                           prefactor * blockAverage(k, l).VanDerWaalsTailCorrection.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first(k, l).VanDerWaalsTailCorrection.energy,
        prefactor * computedAverage.second(k, l).VanDerWaalsTailCorrection.energy);
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
        std::print(stream, "        Block[ {:2d}] {: .6e}\n", i, prefactor * blockAverage(k, l).CoulombicReal.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first(k, l).CoulombicReal.energy,
        prefactor * computedAverage.second(k, l).CoulombicReal.energy);
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
                           i, prefactor * blockAverage(k, l).CoulombicFourier.energy);
      }
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Average  {: .6e} +/- {: .6e} [K]\n",
        prefactor * computedAverage.first(k, l).CoulombicFourier.energy,
        prefactor * computedAverage.second(k, l).CoulombicFourier.energy);
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
  archive >> e.numberOfComponents;
  archive >> e.bookKeepingEnergyStatus;

  return archive;
}
