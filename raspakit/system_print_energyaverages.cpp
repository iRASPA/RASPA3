module;

module system;

import component;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import averages;
import print;
import units;
import property_energy;

import <iostream>;
import <random>;
import <sstream>;

std::string System::writeEnergyAveragesStatistics() const
{
    std::ostringstream stream;

    std::pair<EnergyStatus, EnergyStatus> computedAverage = averageEnergies.averageEnergy();

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
        if(k==l)
        { 
          std::print(stream, "    Inter-molecular energy {}-{} [{}-{}]:\n",k,l,
                  components[k].name, components[l].name);
        }
        else 
        {
          prefactor *= 2.0;
          std::print(stream, "    Inter-molecular energy {}-{} + {}-{} [{}-{}]:\n",k,l,l,k,
                  components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < averageEnergies.bookKeepingEnergyStatus.size(); ++i)
        {
            EnergyStatus blockAverage = averageEnergies.averagedEnergy(i);
            std::print(stream, "        Block[ {:2d}] {}\n", i, prefactor * blockAverage(k,l).totalInter);
        }
            std::print(stream, "        -----------------------------------------------------------------------\n");
            std::print(stream, "        Average  {} +/- {} [K]\n", 
                    prefactor * computedAverage.first(k,l).totalInter,
                    prefactor * computedAverage.second(k,l).totalInter);
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
        if(k==l)
        { 
          std::print(stream, "    Van der Waals energy {}-{} [{}-{}]:\n",k,l,
                  components[k].name, components[l].name);
        }
        else 
        {
          prefactor *= 2.0;
          std::print(stream, "    Van der Waals energy {}-{} + {}-{} [{}-{}]:\n",k,l,l,k,
                  components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < averageEnergies.bookKeepingEnergyStatus.size(); ++i)
        {
            EnergyStatus blockAverage = averageEnergies.averagedEnergy(i);
            std::print(stream, "        Block[ {:2d}] {}\n", i, prefactor * blockAverage(k,l).VanDerWaals);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {} +/- {} [K]\n", 
                prefactor * computedAverage.first(k,l).VanDerWaals,
                prefactor * computedAverage.second(k,l).VanDerWaals);
        std::print(stream, "\n");
      }
    }

    for (size_t k = 0; k < components.size(); k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if(k==l)
        { 
          std::print(stream, "    VDW Tail-Correction energy {}-{} [{}-{}]:\n",k,l,
                  components[k].name, components[l].name);
        }
        else 
        {
          prefactor *= 2.0;
          std::print(stream, "    VDW Tail-Correction energy {}-{} + {}-{} [{}-{}]:\n",k,l,l,k,
                  components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < averageEnergies.bookKeepingEnergyStatus.size(); ++i)
        {
            EnergyStatus blockAverage = averageEnergies.averagedEnergy(i);
            std::print(stream, "        Block[ {:2d}] {}\n", i, prefactor * blockAverage(k,l).VanDerWaalsTailCorrection);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {} +/- {} [K]\n", 
                prefactor * computedAverage.first(k,l).VanDerWaalsTailCorrection,
                prefactor * computedAverage.second(k,l).VanDerWaalsTailCorrection);
        std::print(stream, "\n");
      }
    }

    for (size_t k = 0; k < components.size(); k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if(k==l)
        { 
          std::print(stream, "    Coulomb Real energy {}-{} [{}-{}]:\n",k,l,
                  components[k].name, components[l].name);
        }
        else 
        {
          prefactor *= 2.0;
          std::print(stream, "    Coulomb Real energy {}-{} + {}-{} [{}-{}]:\n",k,l,l,k,
                  components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < averageEnergies.bookKeepingEnergyStatus.size(); ++i)
        {
            EnergyStatus blockAverage = averageEnergies.averagedEnergy(i);
            std::print(stream, "        Block[ {:2d}] {}\n", i, prefactor * blockAverage(k,l).CoulombicReal);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {} +/- {} [K]\n", 
                prefactor * computedAverage.first(k,l).CoulombicReal,
                prefactor * computedAverage.second(k,l).CoulombicReal);
        std::print(stream, "\n");
      }
    }

    for (size_t k = 0; k < components.size(); k++)
    {
      for (size_t l = k; l < components.size(); l++)
      {
        double prefactor = Units::EnergyToKelvin;
        if(k==l)
        { 
          std::print(stream, "    Coulomb Fourier energy {}-{} [{}-{}]:\n",k,l,
                  components[k].name, components[l].name);
        }
        else 
        {
          prefactor *= 2.0;
          std::print(stream, "    Coulomb Fourier energy {}-{} + {}-{} [{}-{}]:\n",k,l,l,k,
                  components[k].name, components[l].name);
        }
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        for (size_t i = 0; i < averageEnergies.bookKeepingEnergyStatus.size(); ++i)
        {
            EnergyStatus blockAverage = averageEnergies.averagedEnergy(i);
            std::print(stream, "        Block[ {:2d}] {}\n", i, prefactor * blockAverage(k,l).CoulombicFourier);
        }
        std::print(stream, "        -----------------------------------------------------------------------\n");
        std::print(stream, "        Average  {} +/- {} [K]\n", 
                prefactor * computedAverage.first(k,l).CoulombicFourier,
                prefactor * computedAverage.second(k,l).CoulombicFourier);
        std::print(stream, "\n");
      }
    }
    std::print(stream, "\n");

    std::print(stream, "Total energy:\n");
    std::print(stream, "-------------------------------------------------------------------------------\n");
    double prefactor = Units::EnergyToKelvin;
    for (size_t i = 0; i < averageEnergies.bookKeepingEnergyStatus.size(); ++i)
    {
        EnergyStatus blockAverage = averageEnergies.averagedEnergy(i);
        std::print(stream, "    Block[ {:2d}] {}\n", i, prefactor * blockAverage.totalEnergy);
    }
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    std::print(stream, "    Average  {} +/- {} [K]\n", prefactor * computedAverage.first.totalEnergy, prefactor * computedAverage.second.totalEnergy);

    std::print(stream, "\n");
    return stream.str();
}

