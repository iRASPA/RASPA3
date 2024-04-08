module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#include <type_traits>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module energy_status;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <iostream>;
import <sstream>;
import <fstream>;
import <vector>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <type_traits>;
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
import component;
import energy_status_intra;
import energy_status_inter;
import energy_factor;


std::string EnergyStatus::printEnergyStatus(const std::vector<Component>& components, const std::string &label)
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Total potential energy:  {: .6e}\n", conv * totalEnergy.energy);
  std::print(stream, "    framework-molecule Van der Waals:        {: .6e}\n", conv * frameworkMoleculeEnergy.VanDerWaals.energy);
  std::print(stream, "    framework-molecule Van der Waals (Tail): {: .6e}\n", conv * frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy);
  std::print(stream, "    framework-molecule Coulombic Real:       {: .6e}\n", conv * frameworkMoleculeEnergy.CoulombicReal.energy);
  std::print(stream, "    framework-molecule Coulombic Fourier:    {: .6e}\n", conv * frameworkMoleculeEnergy.CoulombicFourier.energy);
  std::print(stream, "    molecule-molecule  Van der Waals:        {: .6e}\n", conv * interEnergy.VanDerWaals.energy);
  std::print(stream, "    molecule-molecule  Van der Waals (Tail): {: .6e}\n", conv * interEnergy.VanDerWaalsTailCorrection.energy);
  std::print(stream, "    molecule-molecule  Coulombic Real:       {: .6e}\n", conv * interEnergy.CoulombicReal.energy);
  std::print(stream, "    molecule-molecule  Coulombic Fourier:    {: .6e}\n", conv * interEnergy.CoulombicFourier.energy);
  std::print(stream, "    dU/dlambda:                              {: .6e}\n", conv * totalEnergy.dUdlambda);
  
  for (size_t i = 0; i < components.size(); ++i)
  {
    std::print(stream, "    Component: {} [{}]\n", i, components[i].name);
    std::print(stream, "    ---------------------------------------------------------------------------\n\n");
    std::print(stream, "    Molecule bond:             {: .6e}\n", conv * intraComponentEnergies[i].bond);
    std::print(stream, "    Molecule bend:             {: .6e}\n", conv * intraComponentEnergies[i].bend);
    std::print(stream, "    Molecule inversionBend     {: .6e}\n", conv * intraComponentEnergies[i].inversionBend);
    std::print(stream, "    Molecule ureyBradley       {: .6e}\n", conv * intraComponentEnergies[i].ureyBradley);
    std::print(stream, "    Molecule torsion           {: .6e}\n", conv * intraComponentEnergies[i].torsion);
    std::print(stream, "    Molecule improperTorsion   {: .6e}\n", conv * intraComponentEnergies[i].improperTorsion);
    std::print(stream, "    Molecule bondBond          {: .6e}\n", conv * intraComponentEnergies[i].bondBond);
    std::print(stream, "    Molecule bondBend          {: .6e}\n", conv * intraComponentEnergies[i].bondBend);
    std::print(stream, "    Molecule bondTorsion       {: .6e}\n", conv * intraComponentEnergies[i].bondTorsion);
    std::print(stream, "    Molecule bendBend          {: .6e}\n", conv * intraComponentEnergies[i].bendBend);
    std::print(stream, "    Molecule bendTorsion       {: .6e}\n", conv * intraComponentEnergies[i].bendTorsion);
    std::print(stream, "    Molecule intraVDW          {: .6e}\n", conv * intraComponentEnergies[i].intraVDW);
    std::print(stream, "    Molecule intraChargeCharge {: .6e}\n\n", conv * intraComponentEnergies[i].intraChargeCharge);
    for (size_t j = 0; j < components.size(); ++j)
    {
      std::print(stream, "    Component: {}-{} [{}-{}]\n", i, j, components[i].name, components[j].name);
      std::print(stream, "        Van der Waals:        {: .6e}\n", conv * interComponentEnergies[i * numberOfComponents + j].VanDerWaals.energy);
      std::print(stream, "        Van der Waals (Tail): {: .6e}\n", conv * 
                                                                    interComponentEnergies[i * numberOfComponents + j].VanDerWaalsTailCorrection.energy);
      std::print(stream, "        Coulombic Real:       {: .6e}\n", conv * interComponentEnergies[i * numberOfComponents + j].CoulombicReal.energy);
      std::print(stream, "        Coulombic Fourier:    {: .6e}\n", conv * interComponentEnergies[i * numberOfComponents + j].CoulombicFourier.energy);
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Sum                   {: .6e}\n\n", conv * interComponentEnergies[i * numberOfComponents + j].totalInter.energy);
    }
  }
  std::print(stream, "\n");
  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyStatus &e)
{
  archive << e.versionNumber;

  archive << e.numberOfExternalFields;
  archive << e.numberOfFrameworks;
  archive << e.numberOfComponents;
  archive << e.totalEnergy;
  archive << e.intraEnergy;
  archive << e.interEnergy;
  archive << e.intraComponentEnergies;
  archive << e.interComponentEnergies;
  archive << e.dUdlambda;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyStatus &e)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > e.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyStatus' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfExternalFields;
  archive >> e.numberOfFrameworks;
  archive >> e.numberOfComponents;
  archive >> e.totalEnergy;
  archive >> e.intraEnergy;
  archive >> e.interEnergy;
  archive >> e.intraComponentEnergies;
  archive >> e.interComponentEnergies;
  archive >> e.dUdlambda;

  return archive;
}
