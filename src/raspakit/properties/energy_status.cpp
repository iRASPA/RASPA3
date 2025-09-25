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
#include <source_location>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#endif

module energy_status;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import stringutils;
import units;
import component;
import energy_status_intra;
import energy_status_inter;
import energy_factor;
import json;

std::string EnergyStatus::printEnergyStatus(const std::vector<Component> &components, const std::string &label)
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Total potential energy:  {: .6e}\n", conv * totalEnergy.energy);
  std::print(stream, "    framework-molecule Van der Waals:        {: .6e} [{}]\n",
             conv * frameworkMoleculeEnergy.VanDerWaals.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Van der Waals (Tail): {: .6e} [{}]\n",
             conv * frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Coulombic Real:       {: .6e} [{}]\n",
             conv * frameworkMoleculeEnergy.CoulombicReal.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Coulombic Fourier:    {: .6e} [{}]\n",
             conv * frameworkMoleculeEnergy.CoulombicFourier.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule  Van der Waals:        {: .6e} [{}]\n",
             conv * interEnergy.VanDerWaals.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule  Van der Waals (Tail): {: .6e} [{}]\n",
             conv * interEnergy.VanDerWaalsTailCorrection.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule  Coulombic Real:       {: .6e} [{}/-]\n",
             conv * interEnergy.CoulombicReal.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule  Coulombic Fourier:    {: .6e} [{}/-]\n",
             conv * interEnergy.CoulombicFourier.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    polarization:                            {: .6e} [{}/-]\n", conv * polarizationEnergy.energy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda:                              {: .6e} [{}/-]\n", conv * totalEnergy.dUdlambda,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "    translational kinetic energy:            {: .6e} [{}/-]\n", conv * translationalKineticEnergy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "    rotational kinetic energy:               {: .6e} [{}/-]\n", conv * rotationalKineticEnergy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "    Nose-Hoover energy:                      {: .6e} [{}/-]\n", conv * noseHooverEnergy,
             Units::displayedUnitOfEnergyString);

  for (std::size_t i = 0; i < components.size(); ++i)
  {
    std::print(stream, "    Component: {} [{}]\n", i, components[i].name);
    std::print(stream, "    ---------------------------------------------------------------------------\n\n");
    std::print(stream, "    Molecule bond:              {: .6e} [{}]\n", conv * intraComponentEnergies[i].bond,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule bend:              {: .6e} [{}]\n", conv * intraComponentEnergies[i].bend,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule inversionBend:     {: .6e} [{}]\n", conv * intraComponentEnergies[i].inversionBend,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule ureyBradley:       {: .6e} [{}]\n", conv * intraComponentEnergies[i].ureyBradley,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule torsion:           {: .6e} [{}]\n", conv * intraComponentEnergies[i].torsion,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule improperTorsion:   {: .6e} [{}]\n",
               conv * intraComponentEnergies[i].improperTorsion, Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule bondBond:          {: .6e} [{}]\n", conv * intraComponentEnergies[i].bondBond,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule bondBend:          {: .6e} [{}]\n", conv * intraComponentEnergies[i].bondBend,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule bondTorsion:       {: .6e} [{}]\n", conv * intraComponentEnergies[i].bondTorsion,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule bendBend:          {: .6e} [{}]\n", conv * intraComponentEnergies[i].bendBend,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule bendTorsion:       {: .6e} [{}]\n", conv * intraComponentEnergies[i].bendTorsion,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule intraVDW:          {: .6e} [{}]\n", conv * intraComponentEnergies[i].vanDerWaals,
               Units::displayedUnitOfEnergyString);
    std::print(stream, "    Molecule intraChargeCharge: {: .6e} [{}]\n\n", conv * intraComponentEnergies[i].coulomb,
               Units::displayedUnitOfEnergyString);

    for (std::size_t j = 0; j < components.size(); ++j)
    {
      std::print(stream, "    Component: {}-{} [{}-{}]\n", i, j, components[i].name, components[j].name);
      std::print(stream, "        Van der Waals:        {: .6e} [{}]\n",
                 conv * interComponentEnergies[i * numberOfComponents + j].VanDerWaals.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "        Van der Waals (Tail): {: .6e} [{}]\n",
                 conv * interComponentEnergies[i * numberOfComponents + j].VanDerWaalsTailCorrection.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "        Coulombic Real:       {: .6e} [{}]\n",
                 conv * interComponentEnergies[i * numberOfComponents + j].CoulombicReal.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "        Coulombic Fourier:    {: .6e} [{}]\n",
                 conv * interComponentEnergies[i * numberOfComponents + j].CoulombicFourier.energy,
                 Units::displayedUnitOfEnergyString);
      std::print(stream, "        -----------------------------------------------------------------------\n");
      std::print(stream, "        Sum                   {: .6e} [{}]\n\n",
                 conv * interComponentEnergies[i * numberOfComponents + j].totalInter.energy,
                 Units::displayedUnitOfEnergyString);
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
  archive << e.externalFieldMoleculeEnergy;
  archive << e.frameworkMoleculeEnergy;
  archive << e.interEnergy;
  archive << e.intraComponentEnergies;
  archive << e.externalFieldComponentEnergies;
  archive << e.frameworkComponentEnergies;
  archive << e.interComponentEnergies;
  archive << e.polarizationEnergy;
  archive << e.dUdlambda;
  archive << e.translationalKineticEnergy;
  archive << e.rotationalKineticEnergy;
  archive << e.noseHooverEnergy;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyStatus &e)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyStatus' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfExternalFields;
  archive >> e.numberOfFrameworks;
  archive >> e.numberOfComponents;
  archive >> e.totalEnergy;
  archive >> e.intraEnergy;
  archive >> e.externalFieldMoleculeEnergy;
  archive >> e.frameworkMoleculeEnergy;
  archive >> e.interEnergy;
  archive >> e.intraComponentEnergies;
  archive >> e.externalFieldComponentEnergies;
  archive >> e.frameworkComponentEnergies;
  archive >> e.interComponentEnergies;
  archive >> e.polarizationEnergy;
  archive >> e.dUdlambda;
  archive >> e.translationalKineticEnergy;
  archive >> e.rotationalKineticEnergy;
  archive >> e.noseHooverEnergy;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("EnergyStatus: Error in binary restart\n"));
  }
#endif

  return archive;
}
