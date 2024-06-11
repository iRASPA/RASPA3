module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <ostream>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <complex>
#include <exception>
#include <source_location>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module running_energy;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <map>;
import <functional>;
import <iostream>;
import <sstream>;
import <ostream>;
import <fstream>;
import <vector>;
import <array>;
import <ranges>;
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

import units;
import archive;
import stringutils;

std::string RunningEnergy::printMC() const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Total potential energy:      {: .6e} [K]\n", conv * potentialEnergy());
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "    external field VDW:      {: .6e} [K]\n", conv * externalFieldVDW);
  std::print(stream, "    external field Real:     {: .6e} [K]\n", conv * externalFieldCharge);
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n", conv * frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n", conv * frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n", conv * moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n", conv * moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n", conv * tail);
  std::print(stream, "    Coulombic Ewald:         {: .6e} [K]\n", conv * ewald);
  std::print(stream, "    intra VDW:               {: .6e} [K]\n", conv * intraVDW);
  std::print(stream, "    intra Coulombic:         {: .6e} [K]\n", conv * intraCoul);
  std::print(stream, "    polarization:            {: .6e} [K]\n", conv * polarization);
  std::print(stream, "    dU/dlambda VDW:          {: .6e} [K]\n", conv * dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real:         {: .6e} [K]\n", conv * dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald:        {: .6e} [K]\n", conv * dudlambdaEwald);
  std::print(stream, "\n");

  return stream.str();
}

std::string RunningEnergy::printMD() const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Total potential energy:      {: .6e} [K]\n", conv * potentialEnergy());
  std::print(stream, "    external field VDW:      {: .6e} [K]\n", conv * externalFieldVDW);
  std::print(stream, "    external field Real:     {: .6e} [K]\n", conv * externalFieldCharge);
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n", conv * frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n", conv * frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n", conv * moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n", conv * moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n", conv * tail);
  std::print(stream, "    Coulombic Ewald:         {: .6e} [K]\n", conv * ewald);
  std::print(stream, "    intra VDW:               {: .6e} [K]\n", conv * intraVDW);
  std::print(stream, "    intra Coulombic:         {: .6e} [K]\n", conv * intraCoul);
  std::print(stream, "    polarization:            {: .6e} [K]\n", conv * polarization);
  std::print(stream, "    dU/dlambda VDW:          {: .6e} [K]\n", conv * dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real:         {: .6e} [K]\n", conv * dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald:        {: .6e} [K]\n", conv * dudlambdaEwald);
  std::print(stream, "\n");
  std::print(stream, "Total kinetic energy:        {: .6e} [K]\n", conv * kineticEnergy());
  std::print(stream, "    translational:           {: .6e} [K]\n", conv * translationalKineticEnergy);
  std::print(stream, "    rotational:              {: .6e} [K]\n", conv * rotationalKineticEnergy);
  std::print(stream, "\n");

  return stream.str();
}

std::string RunningEnergy::printMC(const std::string &label) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Total potential energy:      {: .6e} [K]\n", conv * potentialEnergy());
  std::print(stream, "    external field VDW:      {: .6e} [K]\n", conv * externalFieldVDW);
  std::print(stream, "    external field Real:     {: .6e} [K]\n", conv * externalFieldCharge);
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n", conv * frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n", conv * frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n", conv * moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n", conv * moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n", conv * tail);
  std::print(stream, "    Coulombic Ewald:         {: .6e} [K]\n", conv * ewald);
  std::print(stream, "    intra VDW:               {: .6e} [K]\n", conv * intraVDW);
  std::print(stream, "    intra Coulombic:         {: .6e} [K]\n", conv * intraCoul);
  std::print(stream, "    polarization:            {: .6e} [K]\n", conv * polarization);
  std::print(stream, "    dU/dlambda VDW:          {: .6e} [K]\n", conv * dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real:         {: .6e} [K]\n", conv * dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald:        {: .6e} [K]\n", conv * dudlambdaEwald);
  std::print(stream, "\n");

  return stream.str();
}

void RunningEnergy::logMC(HDF5Writer &hdf5, const std::string &label) const
{
  double conv = Units::EnergyToKelvin;
  std::string group = "statistics/" + label;
  hdf5.createGroup(group);
  hdf5.writeMetaInfo(group, "Total potential energy: [K]", conv * potentialEnergy());
  hdf5.writeMetaInfo(group, "external field VDW: [K]", conv * externalFieldVDW);
  hdf5.writeMetaInfo(group, "external field Real: [K]", conv * externalFieldCharge);
  hdf5.writeMetaInfo(group, "framework-molecule VDW: [K]", conv * frameworkMoleculeVDW);
  hdf5.writeMetaInfo(group, "framework-molecule Real [K]", conv * frameworkMoleculeCharge);
  hdf5.writeMetaInfo(group, "molecule-molecule VDW: [K]", conv * moleculeMoleculeVDW);
  hdf5.writeMetaInfo(group, "molecule-molecule Real: [K]", conv * moleculeMoleculeCharge);
  hdf5.writeMetaInfo(group, "Van der Waals (Tail): [K]", conv * tail);
  hdf5.writeMetaInfo(group, "Coulombic Ewald: [K]", conv * ewald);
  hdf5.writeMetaInfo(group, "intra VDW: [K]", conv * intraVDW);
  hdf5.writeMetaInfo(group, "intra Coulombic: [K]", conv * intraCoul);
  hdf5.writeMetaInfo(group, "polarization: [K]", conv * polarization);
  hdf5.writeMetaInfo(group, "dU/dlambda VDW: [K]", conv * dudlambdaVDW);
  hdf5.writeMetaInfo(group, "dU/dlambda Real: [K]", conv * dudlambdaCharge);
  hdf5.writeMetaInfo(group, "dU/dlambda Ewald: [K]", conv * dudlambdaEwald);
}

std::string RunningEnergy::printMD(const std::string &label, double referenceEnergy) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Conserved energy:            {: .6e} [K]\n", conv * conservedEnergy());
  std::print(stream, "Drift:                       {: .6e} [K]\n\n",
             std::abs(conv * (conservedEnergy() - referenceEnergy) / referenceEnergy));
  std::print(stream, "Total potential energy:      {: .6e} [K]\n", conv * potentialEnergy());
  std::print(stream, "    external field VDW:      {: .6e} [K]\n", conv * externalFieldVDW);
  std::print(stream, "    external field Real:     {: .6e} [K]\n", conv * externalFieldCharge);
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n", conv * frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n", conv * frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n", conv * moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n", conv * moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n", conv * tail);
  std::print(stream, "    Coulombic Ewald:         {: .6e} [K]\n", conv * ewald);
  std::print(stream, "    intra VDW:               {: .6e} [K]\n", conv * intraVDW);
  std::print(stream, "    intra Coulombic:         {: .6e} [K]\n", conv * intraCoul);
  std::print(stream, "    polarization:            {: .6e} [K]\n", conv * polarization);
  std::print(stream, "    dU/dlambda VDW:          {: .6e} [K]\n", conv * dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real:         {: .6e} [K]\n", conv * dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald:        {: .6e} [K]\n\n", conv * dudlambdaEwald);
  std::print(stream, "Total kinetic energy:        {: .6e} [K]\n", conv * kineticEnergy());
  std::print(stream, "    translation kinetic:     {: .6e} [K]\n", conv * translationalKineticEnergy);
  std::print(stream, "    rotational kinetic:      {: .6e} [K]\n", conv * rotationalKineticEnergy);
  std::print(stream, "\n");

  return stream.str();
}

void RunningEnergy::logMD(HDF5Writer &hdf5, const std::string &label, double referenceEnergy) const
{
  double conv = Units::EnergyToKelvin;
  std::string group = "statistics/" + label;
  logMC(hdf5, label);
  hdf5.writeMetaInfo(group, "Conserved energy: [K]", conv * conservedEnergy());
  hdf5.writeMetaInfo(group, "Drift: [K]", std::abs(conv * (conservedEnergy() - referenceEnergy) / referenceEnergy));
  hdf5.writeMetaInfo(group, "Total kinetic energy: [K]", conv * kineticEnergy());
  hdf5.writeMetaInfo(group, "translation kinetic: [K]", conv * translationalKineticEnergy);
  hdf5.writeMetaInfo(group, "rotational kinetic: [K]", conv * rotationalKineticEnergy);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RunningEnergy &e)
{
  archive << e.versionNumber;

  archive << e.externalFieldVDW;
  archive << e.frameworkMoleculeVDW;
  archive << e.moleculeMoleculeVDW;
  archive << e.externalFieldCharge;
  archive << e.frameworkMoleculeCharge;
  archive << e.moleculeMoleculeCharge;
  archive << e.ewald;
  archive << e.intraVDW;
  archive << e.intraCoul;
  archive << e.tail;
  archive << e.polarization;
  archive << e.dudlambdaVDW;
  archive << e.dudlambdaCharge;
  archive << e.dudlambdaEwald;
  archive << e.translationalKineticEnergy;
  archive << e.rotationalKineticEnergy;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RunningEnergy &e)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Component' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> e.externalFieldVDW;
  archive >> e.frameworkMoleculeVDW;
  archive >> e.moleculeMoleculeVDW;
  archive >> e.externalFieldCharge;
  archive >> e.frameworkMoleculeCharge;
  archive >> e.moleculeMoleculeCharge;
  archive >> e.ewald;
  archive >> e.intraVDW;
  archive >> e.intraCoul;
  archive >> e.tail;
  archive >> e.polarization;
  archive >> e.dudlambdaVDW;
  archive >> e.dudlambdaCharge;
  archive >> e.dudlambdaEwald;
  archive >> e.translationalKineticEnergy;
  archive >> e.rotationalKineticEnergy;

  return archive;
}
