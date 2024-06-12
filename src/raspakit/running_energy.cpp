module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <map>
#include <functional>
#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>
#include <vector>
#include <array>
#include <ranges>
#include <algorithm>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#include <print>
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
import <print>;
#endif


import units;
import archive;
import stringutils;
import json;

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

std::string RunningEnergy::printMD(const std::string &label, double referenceEnergy) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Conserved energy:            {: .6e} [K]\n",   conv * conservedEnergy());
  std::print(stream, "Drift:                       {: .6e} [K]\n\n", std::abs(conv * (conservedEnergy() - referenceEnergy) / referenceEnergy));
  std::print(stream, "Total potential energy:      {: .6e} [K]\n",   conv * potentialEnergy());
  std::print(stream, "    external field VDW:      {: .6e} [K]\n",   conv * externalFieldVDW);
  std::print(stream, "    external field Real:     {: .6e} [K]\n",   conv * externalFieldCharge);
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n",   conv * frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n",   conv * frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n",   conv * moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n",   conv * moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n",   conv * tail);
  std::print(stream, "    Coulombic Ewald:         {: .6e} [K]\n",   conv * ewald);
  std::print(stream, "    intra VDW:               {: .6e} [K]\n",   conv * intraVDW);
  std::print(stream, "    intra Coulombic:         {: .6e} [K]\n",   conv * intraCoul);
  std::print(stream, "    polarization:            {: .6e} [K]\n",   conv * polarization);
  std::print(stream, "    dU/dlambda VDW:          {: .6e} [K]\n",   conv * dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real:         {: .6e} [K]\n",   conv * dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald:        {: .6e} [K]\n\n", conv * dudlambdaEwald);
  std::print(stream, "Total kinetic energy:        {: .6e} [K]\n",   conv * kineticEnergy());
  std::print(stream, "    translation kinetic:     {: .6e} [K]\n",   conv * translationalKineticEnergy);
  std::print(stream, "    rotational kinetic:      {: .6e} [K]\n",   conv * rotationalKineticEnergy);
  std::print(stream, "\n");

  return stream.str();
}

nlohmann::json RunningEnergy::jsonMC() const
{
  nlohmann::json status;

  double conv = Units::EnergyToKelvin;
  status["Total potential energy [K]"] = conv * potentialEnergy();
  status["external field VDW [K]"] = conv * externalFieldVDW;
  status["external field Real [K]"] = conv * externalFieldCharge;
  status["framework-molecule VDW [K]"] = conv * frameworkMoleculeVDW;
  status["framework-molecule Rea [K]"] = conv * frameworkMoleculeCharge;
  status["molecule-molecule VDW [K]"] = conv * moleculeMoleculeVDW;
  status["molecule-molecule Real [K]"] = conv * moleculeMoleculeCharge;
  status["Van der Waals (Tail) [K]"] = conv * tail;
  status["Coulombic Ewald [K]"] = conv * ewald;
  status["intra VDW [K]"] = conv * intraVDW;
  status["intra Coulombic [K]"] = conv * intraCoul;
  status["polarization [K]"] = conv * polarization;
  status["dU/dlambda VDW [K]"] = conv * dudlambdaVDW;
  status["dU/dlambda Real [K]"] = conv * dudlambdaCharge;
  status["dU/dlambda Ewald [K]"] = conv * dudlambdaEwald;

  return status;
}

nlohmann::json RunningEnergy::jsonMD() const
{
  nlohmann::json status = jsonMC();
  double conv = Units::EnergyToKelvin;
  status["Total kinetic energy [K]"] = conv * kineticEnergy();
  status["translational [K]"] = conv * translationalKineticEnergy;
  status["rotational [K]"] = conv * rotationalKineticEnergy;
  return status;
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
  if(versionNumber > e.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Component' at line {} in file {}\n",
                                         location.line(), location.file_name()));
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
