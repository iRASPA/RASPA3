module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <ostream>
#include <print>
#include <ranges>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
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
  std::print(stream, "Total potential energy{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * potentialEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "    external field VDW{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field Real{}      {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule VDW{}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Real{}  {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule VDW{}    {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule Real{}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Van der Waals (Tail){}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * tail, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald Fourier{}            {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_fourier, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald self{}               {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_self, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald exclusion{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_exclusion, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra VDW{}                {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra Coulombic{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraCoul, Units::displayedUnitOfEnergyString);
  std::print(stream, "    polarization{}             {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * polarization, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda VDW{}           {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Real{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Ewald{}         {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaEwald, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  return stream.str();
}

std::string RunningEnergy::printMCDiff(RunningEnergy &other) const
{
  std::ostringstream stream;
  RunningEnergy drift = *this - other;
  double conv = Units::EnergyToKelvin;

  std::print(stream, "\n");
  std::print(stream, "Energy statistics             |  Energy [{}]   | Recomputed [{}]|  Drift [{}]    |\n",
             Units::displayedUnitOfEnergyString, Units::displayedUnitOfEnergyString,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "===============================================================================\n");
  std::print(stream, "Total potential energy{}     | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * potentialEnergy(), conv * other.potentialEnergy(),
             conv * drift.potentialEnergy());
  std::print(stream, "    external field VDW{}     | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * externalFieldVDW, conv * other.externalFieldVDW,
             conv * drift.externalFieldVDW);
  std::print(stream, "    external field Real{}    | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * externalFieldCharge, conv * other.externalFieldCharge,
             conv * drift.externalFieldCharge);
  std::print(stream, "    framework-molecule VDW{} | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * frameworkMoleculeVDW,
             conv * other.frameworkMoleculeVDW, conv * drift.frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real{}| {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * frameworkMoleculeCharge,
             conv * other.frameworkMoleculeCharge, conv * drift.frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW{}  | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * moleculeMoleculeVDW, conv * other.moleculeMoleculeVDW,
             conv * drift.moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real{} | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * moleculeMoleculeCharge,
             conv * other.moleculeMoleculeCharge, conv * drift.moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail){}   | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * tail, conv * other.tail, conv * drift.tail);
  std::print(stream, "    Ewald Fourier{}          | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * ewald_fourier, conv * other.ewald_fourier,
             conv * drift.ewald_fourier);
  std::print(stream, "    Ewald self{}             | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * ewald_self, conv * other.ewald_self,
             conv * drift.ewald_self);
  std::print(stream, "    Ewald exclusion{}        | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * ewald_exclusion, conv * other.ewald_exclusion,
             conv * drift.ewald_exclusion);
  std::print(stream, "    intra VDW{}              | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * intraVDW, conv * other.intraVDW,
             conv * drift.intraVDW);
  std::print(stream, "    intra Coulombic{}        | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * intraCoul, conv * other.intraCoul,
             conv * drift.intraCoul);
  std::print(stream, "    polarization{}           | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * polarization, conv * other.polarization,
             conv * drift.polarization);
  std::print(stream, "    dU/dlambda VDW{}         | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * dudlambdaVDW, conv * other.dudlambdaVDW,
             conv * drift.dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real{}        | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * dudlambdaCharge, conv * other.dudlambdaCharge,
             conv * drift.dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald{}       | {: 10.6e} | {: 10.6e} | {: 10.6e} |\n",
             Units::displayedUnitOfEnergyConversionString, conv * dudlambdaEwald, conv * other.dudlambdaEwald,
             conv * drift.dudlambdaEwald);
  std::print(stream, "-------------------------------------------------------------------------------\n");

  return stream.str();
}

std::string RunningEnergy::printMD() const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Total potential energy{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * potentialEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field VDW{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field Real{}      {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule VDW{}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Real{}  {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule VDW{}    {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule Real{}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Van der Waals (Tail){}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * tail, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald Fourier{}            {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_fourier, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald self{}               {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_self, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald exclusion{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_exclusion, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra VDW{}                {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra Coulombic{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraCoul, Units::displayedUnitOfEnergyString);
  std::print(stream, "    polarization{}             {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * polarization, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda VDW{}           {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Real{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Ewald{}         {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaEwald, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");
  std::print(stream, "Total kinetic energy{}         {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * kineticEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "    translational{}            {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * translationalKineticEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    rotational{}               {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * rotationalKineticEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "Baro/Thermostat energy{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * NoseHooverEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  return stream.str();
}

std::string RunningEnergy::printMC(const std::string &label) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Total potential energy{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * potentialEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field VDW{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field Real{}      {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule VDW{}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Real{}  {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule VDW{}    {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule Real{}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Van der Waals (Tail){}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * tail, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald Fourier{}            {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_fourier, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald self{}               {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_self, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald exclusion{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_exclusion, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra VDW{}                {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra Coulombic{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraCoul, Units::displayedUnitOfEnergyString);
  std::print(stream, "    polarization{}             {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * polarization, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda VDW{}           {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Real{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Ewald{}         {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaEwald, Units::displayedUnitOfEnergyString);
  std::print(stream, "\n");

  return stream.str();
}

std::string RunningEnergy::printMD(const std::string &label, double referenceEnergy) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Conserved energy{}           {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * conservedEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "Drift{}                      {: .6e} [{}]\n\n", Units::displayedUnitOfEnergyConversionString,
             std::abs(conv * ((conservedEnergy() - referenceEnergy) / referenceEnergy)),
             Units::displayedUnitOfEnergyString);
  std::print(stream, "Total potential energy{}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * potentialEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field VDW{}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    external field Real{}    {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * externalFieldCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule VDW{} {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    framework-molecule Real{}{: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * frameworkMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule VDW{}  {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    molecule-molecule Real{} {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * moleculeMoleculeCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Van der Waals (Tail){}   {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * tail, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald Fourier{}          {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_fourier, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald self{}             {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_self, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Ewald exclusion{}        {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * ewald_exclusion, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra VDW{}              {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    intra Coulombic{}        {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * intraCoul, Units::displayedUnitOfEnergyString);
  std::print(stream, "    polarization{}           {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * polarization, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda VDW{}         {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaVDW, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Real{}        {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaCharge, Units::displayedUnitOfEnergyString);
  std::print(stream, "    dU/dlambda Ewald{}       {: .6e} [{}]\n\n", Units::displayedUnitOfEnergyConversionString,
             conv * dudlambdaEwald, Units::displayedUnitOfEnergyString);
  std::print(stream, "Total kinetic energy{}       {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * kineticEnergy(), Units::displayedUnitOfEnergyString);
  std::print(stream, "    translation kinetic{}    {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * translationalKineticEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    rotational kinetic{}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * rotationalKineticEnergy, Units::displayedUnitOfEnergyString);
  std::print(stream, "Baro/thermostat energy{}     {: .6e} [{}]\n", Units::displayedUnitOfEnergyConversionString,
             conv * NoseHooverEnergy, Units::displayedUnitOfEnergyString);
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
  status["Ewald Fourier [K]"] = conv * ewald_fourier;
  status["Ewald self [K]"] = conv * ewald_self;
  status["Ewald exclusion [K]"] = conv * ewald_exclusion;
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
  archive << e.ewald_fourier;
  archive << e.ewald_self;
  archive << e.ewald_exclusion;
  archive << e.intraVDW;
  archive << e.intraCoul;
  archive << e.tail;
  archive << e.polarization;
  archive << e.dudlambdaVDW;
  archive << e.dudlambdaCharge;
  archive << e.dudlambdaEwald;
  archive << e.translationalKineticEnergy;
  archive << e.rotationalKineticEnergy;
  archive << e.NoseHooverEnergy;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

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
  archive >> e.ewald_fourier;
  archive >> e.ewald_self;
  archive >> e.ewald_exclusion;
  archive >> e.intraVDW;
  archive >> e.intraCoul;
  archive >> e.tail;
  archive >> e.polarization;
  archive >> e.dudlambdaVDW;
  archive >> e.dudlambdaCharge;
  archive >> e.dudlambdaEwald;
  archive >> e.translationalKineticEnergy;
  archive >> e.rotationalKineticEnergy;
  archive >> e.NoseHooverEnergy;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("RunningEnergy: Error in binary restart\n"));
  }
#endif

  return archive;
}
