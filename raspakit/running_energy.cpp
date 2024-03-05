module;

module running_energy;

import units;

import <string>;
import <iostream>;
import <sstream>;
import <vector>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

import archive;
import stringutils;

void RunningEnergy::print(std::ostream &stream, const std::string &label)
{
  double conv = Units::EnergyToKelvin;
  std::print(stream, "Energy status {}\n", label);
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Total potential energy:      {: .6e}\n", conv * total());
  std::print(stream, "    external field VDW:      {: .6e}\n", conv * externalFieldVDW);
  std::print(stream, "    external field Real:     {: .6e}\n", conv * externalFieldCharge);
  std::print(stream, "    framework-molecule VDW:  {: .6e}\n", conv * frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e}\n", conv * frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e}\n", conv * moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e}\n", conv * moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e}\n", conv * tail);
  std::print(stream, "    Coulombic Ewald:         {: .6e}\n", conv * ewald);
  std::print(stream, "    intra VDW:               {: .6e}\n", conv * intraVDW);
  std::print(stream, "    intra Coulombic:         {: .6e}\n", conv * intraCoul);
  std::print(stream, "    polarization:            {: .6e}\n", conv * polarization);
  std::print(stream, "    dU/dlambda VDW:          {: .6e}\n", conv * dudlambdaVDW);
  std::print(stream, "    dU/dlambda Real:         {: .6e}\n", conv * dudlambdaCharge);
  std::print(stream, "    dU/dlambda Ewald:        {: .6e}\n", conv * dudlambdaEwald);
  std::print(stream, "\n");
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

  return archive;
}
