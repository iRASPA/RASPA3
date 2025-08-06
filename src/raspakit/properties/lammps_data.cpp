module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <span>
#include <streambuf>
#include <string>
#include <vector>
#endif

module write_lammps_data;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import atom;
import simulationbox;
import forcefield;
import component;
import lammps_io;
import framework;

WriteLammpsData::WriteLammpsData(std::size_t systemId, std::size_t sampleEvery)
    : sampleEvery(sampleEvery), systemId(systemId)
{
  std::filesystem::create_directory("lammps");
  std::ofstream stream(std::format("lammps/s{}.data", systemId));
}

void WriteLammpsData::update(std::size_t currentCycle, std::span<const Component> components,
                             std::span<const Atom> atomPositions, const SimulationBox simulationBox,
                             const ForceField forceField, std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                             std::optional<Framework> framework)
{
  if (currentCycle % sampleEvery != 0) return;
  ;
  std::ofstream stream(std::format("lammps/s{}.data", systemId), std::ios_base::out);
  stream << IO::WriteLAMMPSDataFile(components, atomPositions, simulationBox, forceField,
                                    numberOfIntegerMoleculesPerComponent, framework)
         << std::endl;
}
