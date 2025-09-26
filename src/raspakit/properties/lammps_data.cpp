module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <optional>
#include <print>
#include <source_location>
#include <span>
#include <streambuf>
#include <string>
#include <vector>
#endif

module write_lammps_data;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double3;
import atom;
import simulationbox;
import forcefield;
import component;
import lammps_io;
import framework;
import molecule;

WriteLammpsData::WriteLammpsData(std::size_t systemId, std::size_t sampleEvery)
    : sampleEvery(sampleEvery), systemId(systemId)
{
  std::filesystem::create_directory("lammps");
  std::ofstream stream(std::format("lammps/s{}.data", systemId));
}

void WriteLammpsData::update(std::size_t currentCycle, std::span<const Component> components,
                             std::span<const Atom> atomData, std::span<const Molecule> moleculeData,
                             const SimulationBox simulationBox, const ForceField forceField,
                             std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                             std::optional<Framework> framework)
{
  if (currentCycle % sampleEvery != 0) return;
  ;
  std::ofstream stream(std::format("lammps/s{}.data", systemId), std::ios_base::out);
  stream << IO::WriteLAMMPSDataFile(components, atomData, moleculeData, simulationBox, forceField,
                                    numberOfIntegerMoleculesPerComponent, framework)
         << std::endl;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const WriteLammpsData &m)
{
  archive << m.versionNumber;

  archive << m.sampleEvery;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, WriteLammpsData &m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'WriteLammpsData' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> m.sampleEvery;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("WriteLammpsData: Error in binary restart\n"));
  }
#endif

  return archive;
}
