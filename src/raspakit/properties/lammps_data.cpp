module;

module write_lammps_data;

import std;

import archive;
import double3;
import atom;
import atom_dynamics;
import simulationbox;
import forcefield;
import component;
import lammps_io;
import lammps_topology;
import framework;
import molecule;

WriteLammpsData::WriteLammpsData(std::size_t systemId, std::size_t sampleEvery)
    : sampleEvery(sampleEvery), systemId(systemId)
{
  std::filesystem::create_directory("lammps");
  std::ofstream stream(std::format("lammps/s{}.data", systemId));
}

void WriteLammpsData::update(std::size_t currentCycle, std::span<const Component> components,
                             std::span<const Atom> atomData, std::span<const AtomDynamics> atomDynamics,
                             std::span<const Molecule> moleculeData,
                             const SimulationBox simulationBox, const ForceField forceField,
                             std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                             std::optional<Framework> framework)
{
  if (currentCycle % sampleEvery != 0) return;

  LAMMPS::ExportOptions options{};
  options.dataFile = std::format("s{}.data", systemId);
  options.tableFile = std::format("s{}.table", systemId);
  options.pairListFile = std::format("s{}.pairs", systemId);
  LAMMPS::ExportFiles files =
      LAMMPS::exportSystem(components, atomData, atomDynamics, moleculeData, simulationBox, forceField,
                           numberOfIntegerMoleculesPerComponent, framework, options);

  std::ofstream(std::format("lammps/{}", options.dataFile), std::ios_base::out) << files.data;
  std::ofstream(std::format("lammps/s{}.in", systemId), std::ios_base::out) << files.input;
  if (!files.table.empty())
  {
    std::ofstream(std::format("lammps/{}", options.tableFile), std::ios_base::out) << files.table;
  }
  if (!files.pairList.empty())
  {
    std::ofstream(std::format("lammps/{}", options.pairListFile), std::ios_base::out) << files.pairList;
  }
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
