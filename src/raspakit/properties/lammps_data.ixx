module;

export module write_lammps_data;

import std;

import archive;
import double3;
import atom;
import atom_dynamics;
import simulationbox;
import forcefield;
import component;
import lammps_io;
import framework;
import molecule;

/**
 * \brief A property writer class for writing LAMMPS data files.
 *
 * Needs IO::WriteLAMMPSDataFile, which does the actual writing, this class itself is just a property holding class
 * managing file names and writing frequency.
 */
export struct WriteLammpsData
{
  std::uint64_t versionNumber{1};

  WriteLammpsData() {};

  WriteLammpsData(std::size_t systemId, std::size_t sampleEvery);

  /**
   * \brief Write data to LAMMPS file.
   *
   * The files are named after 'systemId', the index of the system that writes (and not the 'systemId' the writer was
   * made for: the replicas of a parallel-tempering, TMMC or thermodynamic-integration driver are copies of one system,
   * and wrote to one set of files from as many threads).
   */
  void update(std::size_t systemId, std::size_t currentCycle, std::span<const Component> components,
              std::span<const Atom> atomData, std::span<const std::uint32_t> atomTags,
              std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
              std::span<const GroupState> groupData, std::span<const GroupState> frameworkGroupData,
              const SimulationBox simulationBox, const ForceField forceField,
              std::vector<std::size_t> numberOfIntegerMoleculesPerComponent, std::optional<Framework> framework);

  std::size_t sampleEvery{10};  ///< Writing frequency.
  std::size_t systemId{0};      ///< The system the writer was made for.

  int modelNumber{1};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const WriteLammpsData& box);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, WriteLammpsData& box);
};
