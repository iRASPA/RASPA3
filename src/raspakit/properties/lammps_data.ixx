module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>
#endif

export module write_lammps_data;

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
   */
  void update(std::size_t currentCycle, std::span<const Component> components, std::span<const Atom> atomData,
              std::span<const Molecule> moleculeData, const SimulationBox simulationBox, const ForceField forceField,
              std::vector<std::size_t> numberOfIntegerMoleculesPerComponent, std::optional<Framework> framework);

  std::size_t sampleEvery{10};  ///< Writing frequency.
  std::size_t systemId;         ///< Necessary for determining filename.

  int modelNumber{1};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const WriteLammpsData& box);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, WriteLammpsData& box);
};
