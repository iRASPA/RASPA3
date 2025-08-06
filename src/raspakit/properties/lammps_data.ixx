module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>
#endif

export module write_lammps_data;

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

/**
 * \brief A property writer class for writing LAMMPS data files.
 * 
 * Needs IO::WriteLAMMPSDataFile, which does the actual writing, this class itself is just a property holding class
 * managing file names and writing frequency.
 */
export struct WriteLammpsData
{
  WriteLammpsData(std::size_t systemId, std::size_t sampleEvery);

  /**
   * \brief Write data to LAMMPS file.
   */
  void update(std::size_t currentCycle, std::span<const Component> components, std::span<const Atom> atomPositions,
              const SimulationBox simulationBox, const ForceField forceField,
              std::vector<std::size_t> numberOfIntegerMoleculesPerComponent, std::optional<Framework> framework);

  std::size_t sampleEvery{10}; ///< Writing frequency.
  std::size_t systemId; ///< Necessary for determining filename.

  int modelNumber{1};
};
