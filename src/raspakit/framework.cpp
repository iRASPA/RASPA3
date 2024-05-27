module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <ostream>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

#if !defined(_WIN32)
#include <assert.h>
#endif

module framework;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <ostream>;
import <sstream>;
import <algorithm>;
import <vector>;
import <array>;
import <map>;
import <string>;
import <span>;
import <optional>;
import <filesystem>;
import <fstream>;
import <cstdlib>;
import <exception>;
import <iterator>;
import <chrono>;
import <cstddef>;
#if defined(_WIN32)
import <cassert>;
#endif
import <exception>;
import <source_location>;
import <complex>;
import <type_traits>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import archive;
import int3;
import double3;
import double3x3;
import skposcarparser;
import characterset;
import stringutils;
import skparser;
import skposcarparser;
import skstructure;
import skasymmetricatom;
import skatomcopy;
import skcell;
import skspacegroup;
import forcefield;
import atom;
import multi_site_isotherm;
import simulationbox;
import cif_reader;
import move_statistics;
import bond_potential;
import hdf5;

// default constructor, needed for binary restart-file
Framework::Framework() {}

// create Component in 'inputreader.cpp'
Framework::Framework(size_t currentFramework, const ForceField& forceField, const std::string& componentName,
                     std::optional<const std::string> fileName, int3 numberOfUnitCells) noexcept(false)
    : numberOfUnitCells(numberOfUnitCells), frameworkId(currentFramework), name(componentName), filenameData(fileName)
{
  if (filenameData.has_value())
  {
    readFramework(forceField, filenameData.value());

    for (size_t i = 0; i < unitCellAtoms.size(); ++i)
    {
      unitCellAtoms[i].componentId = static_cast<uint8_t>(currentFramework);
      unitCellAtoms[i].moleculeId = 0;
    }
  }
}

// create programmatically an 'framework' component
Framework::Framework(size_t frameworkId, const ForceField& forceField, std::string fileName,
                     SimulationBox simulationBox, size_t spaceGroupHallNumber, std::vector<Atom> definedAtoms,
                     int3 numberOfUnitCells) noexcept(false)
    : simulationBox(simulationBox),
      spaceGroupHallNumber(spaceGroupHallNumber),
      numberOfUnitCells(numberOfUnitCells),
      frameworkId(frameworkId),
      name(fileName),
      filenameData(fileName),
      definedAtoms(definedAtoms)
{
  expandDefinedAtomsToUnitCell();
  makeSuperCell();

  unitCellMass = 0.0;
  for (const Atom& atom : unitCellAtoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    unitCellMass += forceField.pseudoAtoms[atomType].mass;
  }

  mass = 0.0;
  for (const Atom& atom : atoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    mass += forceField.pseudoAtoms[atomType].mass;
  }

  for (size_t i = 0; i < unitCellAtoms.size(); ++i)
  {
    unitCellAtoms[i].componentId = static_cast<uint8_t>(frameworkId);
    unitCellAtoms[i].moleculeId = 0;
  }
}

void Framework::readFramework(const ForceField& forceField, const std::string& fileName)
{
  const char* env_p = std::getenv("RASPA_DIR");

  const std::string frameworkFileName = fileName + ".cif";

  std::filesystem::path frameworkPathfile = std::filesystem::path(frameworkFileName);
  if (!std::filesystem::exists(frameworkPathfile)) frameworkPathfile = std::filesystem::path(env_p) / frameworkFileName;

  if (!std::filesystem::exists(frameworkPathfile))
    throw std::runtime_error(std::format("File '{}' not found\n", frameworkFileName));

  std::ifstream t(frameworkPathfile);
  std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

  CIFReader parser = CIFReader(fileContent, forceField);
  simulationBox = parser.simulationBox;
  definedAtoms = parser.fractionalAtoms;
  spaceGroupHallNumber = parser._spaceGroupHallNumber.value_or(1);

  // expand the fractional atoms based on the space-group
  expandDefinedAtomsToUnitCell();

  makeSuperCell();

  unitCellMass = 0.0;
  for (const Atom& atom : unitCellAtoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    unitCellMass += forceField.pseudoAtoms[atomType].mass;
  }

  mass = 0.0;
  for (const Atom& atom : atoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    mass += forceField.pseudoAtoms[atomType].mass;
  }

  for (size_t i = 0; i < unitCellAtoms.size(); ++i)
  {
    unitCellAtoms[i].componentId = static_cast<uint8_t>(frameworkId);
    unitCellAtoms[i].moleculeId = static_cast<uint8_t>(frameworkId);
  }
}

void Framework::expandDefinedAtomsToUnitCell()
{
  SKSpaceGroup spaceGroup = SKSpaceGroup(spaceGroupHallNumber);

  // expand the fractional atoms based on the space-group
  std::vector<Atom> expandAtoms{};
  expandAtoms.reserve(definedAtoms.size() * 256);

  for (Atom atomCopy : definedAtoms)
  {
    std::vector<double3> listOfPositions = spaceGroup.listOfSymmetricPositions(atomCopy.position);
    for (const double3& pos : listOfPositions)
    {
      atomCopy.position = simulationBox.cell * pos.fract();
      expandAtoms.push_back(atomCopy);
    }
  }

  // eliminate duplicates
  unitCellAtoms.clear();
  for (size_t i = 0; i < expandAtoms.size(); ++i)
  {
    bool overLap = false;
    for (size_t j = i + 1; j < expandAtoms.size(); ++j)
    {
      double3 dr = expandAtoms[i].position - expandAtoms[j].position;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      if (rr < 0.1)
      {
        overLap = true;
        break;
      }
    }
    if (!overLap)
    {
      unitCellAtoms.push_back(expandAtoms[i]);
    }
  }
}

void Framework::makeSuperCell()
{
  for (int32_t i = 0; i < numberOfUnitCells.x; ++i)
  {
    for (int32_t j = 0; j < numberOfUnitCells.y; ++j)
    {
      for (int32_t k = 0; k < numberOfUnitCells.z; ++k)
      {
        for (const Atom& atom : unitCellAtoms)
        {
          Atom atomCopy = atom;
          atomCopy.position +=
              simulationBox.cell * double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          atoms.push_back(atomCopy);
        }
      }
    }
  }
}

std::string Framework::printStatus(const ForceField& forceField) const
{
  std::ostringstream stream;

  std::print(stream, "Framework {} [{}]\n\n", frameworkId, name);

  std::print(stream, "    number Of Atoms:  {}\n", unitCellAtoms.size());
  std::print(stream, "    Mass:             {} [amu]\n", mass);

  for (size_t i = 0; i != definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);

    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomTypeString,
               definedAtoms[i].position.x, definedAtoms[i].position.y, definedAtoms[i].position.z,
               definedAtoms[i].charge);
  }

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    std::print(stream, "        {}", bonds[i].print());
  }
  std::print(stream, "\n");

  return stream.str();
}

void Framework::logStatus(HDF5Writer& hdf5, const ForceField& forceField) const
{
  std::string group = std::format("{}_{}", frameworkId, name);
  hdf5.createGroup(group);
  hdf5.writeMetaInfo(group, "Number of Atoms", unitCellAtoms.size());
  hdf5.writeMetaInfo(group, "Mass", mass);

  std::vector<std::string> typenames(definedAtoms.size());
  std::vector<double> positions(3 * definedAtoms.size());
  std::vector<double> charges(definedAtoms.size());

  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);
    typenames[i] = forceField.pseudoAtoms[atomType].name;
    positions[3 * i] = definedAtoms[i].position.x;
    positions[3 * i + 1] = definedAtoms[i].position.y;
    positions[3 * i + 2] = definedAtoms[i].position.z;
    charges[i] = definedAtoms[i].charge;
  }

  hdf5.createDataset<double>(group, "positions", {atoms.size(), 3}, {{"dimensions", "(numberOfAtoms, 3)"}});
  hdf5.writeVector(group, "positions", positions);

  hdf5.createDataset<double>(group, "charges", {atoms.size()}, {{"dimensions", "(numberOfAtoms, )"}});
  hdf5.writeVector(group, "charges", charges);

  hdf5.createStringDataset(group, "types", {atoms.size()}, 8);
  hdf5.writeVector<std::string>(group, "types", typenames);

  hdf5.writeMetaInfo(group, "Number of Bonds", bonds.size());
  std::vector<std::string> bondTypes(bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    bondTypes[i] = bonds[i].print();
  }
  hdf5.createStringDataset(group, "bondtypes", {bonds.size()}, 128);
  hdf5.writeVector<std::string>(group, "bondtypes", bondTypes);
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Framework& c)
{
  archive << c.versionNumber;

  archive << c.simulationBox;
  archive << c.spaceGroupHallNumber;
  archive << c.numberOfUnitCells;

  archive << c.frameworkId;
  archive << c.name;
  archive << c.filenameData;
  archive << c.filename;

  archive << c.rigid;

  archive << c.mass;
  archive << c.netCharge;
  archive << c.definedAtoms;
  archive << c.atoms;

  archive << c.chiralCenters;
  archive << c.bonds;
  // std::vector<std::pair<size_t, size_t>> bondDipoles{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  // std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  // std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  // std::vector<std::pair<size_t, size_t>> intraVDW{};
  // std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  // std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  // std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Framework& c)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > c.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Framework' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> c.simulationBox;
  archive >> c.spaceGroupHallNumber;
  archive >> c.numberOfUnitCells;

  archive >> c.frameworkId;
  archive >> c.name;
  archive >> c.filenameData;
  archive >> c.filename;

  archive >> c.rigid;

  archive >> c.mass;
  archive >> c.netCharge;
  archive >> c.definedAtoms;
  archive >> c.atoms;

  archive >> c.chiralCenters;
  archive >> c.bonds;
  // std::vector<std::pair<size_t, size_t>> bondDipoles{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  // std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  // std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  // std::vector<std::pair<size_t, size_t>> intraVDW{};
  // std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  // std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  // std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  return archive;
}

std::string Framework::repr() const { return std::string("Framework test"); }
