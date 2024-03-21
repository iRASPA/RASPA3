module;

#if !defined(_WIN32)
#include <assert.h>
#endif

module framework;

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
import <print>;
#if defined(_WIN32)
import <cassert>;
#endif
import <exception>;
import <source_location>;
import <complex>;

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

// default constructor, needed for binary restart-file
Framework::Framework()
{
}

// create Component in 'inputreader.cpp'
Framework::Framework(size_t currentComponent, const ForceField& forceField, const std::string &componentName,
                     std::optional<const std::string> fileName) noexcept(false) :
                     frameworkId(currentComponent),
                     name(componentName),
                     filenameData(fileName)
{
  if (filenameData.has_value())
  {
    readFramework(forceField, filenameData.value());

    for (size_t i = 0; i < unitCellAtoms.size(); ++i)
    {
      unitCellAtoms[i].componentId = static_cast<uint8_t>(currentComponent);
      unitCellAtoms[i].moleculeId = 0;
    }
  }
}

// create programmatically an 'framework' component
Framework::Framework(size_t componentId, std::string fileName, double mass, SimulationBox simulationBox, 
                     size_t spaceGroupHallNumber, std::vector<Atom> definedAtoms, int3 numberOfUnitCells) noexcept(false) :
    simulationBox(simulationBox),
    spaceGroupHallNumber(spaceGroupHallNumber),
    numberOfUnitCells(numberOfUnitCells),
    frameworkId(componentId),
    name(fileName),
    filenameData(fileName),
    mass(mass),
    definedAtoms(definedAtoms)
{
  // expand the fractional atoms based on the space-group
  SKSpaceGroup spaceGroup = SKSpaceGroup(spaceGroupHallNumber);
  std::vector<Atom> expandedAtoms;
  expandedAtoms.reserve(definedAtoms.size() * 256uz);
  for (const Atom& atom : definedAtoms)
  {
    Atom atomCopy = atom;

    std::vector<double3> listOfPositions = spaceGroup.listOfSymmetricPositions(atom.position);
    for (const double3& pos : listOfPositions)
    {
      atomCopy.position = simulationBox.cell * pos.fract();
      expandedAtoms.push_back(atomCopy);
    }
  }

  // eliminate duplicates
  unitCellAtoms.clear();
  for (size_t i = 0; i < expandedAtoms.size(); ++i)
  {
    bool overLap = false;
    for (size_t j = i + 1; j < expandedAtoms.size(); ++j)
    {
      double3 dr = expandedAtoms[i].position - expandedAtoms[j].position;
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
      unitCellAtoms.push_back(expandedAtoms[i]);
    }
  }

  //for (int32_t i = 0; i < numberOfUnitCells.x; ++i)
  //{
  //  for (int32_t j = 0; j < numberOfUnitCells.y; ++j)
  //  {
  //    for (int32_t k = 0; k < numberOfUnitCells.z; ++k)
  //    {
  //      for (const Atom& atom : unitCellAtoms)
  //      {
  //        Atom atomCopy = atom;
  //        atomCopy.position += simulationBox.unitCell * 
  //                             double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
  //        atoms.push_back(atomCopy);
  //      }
  //    }
  //  }
  //}


  for (size_t i = 0; i < unitCellAtoms.size(); ++i)
  {
    unitCellAtoms[i].componentId = static_cast<uint8_t>(componentId);
    unitCellAtoms[i].moleculeId = 0;
  }
}

std::vector<Atom> Framework::frameworkAtoms() const
{
  std::vector<Atom> atoms1{};

  for (int32_t i = 0; i < numberOfUnitCells.x; ++i)
  {
    for (int32_t j = 0; j < numberOfUnitCells.y; ++j)
    {
      for (int32_t k = 0; k < numberOfUnitCells.z; ++k)
      {
        for (const Atom& atom : unitCellAtoms)
        {
          Atom atomCopy = atom;
          atomCopy.position += simulationBox.cell *
                               double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          atoms1.push_back(atomCopy);
        }
      }
    }
  }
  return atoms1;
}


template<typename T>
std::vector<T> parseListOfParameters(const std::string& arguments, size_t lineNumber)
{
  std::vector<T> list{};

  std::string str;
  std::istringstream ss(arguments);

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read at line: {}\n", lineNumber));
      }
      return list;
    }
    T value;
    std::istringstream s(str);
    if (s >> value)
    {
      list.push_back(value);
    }
    else
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read at line: {}\n", lineNumber));
      }
      return list;
    }
  };

  return list;
}


void Framework::readFramework(const ForceField& forceField, const std::string& fileName)
{
  const char* env_p = std::getenv("RASPA_DIR");

  const std::string frameworkFileName = fileName + ".cif";

  std::filesystem::path frameworkPathfile = std::filesystem::path(frameworkFileName);
  if (!std::filesystem::exists(frameworkPathfile)) 
    frameworkPathfile = std::filesystem::path(env_p) / frameworkFileName;

  if (!std::filesystem::exists(frameworkPathfile)) 
    throw std::runtime_error(std::format("File '{}' not found\n", frameworkFileName));

  std::ifstream t(frameworkPathfile);
  std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
 
  CIFReader parser = CIFReader(fileContent, forceField);
  simulationBox = parser.simulationBox;
  definedAtoms = parser.fractionalAtoms;

  // expand the fractional atoms based on the space-group
  size_t number = parser._spaceGroupHallNumber.value_or(1);
  SKSpaceGroup spaceGroup = SKSpaceGroup(number);
  std::vector<Atom> expandedAtoms;
  expandedAtoms.reserve(parser.fractionalAtoms.size() * 256);
  for (const Atom &atom : definedAtoms)
  {
    Atom atomCopy = atom;
    std::vector<double3> listOfPositions = spaceGroup.listOfSymmetricPositions(atom.position);
    for (const double3& pos : listOfPositions)
    {
      atomCopy.position = parser.simulationBox.cell * pos.fract();
      expandedAtoms.push_back(atomCopy);
    }
  }

  // eliminate duplicates
  unitCellAtoms.clear();
  for (size_t i = 0; i < expandedAtoms.size(); ++i)
  {
    bool overLap = false;
    for (size_t j = i + 1; j < expandedAtoms.size(); ++j)
    {
      double3 dr = expandedAtoms[i].position - expandedAtoms[j].position;
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
      unitCellAtoms.push_back(expandedAtoms[i]);
    }
  }

  //for (int32_t i = 0; i < numberOfUnitCells.x; ++i)
  //{
  //  for (int32_t j = 0; j < numberOfUnitCells.y; ++j)
  //  {
  //    for (int32_t k = 0; k < numberOfUnitCells.z; ++k)
  //    {
  //      for (const Atom& atom : unitCellAtoms)
  //      {
  //        Atom atomCopy = atom;
  //        atomCopy.position += simulationBox.unitCell * 
  //          double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
  //        atoms.push_back(atomCopy);
  //      }
  //    }
  //  }
  //}

  mass = 0.0;
  for (const Atom& atom : unitCellAtoms)
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

std::string Framework::printStatus(const ForceField& forceField) const
{
  std::ostringstream stream;

  std::print(stream, "Framework {} [{}]\n\n", frameworkId, name);

  std::print(stream, "    number Of Atoms:  {}\n", unitCellAtoms.size());
  std::print(stream, "    Mass:             {} [amu]\n", mass);
  
  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);
    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", 
                   i, atomTypeString, definedAtoms[i].position.x, definedAtoms[i].position.y, 
                   definedAtoms[i].position.z, definedAtoms[i].charge);
  }

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    std::print(stream, "        {}", bonds[i].print());
  }
  std::print(stream, "\n");

  return stream.str();
}


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Framework &c)
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
  //std::vector<std::pair<size_t, size_t>> bondDipoles{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  //std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  //std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  //std::vector<std::pair<size_t, size_t>> intraVDW{};
  //std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  //std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  //std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Framework &c)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > c.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Framework' at line {} in file {}\n",
                                         location.line(), location.file_name()));
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
  //std::vector<std::pair<size_t, size_t>> bondDipoles{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  //std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  //std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  //std::vector<std::pair<size_t, size_t>> intraVDW{};
  //std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  //std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  //std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  return archive;
}

std::string Framework::repr() const
{
  return std::string("Framework test");
}
