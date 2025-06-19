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
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include <numbers>
#include <limits>
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
import <print>;
import <limits>;
#endif

import archive;
import int3;
import double3;
import double2;
import double3x3;
import skposcarparser;
import skelement;
import characterset;
import stringutils;
import skparser;
import skposcarparser;
import skstructure;
import skatom;
import skcell;
import skspacegroup;
import forcefield;
import atom;
import multi_site_isotherm;
import simulationbox;
import cif_reader;
import move_statistics;
import bond_potential;
import json;
import charge_equilibration_wilmer_snurr;

// default constructor, needed for binary restart-file
Framework::Framework() {}

// create Component in 'inputreader.cpp'
Framework::Framework(size_t currentFramework, const ForceField& forceField, const std::string& componentName,
                     std::optional<const std::string> fileName, std::optional<int3> numberOfUnitCells,
                     Framework::UseChargesFrom useChargesFrom) noexcept(false)
    : frameworkId(currentFramework),
      name(componentName),
      filenameData(fileName),
      useChargesFrom(useChargesFrom)
{
  if (filenameData.has_value())
  {
    readFramework(forceField, filenameData.value());

    double cutOff = forceField.cutOffFrameworkVDW;
    const int3 minimum_number_of_unit_cells = simulationBox.smallestNumberOfUnitCellsForMinimumImagesConvention(cutOff);

    this->numberOfUnitCells = numberOfUnitCells.value_or(minimum_number_of_unit_cells);

    makeSuperCell();

    unitCellMass = 0.0;
    for (const Atom& atom : unitCellAtoms)
    {
      size_t atomType = static_cast<size_t>(atom.type);
      unitCellMass += forceField.pseudoAtoms[atomType].mass;
    }

    mass = 0.0;
    netCharge = 0.0;
    smallestCharge = std::numeric_limits<double>::max();
    largestCharge = std::numeric_limits<double>::lowest();
    for (const Atom& atom : atoms)
    {
      size_t atomType = static_cast<size_t>(atom.type);
      mass += forceField.pseudoAtoms[atomType].mass;
      netCharge += atom.charge;
      if (atom.charge > largestCharge) largestCharge = atom.charge;
      if (atom.charge < smallestCharge) smallestCharge = atom.charge;
    }

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
      useChargesFrom(UseChargesFrom::PseudoAtoms),
      definedAtoms(definedAtoms)
{
  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
    definedAtoms[i].moleculeId = static_cast<uint32_t>(i);
  }

  if (useChargesFrom == UseChargesFrom::PseudoAtoms)
  {
    for (Atom& atom : definedAtoms)
    {
      atom.charge = forceField.pseudoAtoms[atom.type].charge;
    }
  }

  expandDefinedAtomsToUnitCell();

  if (useChargesFrom == UseChargesFrom::ChargeEquilibration)
  {
    ChargeEquilibration::computeChargeEquilibration(forceField, simulationBox, unitCellAtoms,
                                                    ChargeEquilibration::Type::PeriodicEwaldSum);

    std::vector<size_t> countCharge(definedAtoms.size());
    std::vector<double> sumCharge(definedAtoms.size());
    for (const Atom& atom : unitCellAtoms)
    {
      ++countCharge[atom.moleculeId];
      sumCharge[atom.moleculeId] += atom.charge;
    }
    for (size_t i = 0; i < definedAtoms.size(); ++i)
    {
      definedAtoms[i].charge = sumCharge[i] / static_cast<double>(countCharge[i]);
    }
    for (Atom& atom : unitCellAtoms)
    {
      atom.charge = definedAtoms[atom.moleculeId].charge;
    }
  }

  makeSuperCell();

  unitCellMass = 0.0;
  for (const Atom& atom : unitCellAtoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    unitCellMass += forceField.pseudoAtoms[atomType].mass;
  }

  mass = 0.0;
  netCharge = 0.0;
  smallestCharge = std::numeric_limits<double>::max();
  largestCharge = std::numeric_limits<double>::lowest();
  for (const Atom& atom : atoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    mass += forceField.pseudoAtoms[atomType].mass;
    netCharge += atom.charge;

    if (atom.charge > largestCharge) largestCharge = atom.charge;
    if (atom.charge < smallestCharge) smallestCharge = atom.charge;
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

  std::string frameworkFileName = addExtension(fileName, ".cif");

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

  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
    definedAtoms[i].moleculeId = static_cast<uint32_t>(i);
  }

  if (useChargesFrom == UseChargesFrom::PseudoAtoms)
  {
    for (Atom& atom : definedAtoms)
    {
      atom.charge = forceField.pseudoAtoms[atom.type].charge;
    }
  }

  // expand the fractional atoms based on the space-group
  expandDefinedAtomsToUnitCell();

  if (useChargesFrom == UseChargesFrom::ChargeEquilibration)
  {
    ChargeEquilibration::computeChargeEquilibration(forceField, simulationBox, unitCellAtoms,
                                                    ChargeEquilibration::Type::PeriodicEwaldSum);

    std::vector<size_t> countCharge(definedAtoms.size());
    std::vector<double> sumCharge(definedAtoms.size());
    for (const Atom& atom : unitCellAtoms)
    {
      ++countCharge[atom.moleculeId];
      sumCharge[atom.moleculeId] += atom.charge;
    }
    for (size_t i = 0; i < definedAtoms.size(); ++i)
    {
      definedAtoms[i].charge = sumCharge[i] / static_cast<double>(countCharge[i]);
    }
    for (Atom& atom : unitCellAtoms)
    {
      atom.charge = definedAtoms[atom.moleculeId].charge;
    }
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

std::vector<Atom> Framework::makeSuperCell(int3 numberOfCells) const
{
  std::vector<Atom> superCellAtoms{};

  for (int32_t i = 0; i < numberOfCells.x; ++i)
  {
    for (int32_t j = 0; j < numberOfCells.y; ++j)
    {
      for (int32_t k = 0; k < numberOfCells.z; ++k)
      {
        for (const Atom& atom : unitCellAtoms)
        {
          Atom atomCopy = atom;
          atomCopy.position +=
              simulationBox.cell * double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          superCellAtoms.push_back(atomCopy);
        }
      }
    }
  }

  return superCellAtoms;
}

std::optional<double> Framework::computeLargestNonOverlappingFreeRadius(const ForceField &forceField, double3 probe_position, double well_depth_factor) const
{
  double smallest_radius = std::numeric_limits<double>::max();


  // if inside blockingpocket, then return
  //    return std::nullopt;

  for(const Atom& atom: unitCellAtoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    double size_parameter = forceField[atomType].sizeParameter();
    double3 dr = probe_position - atom.position;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    double mixing_radius = 0.5 * well_depth_factor * size_parameter;
    if(rr < mixing_radius * mixing_radius)
    {
      return std::nullopt;
    }

    double radius = std::sqrt(rr) - mixing_radius;
    smallest_radius = std::min(smallest_radius, radius);
  }

  return smallest_radius;
}

bool Framework::computeVanDerWaalsRadiusOverlap(const ForceField &forceField, double3 probe_position) const
{
  for(const Atom& atom: unitCellAtoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    size_t atomicNumber = forceField.pseudoAtoms[atomType].atomicNumber;
    double radius = PredefinedElements::predefinedElements[atomicNumber]._VDWRadius;
    double3 dr = probe_position - atom.position;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if(rr < radius * radius)
    {
      return true;
    }
  }

  return false;
}
bool Framework::computeOverlap(const ForceField &forceField, double3 probe_position, double well_depth_factor, size_t probe_type, std::make_signed_t<std::size_t> skip) const
{
  for(std::make_signed_t<std::size_t> atom_index = 0; const Atom& atom: unitCellAtoms)
  {
    if(atom_index != skip)
    {
      size_t atomType = static_cast<size_t>(atom.type);
      double size_parameter = forceField(probe_type, atomType).sizeParameter();
      double3 dr = probe_position - atom.position;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);

      double radius = well_depth_factor * size_parameter;
      if(rr < radius * radius)
      {
        return true;
      }
    }

    ++atom_index;
  }

  return false;
}

std::vector<double3> Framework::fractionalAtomPositionsUnitCell() const
{
  std::vector<double3> positions;
  positions.reserve(unitCellAtoms.size());

  double3x3 inverseCell = simulationBox.inverseCell;
  for(const Atom& atom: unitCellAtoms)
  {
    double3 s = (inverseCell * atom.position).fract();
    positions.push_back(s);
  }

  return positions;
}

std::vector<double2> Framework::atomUnitCellLennardJonesPotentialParameters(const ForceField& forceField) const
{
  std::vector<double2> parameters;
  parameters.reserve(unitCellAtoms.size());

  for(const Atom& atom: unitCellAtoms)
  {
    size_t type = atom.type; 
    double2 parameter = double2(forceField[type].parameters.x, forceField[type].parameters.y);
    parameters.push_back(parameter);
  }


  return parameters;
}

std::string Framework::printStatus(const ForceField& forceField) const
{
  std::ostringstream stream;

  std::print(stream, "Framework {} [{}]\n\n", frameworkId, name);

  std::print(stream, "    Box:     {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.cell.ax, simulationBox.cell.bx, simulationBox.cell.cx);
  std::print(stream, "             {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.cell.ay, simulationBox.cell.by, simulationBox.cell.cy);
  std::print(stream, "             {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.cell.az, simulationBox.cell.bz, simulationBox.cell.cz);
  std::print(stream, "    Lengths: {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.lengthA, simulationBox.lengthB, simulationBox.lengthC);
  double conv = 180.0 / std::numbers::pi;
  std::print(stream, "    Angles:  {:9.5f} {:9.5f} {:9.5f}\n", conv * simulationBox.angleAlpha, conv * simulationBox.angleBeta, conv * simulationBox.angleGamma);
  double3 widths = simulationBox.perpendicularWidths();
  std::print(stream, "    Perpendicular widths:  {:9.5f} {:9.5f} {:9.5f}\n\n", widths.x, widths.y, widths.z);

  std::print(stream, "    number Of Atoms:          {:>12d} [-]\n", unitCellAtoms.size());
  std::print(stream, "    mass:                     {:>12.5f} [amu]\n", mass);

  for (size_t i = 0; i != definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);

    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomTypeString,
               definedAtoms[i].position.x, definedAtoms[i].position.y, definedAtoms[i].position.z,
               definedAtoms[i].charge);
  }
  switch (useChargesFrom)
  {
    case UseChargesFrom::PseudoAtoms:
      std::print(stream, "    use charge from: pseudo-atoms\n");
      break;
    case UseChargesFrom::CIF_File:
      std::print(stream, "    use charge from: CIF-file\n");
      break;
    case UseChargesFrom::ChargeEquilibration:
      std::print(stream, "    use charge from: charge-equilibration method\n");
      break;
  }
  std::print(stream, "    net charge:                 {:>12.5f} [e]\n", netCharge);
  std::print(stream, "    smallest charge:            {:>12.5f} [e]\n", smallestCharge);
  std::print(stream, "    largest charge:             {:>12.5f} [e]\n", largestCharge);

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    std::print(stream, "        {}", bonds[i].print());
  }
  std::print(stream, "\n");

  return stream.str();
}

nlohmann::json Framework::jsonStatus() const
{
  nlohmann::json status;
  status["name"] = name;
  status["id"] = frameworkId;
  status["mass"] = mass;

  // TODO I feel that the masses, positions and charges belong in the hdf5.
  status["n_bonds"] = bonds.size();

  std::vector<std::string> bondTypes(bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    bondTypes[i] = bonds[i].print();
  }
  status["bondTypes"] = bondTypes;

  return status;
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

  archive << c.useChargesFrom;
  archive << c.netCharge;
  archive << c.smallestCharge;
  archive << c.largestCharge;

  archive << c.definedAtoms;
  archive << c.unitCellAtoms;
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

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

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

  archive >> c.useChargesFrom;
  archive >> c.netCharge;
  archive >> c.smallestCharge;
  archive >> c.largestCharge;

  archive >> c.definedAtoms;
  archive >> c.unitCellAtoms;
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

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Framework: Error in binary restart\n"));
  }
#endif

  return archive;
}

std::string Framework::repr() const
{
  std::ostringstream stream;

  std::print(stream, "Framework {} [{}]\n\n", frameworkId, name);

  std::print(stream, "    number Of Atoms:  {}\n", unitCellAtoms.size());
  std::print(stream, "    net charge:       {:12.5f} [e]\n", netCharge);
  std::print(stream, "    mass:             {:12.5f} [amu]\n", mass);

  for (size_t i = 0; i != definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);

    std::print(stream, "    {:3d}: type {:3d} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomType,
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
