module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <exception>
#include <format>
#include <iostream>
#include <map>
#include <numbers>
#include <optional>
#include <print>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#endif

module cif_reader;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import skspacegroup;
import skelement;
import scanner;
import characterset;
import atom;
import forcefield;
import simulationbox;
import charge_equilibration_wilmer_snurr;

CIFReader::CIFReader(const std::string& content)
    :scanner(content, CharacterSet::whitespaceAndNewlineCharacterSet())
{
}

std::tuple<SimulationBox, std::size_t, std::vector<Atom>, std::vector<Atom>> 
  CIFReader::readString(const std::string& content, const ForceField& forceField, CIFReader::UseChargesFrom useChargesFrom)
{
  CIFReader cif_reader(content);

  while (!cif_reader.scanner.isAtEnd())
  {
    std::string tempString;

    // scan to first keyword
    if (cif_reader.scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString))
    {
      // FIX: cast to lower case
      std::string keyword = tempString;

      if (keyword.starts_with(std::string("_audit")))
      {
        cif_reader.parseAudit(keyword);
      }
      else if (keyword.starts_with(std::string("_chemical")))
      {
        cif_reader.parseChemical(keyword);
      }
      else if (keyword.starts_with("_cell"))
      {
        cif_reader.parseCell(keyword);
      }
      else if (keyword.starts_with(std::string("_symmetry")))
      {
        cif_reader.parseSymmetry(keyword);
      }
      else if (keyword.starts_with(std::string("_space_group")))
      {
        cif_reader.parseSymmetry(keyword);
      }
      else if (keyword.starts_with(std::string("data_")))
      {
        cif_reader.parseName(keyword);
      }
      else if (keyword.starts_with(std::string("loop_")))
      {
        cif_reader.parseLoop(keyword, forceField);
      }
      else if (keyword.starts_with("#"))
      {
        cif_reader.skipComment();
      };
    }
  }

  SimulationBox::Type type =
      (std::abs(cif_reader.alpha - 90.0) > 1.0e-3) || (std::abs(cif_reader.beta - 90.0) > 1.0e-3) || (std::abs(cif_reader.gamma - 90.0) > 1.0e-3)
          ? SimulationBox::Type::Triclinic
          : SimulationBox::Type::Rectangular;
  SimulationBox simulation_box(cif_reader.a, cif_reader.b, cif_reader.c, 
                               cif_reader.alpha * std::numbers::pi / 180.0, 
                               cif_reader.beta * std::numbers::pi / 180.0,
                               cif_reader.gamma * std::numbers::pi / 180.0, type);

  if (useChargesFrom == UseChargesFrom::PseudoAtoms)
  {
    for (Atom& atom : cif_reader.fractionalAtoms)
    {
      atom.charge = forceField.pseudoAtoms[atom.type].charge;
    }
  }

  std::vector<Atom> atoms = CIFReader::expandDefinedAtomsToUnitCell(simulation_box, cif_reader.spaceGroupHallNumber.value_or(1), cif_reader.fractionalAtoms);

  //if (useChargesFrom == UseChargesFrom::ChargeEquilibration)
  //{
  //  ChargeEquilibration::computeChargeEquilibration(forceField, simulationBox, unitCellAtoms,
  //                                                  ChargeEquilibration::Type::PeriodicEwaldSum);

  //  std::vector<std::size_t> countCharge(definedAtoms.size());
  //  std::vector<double> sumCharge(definedAtoms.size());
  //  for (const Atom& atom : unitCellAtoms)
  //  {
  //    ++countCharge[atom.moleculeId];
  //    sumCharge[atom.moleculeId] += atom.charge;
  //  }
  //  for (std::size_t i = 0; i < definedAtoms.size(); ++i)
  //  {
  //    definedAtoms[i].charge = sumCharge[i] / static_cast<double>(countCharge[i]);
  //  }
  //  for (Atom& atom : unitCellAtoms)
  //  {
  //    atom.charge = definedAtoms[atom.moleculeId].charge;
  //  }
  //}


  return {simulation_box, cif_reader.spaceGroupHallNumber.value_or(1), cif_reader.fractionalAtoms, atoms};
}

void CIFReader::parseLine([[maybe_unused]] std::string& string) {}

void CIFReader::parseAudit([[maybe_unused]] std::string& string) {}

void CIFReader::parseChemical([[maybe_unused]] std::string& string) {}

void CIFReader::parseCell(std::string& string)
{
  if (string == std::string("_cell_length_a") || string == std::string("_cell.length_a"))
  {
    a = scanDouble();
  }
  if (string == std::string("_cell_length_b") || string == std::string("_cell.length_b"))
  {
    b = scanDouble();
  }
  if (string == std::string("_cell_length_c") || string == std::string("_cell.length_c"))
  {
    c = scanDouble();
  }

  if (string == std::string("_cell_angle_alpha") || string == std::string("_cell.angle_alpha"))
  {
    alpha = scanDouble();
  }
  if (string == std::string("_cell_angle_beta") || string == std::string("_cell.angle_beta"))
  {
    beta = scanDouble();
  }
  if (string == std::string("_cell_angle_gamma") || string == std::string("_cell.angle_gamma"))
  {
    gamma = scanDouble();
  }
}

void CIFReader::parseSymmetry(std::string& string)
{
  if (string == std::string("_symmetry_cell_settings"))
  {
    return;
  }

  // prefer setting spacegroup based on Hall-symbol
  if ((string == std::string("_symmetry_space_group_name_Hall")) ||
      (string == std::string("_symmetry.space_group_name_Hall")))
  {
    std::optional<std::string> possibleString = scanString();

    if (possibleString)
    {
      spaceGroupHallNumber = SKSpaceGroup::HallNumber(*possibleString);
    }
  }

  if (!spaceGroupHallNumber)
  {
    if ((string == std::string("_space_group_name_H-M_alt")) ||
        (string == std::string("_symmetry_space_group_name_H-M")) ||
        (string == std::string("_symmetry.pdbx_full_space_group_name_H-M")))
    {
      std::optional<std::string> possibleString = scanString();
      if (possibleString)
      {
        spaceGroupHallNumber = SKSpaceGroup::HallNumberFromHMString(*possibleString);
      }
    }
  }

  if (!spaceGroupHallNumber)
  {
    if ((string == std::string("_space_group_IT_number")) || (string == std::string("_symmetry_Int_Tables_number")) ||
        (string == std::string("_symmetry.Int_Tables_number")))
    {
      std::size_t spaceGroupNumber = scanInt();
      spaceGroupHallNumber = SKSpaceGroup::HallNumberFromSpaceGroupNumber(spaceGroupNumber);
    }
  }
}

void CIFReader::parseName([[maybe_unused]] std::string& string) {}

void CIFReader::skipComment()
{
  std::string tempString;
  scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), tempString);
}

void CIFReader::parseLoop([[maybe_unused]] std::string& string, const ForceField& forceField)
{
  std::string tempString;
  std::string::const_iterator previousScanLocation;
  std::vector<std::string> tags;

  // part 1: read the 'tags'
  previousScanLocation = scanner.scanLocation();
  while (scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString) &&
         (tempString.size() > 0) &&
         (tempString.starts_with(std::string("_")) || (tempString.starts_with(std::string("#")))))
  {
    // std::string tag = tempString.toLower();
    std::string tag = tempString;

    if (tag.starts_with(std::string("_")))
    {
      tags.push_back(tag);
      previousScanLocation = scanner.scanLocation();
    }
  }

  // set scanner back to the first <value>
  scanner.setScanLocation(previousScanLocation);

  std::optional<std::string> value1 = std::nullopt;
  do
  {
    std::map<std::string, std::string> dictionary{};

    for (const std::string& tag : tags)
    {
      if ((value1 = parseValue()))
      {
        dictionary[tag] = *value1;
      }
    }

    if (value1)
    {
      if (std::map<std::string, std::string>::iterator index = dictionary.find(std::string("_atom_site_type_symbol"));
          (index != dictionary.end()))
      {
        // at least _atom_site_type_symbol is present
        std::string chemicalElement = index->second;

        if (chemicalElement.size() > 0)
        {
          std::replace_if(chemicalElement.begin(), chemicalElement.end(), [](char c) { return std::isdigit(c); }, ' ');

          // First character to uppercase
          chemicalElement[0] = static_cast<char>(std::toupper(chemicalElement[0]));
        }

        Atom atom = Atom();

        if (std::map<std::string, std::string>::iterator atomSiteIndex =
                dictionary.find(std::string("_atom_site_label"));
            (atomSiteIndex != dictionary.end()))
        {
          std::string tempString1 = atomSiteIndex->second;
          // std::replace_if(tempString1.begin(), tempString1.end(), [](char c) { return std::isdigit(c); }, ' ');
          std::istringstream ss(tempString1);
          std::string value2;
          if (ss >> value2)
          {
            std::optional<std::size_t> index1 = forceField.findPseudoAtom(value2);

            // find by stripping of numbers
            if (!index1.has_value())
            {
              std::replace_if(tempString1.begin(), tempString1.end(), [](char c) { return std::isdigit(c); }, ' ');
              std::istringstream ss2(tempString1);
              ss2 >> value2;
              index1 = forceField.findPseudoAtom(value2);
            }

            // TODO: add pseudoAtom if not found
            if (!index1.has_value())
            {
              throw std::runtime_error(std::format("[cif reader]: atom type {} not recognized\n", value2));
            }

            if (index1.has_value())
            {
              atom.type = static_cast<std::uint16_t>(index1.value());
            }
          }

          //  atom->setDisplayName(atomSiteIndex->second);
          //  atom->setUniqueForceFieldName(atomSiteIndex->second);
        }

        if (std::map<std::string, std::string>::iterator atomSiteForceFieldIndex =
                dictionary.find(std::string("_atom_site_forcefield_label"));
            (atomSiteForceFieldIndex != dictionary.end()))
        {
          //  atom->setUniqueForceFieldName(atomSiteForceFieldIndex->second);
        }

        std::map<std::string, std::string>::iterator atom_site_x = dictionary.find(std::string("_atom_site_fract_x"));
        std::map<std::string, std::string>::iterator atom_site_y = dictionary.find(std::string("_atom_site_fract_y"));
        std::map<std::string, std::string>::iterator atom_site_z = dictionary.find(std::string("_atom_site_fract_z"));
        if ((atom_site_x != dictionary.end()) && (atom_site_y != dictionary.end()) && (atom_site_z != dictionary.end()))
        {
          double3 position;
          // bool succes = false;
          position.x = scanDouble(atom_site_x->second);
          position.y = scanDouble(atom_site_y->second);
          position.z = scanDouble(atom_site_z->second);
          atom.position = position;
        }

        std::map<std::string, std::string>::iterator atom_charge = dictionary.find(std::string("_atom_site_charge"));
        if (atom_charge != dictionary.end())
        {
          double q = scanDouble(atom_charge->second);
          bool succes = true;
          // charge = atom_charge->second.split('(').at(0).toDouble(&succes);
          if (succes)
          {
            atom.charge = q;
          }
        }

        std::map<std::string, std::string>::iterator atom_occupancy =
            dictionary.find(std::string("_atom_site_occupancy"));
        if (atom_occupancy != dictionary.end())
        {
          // double occpuancy = 0.0;
          // bool succes = false;
          // occpuancy = atom_occupancy->second.split('(').at(0).toDouble(&succes);
          // if (succes)
          //{
          //   atom->setOccupancy(occpuancy);
          //}
        }

        if (std::map<std::string, std::size_t>::iterator chemicalElementIndex =
                PredefinedElements::atomicNumberData.find(chemicalElement);
            chemicalElementIndex != PredefinedElements::atomicNumberData.end())
        {
          // atom.type = chemicalElementIndex->second;
        }

        fractionalAtoms.push_back(atom);
      }
      else
      {
        if (std::map<std::string, std::string>::iterator index2 =
                dictionary.find(std::string("_atom_site.type_symbol"));
            (index2 != dictionary.end()))
        {
          Atom atom = Atom();

          if (std::map<std::string, std::string>::iterator atomSiteIndex =
                  dictionary.find(std::string("_atom_site.id"));
              (atomSiteIndex != dictionary.end()))
          {
            // atom->setDisplayName(atomSiteIndex->second);
            // atom->setUniqueForceFieldName(atomSiteIndex->second);
          }

          if (std::map<std::string, std::string>::iterator atomSiteForceFieldIndex =
                  dictionary.find(std::string("_atom_site.forcefield_label"));
              (atomSiteForceFieldIndex != dictionary.end()))
          {
            // atom->setUniqueForceFieldName(atomSiteForceFieldIndex->second);
          }

          std::map<std::string, std::string>::iterator atom_site_x = dictionary.find(std::string("_atom_site.fract_x"));
          std::map<std::string, std::string>::iterator atom_site_y = dictionary.find(std::string("_atom_site.fract_y"));
          std::map<std::string, std::string>::iterator atom_site_z = dictionary.find(std::string("_atom_site.fract_z"));
          if ((atom_site_x != dictionary.end()) && (atom_site_y != dictionary.end()) &&
              (atom_site_z != dictionary.end()))
          {
            // double3 position;
            // bool succes = false;
            // position.x = atom_site_x->second.split('(').at(0).toDouble(&succes);
            // position.y = atom_site_y->second.split('(').at(0).toDouble(&succes);
            // position.z = atom_site_z->second.split('(').at(0).toDouble(&succes);
            // atom->setPosition(position);
          }

          std::map<std::string, std::string>::iterator atom_charge = dictionary.find(std::string("_atom_site.charge"));
          if (atom_charge != dictionary.end())
          {
            // double charge = 0.0;
            // bool succes = false;
            // charge = atom_charge->second.split('(').at(0).toDouble(&succes);
            // atom->setCharge(charge);
          }

          fractionalAtoms.push_back(atom);
        }
      }
    }
  } while (value1);
}

std::optional<std::string> CIFReader::parseValue()
{
  if (scanner.isAtEnd())
  {
    return std::nullopt;
  }

  std::string::const_iterator previousScanLocation = scanner.scanLocation();

  std::string tempString;
  while (scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString) &&
         tempString.starts_with(std::string("#")))
  {
    skipComment();
  };

  // std::string keyword = tempString.toLower();
  std::string keyword = tempString;

  if ((keyword.starts_with(std::string("_")) || keyword.starts_with(std::string("loop_"))))
  {
    scanner.setScanLocation(previousScanLocation);

    return std::nullopt;
  }
  else
  {
    return tempString;
  }
}

std::size_t CIFReader::scanInt()
{
  std::string tempString;
  if (scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString))
  {
    std::replace_if(tempString.begin(), tempString.end(), [](char c) { return !std::isalnum(c); }, ' ');

    std::istringstream ss(tempString);
    std::size_t value;
    if (ss >> value)
    {
      return value;
    }
  }
  return 0;
}

double CIFReader::scanDouble()
{
  std::string tempString;
  if (scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString))
  {
    // std::replace_if(tempString.begin(), tempString.end(), [](char c) {return !std::isalnum(c); }, ' ');

    std::istringstream ss(tempString);
    double value;
    if (ss >> value)
    {
      return value;
    }
  }
  return 0.0;
}

double CIFReader::scanDouble(std::string tempString)
{
  // std::replace_if(tempString.begin(), tempString.end(), [](char c) {return !std::isdigit(c); }, ' ');

  std::istringstream ss(tempString);
  double value;
  if (ss >> value)
  {
    return value;
  }
  return 0.0;
}

std::optional<std::string> CIFReader::scanString()
{
  std::string tempString;
  if (scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), tempString))
  {
    return tempString;
  }

  return std::nullopt;
}

std::vector<Atom> CIFReader::expandDefinedAtomsToUnitCell(const SimulationBox &simulation_box, std::size_t spaceGroupHallNumber, const std::vector<Atom> &definedAtoms)
{
  SKSpaceGroup spaceGroup = SKSpaceGroup(spaceGroupHallNumber);

  // expand the fractional atoms based on the space-group
  std::vector<Atom> fractional_expanded_atoms{};
  fractional_expanded_atoms.reserve(definedAtoms.size() * 256);

  for (Atom atomCopy : definedAtoms)
  {
    std::vector<double3> listOfPositions = spaceGroup.listOfSymmetricPositions(atomCopy.position);
    for (const double3& pos : listOfPositions)
    {
      atomCopy.position = pos.fract();
      fractional_expanded_atoms.push_back(atomCopy);
    }
  }

  // eliminate duplicates
  std::vector<Atom> fractionalUnitCellAtoms{};
  for (std::size_t i = 0; i < fractional_expanded_atoms.size(); ++i)
  {
    bool overLap = false;
    double3 cartesian_position_A = simulation_box.cell * fractional_expanded_atoms[i].position;
    for (std::size_t j = i + 1; j < fractional_expanded_atoms.size(); ++j)
    {
      double3 cartesian_position_B = simulation_box.cell * fractional_expanded_atoms[j].position;
      double3 dr = cartesian_position_A - cartesian_position_B;
      dr = simulation_box.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      if (rr < 0.1)
      {
        overLap = true;
        break;
      }
    }
    if (!overLap)
    {
      fractionalUnitCellAtoms.push_back(fractional_expanded_atoms[i]);
    }
  }

  return fractionalUnitCellAtoms;
}

