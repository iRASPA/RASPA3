module;

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
#endif

module cif_reader;

#ifndef USE_LEGACY_HEADERS
import <map>;
import <vector>;
import <string>;
import <optional>;
import <algorithm>;
import <sstream>;
import <cmath>;
import <cctype>;
import <numbers>;
import <iostream>;
import <exception>;
import <format>;
import <print>;
#endif

import double3;
import skspacegroup;
import skelement;
import scanner;
import characterset;
import atom;
import forcefield;
import simulationbox;

CIFReader::CIFReader(const std::string& content, const ForceField& forceField)
    : _scanner(content, CharacterSet::whitespaceAndNewlineCharacterSet())
{
  while (!_scanner.isAtEnd())
  {
    std::string tempString;

    // scan to first keyword
    _previousScanLocation = _scanner.scanLocation();
    if (_scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString))
    {
      // FIX: cast to lower case
      std::string keyword = tempString;

      if (keyword.starts_with(std::string("_audit")))
      {
        parseAudit(keyword);
      }
      else if (keyword.starts_with(std::string("_chemical")))
      {
        parseChemical(keyword);
      }
      else if (keyword.starts_with("_cell"))
      {
        parseCell(keyword);
      }
      else if (keyword.starts_with(std::string("_symmetry")))
      {
        parseSymmetry(keyword);
      }
      else if (keyword.starts_with(std::string("_symmetry_space_group")))
      {
        parseSymmetry(keyword);
      }
      else if (keyword.starts_with(std::string("data_")))
      {
        parseName(keyword);
      }
      else if (keyword.starts_with(std::string("loop_")))
      {
        parseLoop(keyword, forceField);
      }
      else if (keyword.starts_with("#"))
      {
        skipComment();
      };
    }
  }

  SimulationBox::Type type =
      (std::abs(_alpha - 90.0) > 1.0e-3) || (std::abs(_beta - 90.0) > 1.0e-3) || (std::abs(_gamma - 90.0) > 1.0e-3)
          ? SimulationBox::Type::Triclinic
          : SimulationBox::Type::Rectangular;
  simulationBox = SimulationBox(_a, _b, _c, _alpha * std::numbers::pi / 180.0, _beta * std::numbers::pi / 180.0,
                                _gamma * std::numbers::pi / 180.0, type);
}

void CIFReader::parseLine([[maybe_unused]] std::string& string) {}

void CIFReader::parseAudit([[maybe_unused]] std::string& string) {}

void CIFReader::parseChemical([[maybe_unused]] std::string& string) {}

void CIFReader::parseCell(std::string& string)
{
  if (string == std::string("_cell_length_a") || string == std::string("_cell.length_a"))
  {
    _a = scanDouble();
  }
  if (string == std::string("_cell_length_b") || string == std::string("_cell.length_b"))
  {
    _b = scanDouble();
  }
  if (string == std::string("_cell_length_c") || string == std::string("_cell.length_c"))
  {
    _c = scanDouble();
  }

  if (string == std::string("_cell_angle_alpha") || string == std::string("_cell.angle_alpha"))
  {
    _alpha = scanDouble();
  }
  if (string == std::string("_cell_angle_beta") || string == std::string("_cell.angle_beta"))
  {
    _beta = scanDouble();
  }
  if (string == std::string("_cell_angle_gamma") || string == std::string("_cell.angle_gamma"))
  {
    _gamma = scanDouble();
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
      _spaceGroupHallNumber = SKSpaceGroup::HallNumber(*possibleString);
    }
  }

  if (!_spaceGroupHallNumber)
  {
    if ((string == std::string("_space_group_name_H-M_alt")) ||
        (string == std::string("_symmetry_space_group_name_H-M")) ||
        (string == std::string("_symmetry.pdbx_full_space_group_name_H-M")))
    {
      std::optional<std::string> possibleString = scanString();
      if (possibleString)
      {
        _spaceGroupHallNumber = SKSpaceGroup::HallNumberFromHMString(*possibleString);
      }
    }
  }

  if (!_spaceGroupHallNumber)
  {
    if ((string == std::string("_space_group_IT_number")) || (string == std::string("_symmetry_Int_Tables_number")) ||
        (string == std::string("_symmetry.Int_Tables_number")))
    {
      size_t spaceGroupNumber = scanInt();
      _spaceGroupHallNumber = SKSpaceGroup::HallNumberFromSpaceGroupNumber(spaceGroupNumber);
    }
  }
}

void CIFReader::parseName([[maybe_unused]] std::string& string) {}

void CIFReader::skipComment()
{
  std::string tempString;
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), tempString);
}

void CIFReader::parseLoop([[maybe_unused]] std::string& string, const ForceField& forceField)
{
  std::string tempString;
  std::string::const_iterator previousScanLocation;
  std::vector<std::string> tags;

  // part 1: read the 'tags'
  previousScanLocation = _scanner.scanLocation();
  while (_scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString) &&
         (tempString.size() > 0) &&
         (tempString.starts_with(std::string("_")) || (tempString.starts_with(std::string("#")))))
  {
    // std::string tag = tempString.toLower();
    std::string tag = tempString;

    if (tag.starts_with(std::string("_")))
    {
      tags.push_back(tag);
      previousScanLocation = _scanner.scanLocation();
    }
  }

  // set scanner back to the first <value>
  _scanner.setScanLocation(previousScanLocation);

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
          chemicalElement[0] = static_cast<char>(toupper(chemicalElement[0]));
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
            std::optional<size_t> index1 = forceField.findPseudoAtom(value2);

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
              atom.type = static_cast<uint16_t>(index1.value());
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
          double c = scanDouble(atom_charge->second);
          bool succes = true;
          // charge = atom_charge->second.split('(').at(0).toDouble(&succes);
          if (succes)
          {
            atom.charge = c;
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

        if (std::map<std::string, size_t>::iterator chemicalElementIndex =
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
  if (_scanner.isAtEnd())
  {
    return std::nullopt;
  }

  std::string::const_iterator previousScanLocation = _scanner.scanLocation();

  std::string tempString;
  while (_scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString) &&
         tempString.starts_with(std::string("#")))
  {
    skipComment();
  };

  // std::string keyword = tempString.toLower();
  std::string keyword = tempString;

  if ((keyword.starts_with(std::string("_")) || keyword.starts_with(std::string("loop_"))))
  {
    _scanner.setScanLocation(previousScanLocation);

    return std::nullopt;
  }
  else
  {
    return tempString;
  }
}

size_t CIFReader::scanInt()
{
  std::string tempString;
  if (_scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString))
  {
    std::replace_if(tempString.begin(), tempString.end(), [](char c) { return !std::isalnum(c); }, ' ');

    std::istringstream ss(tempString);
    size_t value;
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
  if (_scanner.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), tempString))
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
  if (_scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), tempString))
  {
    return tempString;
  }

  return std::nullopt;
}
