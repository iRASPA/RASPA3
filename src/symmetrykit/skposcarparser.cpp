module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <iostream>
#include <locale>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#endif

module skposcarparser;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import skelement;
import skatom;
import double3;
import double3x3;
import scanner;
import characterset;
import skparser;
import skstructure;
import skcell;

static std::string& tolower(std::string& s)
{
  for (auto& c : s)
  {
    [[maybe_unused]] auto t = std::tolower(static_cast<unsigned char>(c));
  }

  return s;
}

static std::vector<std::string> split(const std::string txt, char ch)
{
  std::size_t pos = txt.find(ch);
  std::size_t initialPos = 0;
  std::vector<std::string> strs{};

  // Decompose statement
  while (pos != std::string::npos)
  {
    std::string s = txt.substr(initialPos, pos - initialPos);
    if (!s.empty())
    {
      strs.push_back(s);
    }
    initialPos = pos + 1;

    pos = txt.find(ch, initialPos);
  }

  // Add the last one
  std::string s = txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1);
  if (!s.empty())
  {
    strs.push_back(s);
  }

  return strs;
}

inline bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
  return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(),
                                                  [](auto a, auto b) { return std::tolower(a) == std::tolower(b); });
}

inline static std::string trim2(const std::string& s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start))
  {
    start++;
  }

  auto end = s.end();
  do
  {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}

SKPOSCARParser::SKPOSCARParser(const std::string& content, bool proteinOnlyAsymmetricUnitCell, bool asMolecule,
                               CharacterSet charactersToBeSkipped)
    : SKParser(),
      _scanner(content, charactersToBeSkipped),
      _proteinOnlyAsymmetricUnitCell(proteinOnlyAsymmetricUnitCell),
      _asMolecule(asMolecule),
      _frame(std::make_shared<SKStructure>()),
      _spaceGroupHallNumber(1)
{
  _frame->kind = SKStructure::Kind::crystal;
}

void SKPOSCARParser::startParsing() noexcept(false)
{
  double3x3 unitCell{};
  double3x3 inverseUnitCell{};

  std::string scannedLine;
  std::vector<std::string> termsScannedLined{};

  // skip first line
  _scanner.scanLine(scannedLine);
  if (scannedLine.empty())
  {
    throw std::runtime_error("Empty file");
  }

  // skip first line
  _scanner.scanLine(scannedLine);
  if (scannedLine.empty())
  {
    throw std::runtime_error("Empty file");
  }

  // read first lattice vector
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

  termsScannedLined = split(scannedLine, ' ');

  if (termsScannedLined.size() < 3)
  {
    throw std::runtime_error("Missing first lattice vector in POSCAR");
  }

  unitCell.ax = std::stod(termsScannedLined[0]);
  unitCell.ay = std::stod(termsScannedLined[1]);
  unitCell.az = std::stod(termsScannedLined[2]);

  // read second lattice vector
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

  termsScannedLined = split(scannedLine, ' ');

  if (termsScannedLined.size() < 3)
  {
    throw std::runtime_error("Missing second lattice vector in POSCAR");
  }

  unitCell.bx = std::stod(termsScannedLined[0]);
  unitCell.by = std::stod(termsScannedLined[1]);
  unitCell.bz = std::stod(termsScannedLined[2]);

  // read third lattice vector
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

  termsScannedLined = split(scannedLine, ' ');

  if (termsScannedLined.size() < 3)
  {
    throw std::runtime_error("Missing third lattice vector in POSCAR");
  }

  unitCell.cx = std::stod(termsScannedLined[0]);
  unitCell.cy = std::stod(termsScannedLined[1]);
  unitCell.cz = std::stod(termsScannedLined[2]);

  _frame->cell = std::make_shared<SKCell>(unitCell);
  inverseUnitCell = unitCell.inverse();

  // read amount of atoms per element
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

  std::vector<std::string> elementList = split(scannedLine, ' ');

  if (elementList.empty())
  {
    throw std::runtime_error("List of types of atoms is empty");
  }

  // read amount of atoms per element
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

  std::vector<std::string> amountList = split(scannedLine, ' ');

  if (amountList.empty())
  {
    throw std::runtime_error("List of amount of atoms is empty");
  }

  // skip first line
  _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

  std::vector<std::string> directOrCartesian = split(scannedLine, ' ');

  bool cartesian = false;
  if (tolower(directOrCartesian[0]) == "cartesian")
  {
    cartesian = true;
  }

  for (std::size_t k = 0; k < amountList.size(); k++)
  {
    int numberOfAtoms = std::stoi(amountList[k]);
    std::string chemicalElementString = trim2(elementList[k]);

    for (int i = 0; i < numberOfAtoms; i++)
    {
      std::shared_ptr<SKAsymmetricAtom> atom = std::make_shared<SKAsymmetricAtom>();

      // read atom
      _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

      termsScannedLined = split(scannedLine, ' ');

      if (termsScannedLined.empty())
      {
        throw std::runtime_error("Error reading atoms");
      }

      atom->setDisplayName(chemicalElementString);
      if (std::map<std::string, std::size_t>::iterator index =
              PredefinedElements::atomicNumberData.find(chemicalElementString);
          index != PredefinedElements::atomicNumberData.end())
      {
        atom->setElementIdentifier(index->second);
      }

      double3 position;
      position.x = std::stod(termsScannedLined[0]);
      position.y = std::stod(termsScannedLined[1]);
      position.z = std::stod(termsScannedLined[2]);

      // convert to fractional position if Cartesian coordinates are specified
      if (cartesian)
      {
        atom->setPosition(inverseUnitCell * position);
      }
      else
      {
        atom->setPosition(position);
      }

      _frame->atoms.push_back(atom);
    }
  }

  _movies.push_back({_frame});
}
