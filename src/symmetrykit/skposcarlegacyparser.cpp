module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <locale>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#endif

module skposcarlegacyparser;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <locale>;
import <optional>;
#endif

import skstructure;
import double3;
import double3x3;
import skatom;
import characterset;
import skcell;

static std::string& tolower(std::string& s)
{
  for (auto& c : s)
  {
    [[maybe_unused]] int result = std::tolower(static_cast<unsigned char>(c));
  }

  return s;
}

static std::vector<std::string> split(const std::string txt, char ch)
{
  size_t pos = txt.find(ch);
  size_t initialPos = 0;
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

SKPOSCARLegacyParser::SKPOSCARLegacyParser(const std::string& content, bool proteinOnlyAsymmetricUnitCell,
                                           bool asMolecule, CharacterSet charactersToBeSkipped)
    : SKParser(),
      _scanner(content, charactersToBeSkipped),
      _proteinOnlyAsymmetricUnitCell(proteinOnlyAsymmetricUnitCell),
      _asMolecule(asMolecule),
      _frame(std::make_shared<SKStructure>()),
      _spaceGroupHallNumber(1)
{
  _frame->kind = SKStructure::Kind::crystal;
}

void SKPOSCARLegacyParser::startParsing() noexcept(false)
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

  for (size_t k = 0; k < amountList.size(); k++)
  {
    int numberOfAtoms = std::stoi(amountList[k]);

    for (size_t i = 0; i < static_cast<size_t>(numberOfAtoms); i++)
    {
      std::shared_ptr<SKAsymmetricAtom> atom = std::make_shared<SKAsymmetricAtom>();

      // read atom
      _scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), scannedLine);

      termsScannedLined = split(scannedLine, ' ');

      if (termsScannedLined.empty())
      {
        throw std::runtime_error("Error reading atoms");
      }

      atom->setElementIdentifier(k + 1);
      atom->setDisplayName(std::to_string(k + 1));

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
