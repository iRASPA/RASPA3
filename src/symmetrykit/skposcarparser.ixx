module;

export module skposcarparser;

import std;

import scanner;
import characterset;

import skparser;
import skstructure;
import skcell;
import skatom;

export class SKPOSCARParser : public SKParser
{
 public:
  SKPOSCARParser(const std::string& content, bool onlyAsymmetricUnitCell = false, bool asMolecule = false,
                 CharacterSet charactersToBeSkipped = CharacterSet::whitespaceAndNewlineCharacterSet());
  void startParsing() noexcept(false) override final;

 private:
  Scanner _scanner;
  [[maybe_unused]] bool _proteinOnlyAsymmetricUnitCell;
  [[maybe_unused]] bool _asMolecule;
  [[maybe_unused]] std::string::const_iterator _previousScanLocation;

  [[maybe_unused]] int _numberOfAtoms = 0;
  [[maybe_unused]] int _numberOfAminoAcidAtoms = 0;
  std::shared_ptr<SKStructure> _frame;
  std::optional<SKCell> _cell;
  [[maybe_unused]] int _spaceGroupHallNumber;
};
