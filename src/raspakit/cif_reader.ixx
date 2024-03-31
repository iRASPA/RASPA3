export module cif_reader;

import <string>;
import <optional>;
import <vector>;

import scanner;
import characterset;
import atom;
import simulationbox;
import forcefield;

export struct CIFReader
{
  CIFReader(const std::string& content, const ForceField& forceField);
 
  void parseLine(std::string &string);
  void parseAudit(std::string& string);
  void parseChemical(std::string& string);
  void parseCell(std::string& string);
  void parseSymmetry(std::string& string);
  void parseName(std::string& string);
  void parseLoop([[maybe_unused]] std::string& string, const ForceField& forceField);
  void skipComment();

  std::optional<std::string> parseValue();
  size_t scanInt();
  double scanDouble();
  double scanDouble(std::string tempString);
  std::optional<std::string> scanString();

  Scanner _scanner;
  std::string::const_iterator _previousScanLocation;

  std::vector<Atom> fractionalAtoms;
  SimulationBox simulationBox;
  std::optional<size_t> _spaceGroupHallNumber;
  double _a;
  double _b;
  double _c;
  double _alpha;
  double _beta;
  double _gamma;
};
