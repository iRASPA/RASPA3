module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#endif

class SKAtomCopy;

export module skatom:skasymmetricatom;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import bool3;
import double3;
import double2;
import double4;

import skelement;

export class SKAsymmetricAtom
{
 public:
  SKAsymmetricAtom() : _displayName("C"), _elementIdentifier(6) {}
  SKAsymmetricAtom(const SKAsymmetricAtom& assymetricAtom);
  SKAsymmetricAtom(std::string displayName, std::size_t elementIdentifier);
  SKAsymmetricAtom(std::string displayName, std::size_t elementIdentifier, double occupancy);

  enum class SKAsymmetricAtomType : std::size_t
  {
    container = 0,
    asymmetric = 1
  };

  enum class Hybridization : std::size_t
  {
    untyped = 0,
    sp_linear = 1,
    sp2_trigonal = 2,
    sp3_tetrahedral = 3,
    square_planar = 4,
    trigonal_bipyramidal = 5,
    square_pyramidal = 6,
    octahedral = 7
  };

  std::size_t asymmetricIndex() { return _asymmetricIndex; }
  void setAsymmetricIndex(std::size_t value) { _asymmetricIndex = value; }

  std::string displayName() const { return _displayName; }
  void setDisplayName(std::string newValue) { _displayName = newValue; }
  double3 position() const { return _position; }
  void setPosition(double3 newValue) { _position = newValue; }
  void setPositionX(double newValue) { _position.x = newValue; }
  void setPositionY(double newValue) { _position.y = newValue; }
  void setPositionZ(double newValue) { _position.z = newValue; }
  double charge() const { return _charge; }
  void setCharge(double newValue) { _charge = newValue; }

  std::size_t tag() { return _tag; }
  void setTag(std::size_t tag) { _tag = tag; }
  bool isVisible() { return _isVisible; }
  void toggleVisibility();
  void setVisibility(bool visibility) { _isVisible = visibility; }

  double4 color() { return _color; }
  void setColor(double4 color) { _color = color; }
  double drawRadius() { return _drawRadius; }
  void setDrawRadius(double radius) { _drawRadius = radius; }

  double2 potentialParameters() { return _potentialParameters; }
  void setPotentialParameters(double2 value) { _potentialParameters = value; }

  std::string uniqueForceFieldName() const { return _uniqueForceFieldName; }
  void setUniqueForceFieldName(std::string newValue) { _uniqueForceFieldName = newValue; }
  std::size_t elementIdentifier() const { return _elementIdentifier; }
  void setElementIdentifier(std::size_t newValue) { _elementIdentifier = newValue; }

  double bondDistanceCriteria() const { return _bondDistanceCriteria; }
  void setBondDistanceCriteria(double bondDistanceCriteria) { _bondDistanceCriteria = bondDistanceCriteria; }

  bool3 isFixed() const { return _isFixed; }
  void setIsFixed(bool3 newValue) { _isFixed = newValue; }
  std::size_t serialNumber() const { return _serialNumber; }
  void setSerialNumber(std::size_t newValue) { _serialNumber = newValue; }
  char remotenessIndicator() const { return _remotenessIndicator; }
  void setRemotenessIndicator(char newValue) { _remotenessIndicator = newValue; }
  char branchDesignator() const { return _branchDesignator; }
  void setBranchDesignator(char newValue) { _branchDesignator = newValue; }
  char alternateLocationIndicator() const { return _alternateLocationIndicator; }
  void setAlternateLocationIndicator(char newValue) { _alternateLocationIndicator = newValue; }
  std::string residueName() const { return _residueName; }
  void setResidueName(std::string newValue) { _residueName = newValue; }
  char chainIdentifier() const { return _chainIdentifier; }
  void setChainIdentifier(char newValue) { _chainIdentifier = newValue; }
  std::size_t residueSequenceNumber() const { return _residueSequenceNumber; }
  void setResidueSequenceNumber(std::size_t newValue) { _residueSequenceNumber = newValue; }
  char codeForInsertionOfResidues() const { return _codeForInsertionOfResidues; }
  void setCodeForInsertionOfResidues(char newValue) { _codeForInsertionOfResidues = newValue; }
  double occupancy() const { return _occupancy; }
  void setOccupancy(double newValue) { _occupancy = newValue; }
  double temperaturefactor() const { return _temperaturefactor; }
  void setTemperaturefactor(double newValue) { _temperaturefactor = newValue; }
  std::string segmentIdentifier() const { return _segmentIdentifier; }
  void setSegmentIdentifier(std::string newValue) { _segmentIdentifier = newValue; }
  std::size_t asymetricID() const { return _asymetricID; }
  // void setAsymetricID(int newValue)  {_asymetricID = newValue;}

  bool ligandAtom() const { return _ligandAtom; }
  void setLigandAtom(bool newValue) { _ligandAtom = newValue; }
  bool backBoneAtom() const { return _backBoneAtom; }
  void backBoneAtom(bool newValue) { _backBoneAtom = newValue; }
  bool fractional() const { return _fractional; }
  void fractional(bool newValue) { _fractional = newValue; }
  bool solvent() const { return _solvent; }
  void setSolvent(bool newValue) { _solvent = newValue; }

  std::vector<std::shared_ptr<SKAtomCopy>>& copies() { return _copies; }
  void setCopies(std::vector<std::shared_ptr<SKAtomCopy>> copies) { _copies = copies; }

 private:
  std::size_t _versionNumber{2};
  std::size_t _asymmetricIndex;
  std::string _displayName = std::string("Default");
  double3 _position = double3(0, 0, 0);
  double _charge = 0;

  std::string _uniqueForceFieldName;
  std::size_t _elementIdentifier = 0;
  double4 _color = double4(0, 1.0, 0, 1.0);
  double _drawRadius = 1.0;
  double _bondDistanceCriteria = 1.0;
  double2 _potentialParameters = double2(0.0, 0.0);

  std::size_t _tag = 0;
  [[maybe_unused]] SKAsymmetricAtomType _symmetryType = SKAsymmetricAtomType::asymmetric;
  Hybridization _hybridization = Hybridization::untyped;

  // atom properties (bonds are visible depending on whether the atoms of the bonds are visible)
  bool3 _isFixed = bool3(false, false, false);
  bool _isVisible = true;
  bool _isVisibleEnabled = true;

  std::size_t _serialNumber = 0;
  char _remotenessIndicator = ' ';             // character 'A','B','C','D',...
  char _branchDesignator = ' ';                // character '1','2','3',...
  std::size_t _asymetricID = 0;                // positive integer
  char _alternateLocationIndicator = ' ';      // character ' ' or 'A','B',...
  std::string _residueName = std::string("");  // empty or 3 characters
  char _chainIdentifier = ' ';                 // empty or 'A','B','C',...
  std::size_t _residueSequenceNumber = 0;      // positive integer
  char _codeForInsertionOfResidues = ' ';      // empty or 'A','B','C',...
  double _occupancy = 1.0;
  double _temperaturefactor = 0.0;
  std::string _segmentIdentifier = std::string("");  // empty or 4 characters

  bool _ligandAtom = false;
  bool _backBoneAtom = false;
  bool _fractional = false;
  bool _solvent = false;

  double3 _displacement = double3();

  // the crystallographic copies of the atom
  std::vector<std::shared_ptr<SKAtomCopy>> _copies;
};

SKAsymmetricAtom::SKAsymmetricAtom(std::string displayName, std::size_t elementIdentifier)
    : _displayName(displayName), _elementIdentifier(elementIdentifier)
{
  _uniqueForceFieldName = PredefinedElements::predefinedElements[elementIdentifier]._chemicalSymbol;
}

SKAsymmetricAtom::SKAsymmetricAtom(std::string displayName, std::size_t elementIdentifier, double occupancy)
    : _displayName(displayName), _elementIdentifier(elementIdentifier), _occupancy(occupancy)
{
  _uniqueForceFieldName = PredefinedElements::predefinedElements[elementIdentifier]._chemicalSymbol;
}

SKAsymmetricAtom::SKAsymmetricAtom(const SKAsymmetricAtom& asymmetricAtom)
{
  _versionNumber = asymmetricAtom._versionNumber;
  _asymmetricIndex = asymmetricAtom._asymmetricIndex;
  _displayName = asymmetricAtom._displayName;
  _position = asymmetricAtom._position;
  _charge = asymmetricAtom._charge;
  _hybridization = asymmetricAtom._hybridization;
  _uniqueForceFieldName = asymmetricAtom._uniqueForceFieldName;
  _elementIdentifier = asymmetricAtom._elementIdentifier;
  _color = asymmetricAtom._color;
  _drawRadius = asymmetricAtom._drawRadius;

  _bondDistanceCriteria = asymmetricAtom._bondDistanceCriteria;
  _potentialParameters = asymmetricAtom._potentialParameters;
  _tag = asymmetricAtom._tag;
  _isFixed = asymmetricAtom._isFixed;
  _isVisible = asymmetricAtom._isVisible;
  _isVisibleEnabled = asymmetricAtom._isVisibleEnabled;

  _serialNumber = asymmetricAtom._serialNumber;
  _remotenessIndicator = asymmetricAtom._remotenessIndicator;
  _branchDesignator = asymmetricAtom._branchDesignator;
  _asymetricID = asymmetricAtom._asymetricID;
  _alternateLocationIndicator = asymmetricAtom._alternateLocationIndicator;
  _residueName = asymmetricAtom._residueName;
  _chainIdentifier = asymmetricAtom._chainIdentifier;
  _residueSequenceNumber = asymmetricAtom._residueSequenceNumber;
  _codeForInsertionOfResidues = asymmetricAtom._codeForInsertionOfResidues;
  _occupancy = asymmetricAtom._occupancy;
  _temperaturefactor = asymmetricAtom._temperaturefactor;
  _segmentIdentifier = asymmetricAtom._segmentIdentifier;

  _ligandAtom = asymmetricAtom._ligandAtom;
  _backBoneAtom = asymmetricAtom._backBoneAtom;
  _fractional = asymmetricAtom._fractional;
  _solvent = asymmetricAtom._solvent;

  _copies = std::vector<std::shared_ptr<SKAtomCopy>>{};
}

void SKAsymmetricAtom::toggleVisibility() { _isVisible = !_isVisible; }
