module;

module skasymmetricatom;

import bool3;
import double3;
import double2;
import double4;
import skelement;
import skatomcopy;

import <vector>;
import <string>;
import <memory>;

SKAsymmetricAtom::SKAsymmetricAtom(std::string displayName, size_t elementIdentifier) : _displayName(displayName), _elementIdentifier(elementIdentifier)
{
    _uniqueForceFieldName = PredefinedElements::predefinedElements[elementIdentifier]._chemicalSymbol;
}

SKAsymmetricAtom::SKAsymmetricAtom(std::string displayName, size_t elementIdentifier, double occupancy) : _displayName(displayName), _elementIdentifier(elementIdentifier), _occupancy(occupancy)
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

void SKAsymmetricAtom::toggleVisibility()
{
    _isVisible = !_isVisible;
}
