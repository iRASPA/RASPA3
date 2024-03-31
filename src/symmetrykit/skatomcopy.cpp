module;


module skatomcopy;

import <memory>;
import double3;

SKAtomCopy::SKAtomCopy(const SKAtomCopy& atomCopy)
{
    this->_position = atomCopy._position;
    this->_type = atomCopy._type;
    this->_tag = atomCopy._tag;
    this->_asymmetricIndex = atomCopy._asymmetricIndex;
    //this->_bonds = {};
}
