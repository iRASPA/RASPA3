module;

#ifdef USE_LEGACY_HEADERS
#include <memory>
#endif

module skatomcopy;

#ifndef USE_LEGACY_HEADERS
import <memory>;
#endif

import double3;

SKAtomCopy::SKAtomCopy(const SKAtomCopy& atomCopy)
{
  this->_position = atomCopy._position;
  this->_type = atomCopy._type;
  this->_tag = atomCopy._tag;
  this->_asymmetricIndex = atomCopy._asymmetricIndex;
  // this->_bonds = {};
}
