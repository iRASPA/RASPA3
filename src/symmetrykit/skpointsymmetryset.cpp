module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#endif

module skpointsymmetryset;

#ifndef USE_LEGACY_HEADERS
import <vector>;
#endif

import skrotationmatrix;


SKPointSymmetrySet::SKPointSymmetrySet(std::vector<SKRotationMatrix> rotations) : _rotations(rotations)
{
}
