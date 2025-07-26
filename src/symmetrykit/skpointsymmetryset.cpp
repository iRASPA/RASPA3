module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

module skpointsymmetryset;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import skrotationmatrix;

SKPointSymmetrySet::SKPointSymmetrySet(std::vector<SKRotationMatrix> rotations) : _rotations(rotations) {}
