module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

module skpointsymmetryset;

#ifdef USE_STD_IMPORT
import std;
#endif

import skrotationmatrix;

SKPointSymmetrySet::SKPointSymmetrySet(std::vector<SKRotationMatrix> rotations) : _rotations(rotations) {}
