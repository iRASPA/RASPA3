module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

export module skpointsymmetryset;

#ifdef USE_STD_IMPORT
import std;
#endif

import skrotationmatrix;

export class SKPointSymmetrySet
{
 public:
  SKPointSymmetrySet(std::vector<SKRotationMatrix> rotations);
  const std::vector<SKRotationMatrix>& rotations() { return _rotations; }

 private:
  std::vector<SKRotationMatrix> _rotations;
};
