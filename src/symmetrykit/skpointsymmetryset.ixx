module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

export module skpointsymmetryset;

#ifndef USE_LEGACY_HEADERS
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
