module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

export module skrotationalchangeofbasis;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import skrotationmatrix;
import skseitzintegermatrix;

export struct SKRotationalChangeOfBasis
{
  SKRotationMatrix rotationMatrix;
  SKRotationMatrix inverseRotationMatrix;

  SKRotationalChangeOfBasis(SKRotationMatrix rotationMatrix);
  SKRotationalChangeOfBasis(SKRotationMatrix rotationMatrix, SKRotationMatrix inverseRotationMatrix);

  static SKRotationalChangeOfBasis identity;
  static std::vector<SKRotationalChangeOfBasis> changeOfMonoclinicCentering;
  static std::vector<SKRotationalChangeOfBasis> changeOfOrthorhombicCentering;

  SKRotationalChangeOfBasis inverse();
};

export inline SKSeitzIntegerMatrix operator*(const SKRotationalChangeOfBasis& a, const SKSeitzIntegerMatrix& b)
{
  SKRotationMatrix rotationMatrix = SKRotationMatrix(a.inverseRotationMatrix * b.rotation * a.rotationMatrix);
  int3 translationVector = a.inverseRotationMatrix * b.translation;
  return SKSeitzIntegerMatrix(rotationMatrix, translationVector);
}

export inline double3 operator*(const SKRotationalChangeOfBasis& a, const double3& b)
{
  return a.inverseRotationMatrix * b;
}

export inline int3 operator*(const SKRotationalChangeOfBasis& a, const int3& b) { return a.inverseRotationMatrix * b; }
