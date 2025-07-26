module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <string>
#include <vector>
#endif

export module skseitzintegermatrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import skrotationmatrix;
import int3;
import int3x3;
import double3;
import double3x3;
import skonethirdseitzmatrix;
import hashcombine;

export struct SKSeitzIntegerMatrix
{
  SKRotationMatrix rotation;
  int3 translation;  // denominator = 24

  SKSeitzIntegerMatrix();
  SKSeitzIntegerMatrix(SKRotationMatrix rotation, int3 translation);
  SKSeitzIntegerMatrix(char xvalue, char yvalue, char zvalue);
  static std::vector<SKSeitzIntegerMatrix> SeitzMatrices(std::string encoding);

  int3 normalizedTranslation() const { return int3(translation.x % 24, translation.y % 24, translation.z % 24); }

  bool operator==(SKSeitzIntegerMatrix const& rhs) const
  {
    return (this->rotation == rhs.rotation) && ((this->translation.x % 24) == (rhs.translation.x % 24)) &&
           ((this->translation.y % 24) == (rhs.translation.y % 24)) &&
           ((this->translation.z % 24) == (rhs.translation.z % 24));
  }

  struct hashFunction
  {
    std::size_t operator()(const SKSeitzIntegerMatrix& k) const
    {
      std::size_t h = 0;
      hash_combine(h, k.rotation.int3x3_m.m11);
      hash_combine(h, k.rotation.int3x3_m.m12);
      hash_combine(h, k.rotation.int3x3_m.m13);
      hash_combine(h, k.rotation.int3x3_m.m21);
      hash_combine(h, k.rotation.int3x3_m.m22);
      hash_combine(h, k.rotation.int3x3_m.m23);
      hash_combine(h, k.rotation.int3x3_m.m31);
      hash_combine(h, k.rotation.int3x3_m.m32);
      hash_combine(h, k.rotation.int3x3_m.m33);
      hash_combine(h, k.translation.x);
      hash_combine(h, k.translation.y);
      hash_combine(h, k.translation.z);
      return h;
    }
  };

  static std::vector<SKOneThirdSeitzMatrix> SeitzData;
};

export inline double3 operator*(const SKSeitzIntegerMatrix& a, const double3& b)
{
  double3 v = a.rotation * b;
  double3 trans = double3(static_cast<double>(a.translation.x) / 24.0, static_cast<double>(a.translation.y) / 24.0,
                          static_cast<double>(a.translation.z) / 24.0);
  return v + trans;
}

export inline SKSeitzIntegerMatrix operator*(const SKSeitzIntegerMatrix& a, const SKSeitzIntegerMatrix& b)
{
  SKRotationMatrix rotationMatrix = a.rotation * b.rotation;
  int3 translationVector = a.translation + a.rotation * b.translation;
  return SKSeitzIntegerMatrix(rotationMatrix, translationVector);
}
