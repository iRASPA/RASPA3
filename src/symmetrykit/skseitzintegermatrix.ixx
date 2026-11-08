module;

export module skseitzintegermatrix;

import std;

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
  static std::vector<SKSeitzIntegerMatrix> SeitzMatricesFromGenerators(std::string encoding);

  // the translation in units of 1/24 wrapped into [0, 24): the signed remainder of C++ would leave a negative
  // component negative (-6 % 24 is -6, not 18), so that the products of operations with -1 rotation entries compared
  // unequal to, and hashed apart from, the same operation written with a positive translation
  int3 normalizedTranslation() const
  {
    return int3(((translation.x % 24) + 24) % 24, ((translation.y % 24) + 24) % 24, ((translation.z % 24) + 24) % 24);
  }

  bool operator==(SKSeitzIntegerMatrix const& rhs) const
  {
    return (this->rotation == rhs.rotation) && (this->normalizedTranslation() == rhs.normalizedTranslation());
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
      // the translation as operator== compares it (wrapped into [0, 24)), so that equal matrices hash equally
      const int3 translation = k.normalizedTranslation();
      hash_combine(h, translation.x);
      hash_combine(h, translation.y);
      hash_combine(h, translation.z);
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
