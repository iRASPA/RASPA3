module;

#ifdef USE_LEGACY_HEADERS
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

export module skrotationmatrix;

#ifndef USE_LEGACY_HEADERS
import <ostream>;
import <vector>;
import <type_traits>;
import <tuple>;
#endif

import int3;
import double3;
import int3x3;
import double3x3;
import hashcombine;

export struct SKRotationMatrix
{
  int3x3 int3x3_m;

  enum class RotationType : std::make_signed_t<std::size_t>
  {
    axis_6m = -6,
    axis_4m = -4,
    axis_3m = -3,
    axis_2m = -2,
    axis_1m = -1,
    none = 0,
    axis_1 = 1,
    axis_2 = 2,
    axis_3 = 3,
    axis_4 = 4,
    axis_6 = 6
  };

  enum class SymmetryType : std::size_t
  {
    unknown = 0,
    identity = 1,
    translation = 2,
    inversion = 3,
    pure_rotation = 4,
    pure_reflection = 5,
    screw_rotation = 6,
    glide_reflection = 7
  };

  SKRotationMatrix();
  SKRotationMatrix(int3x3 m) { this->int3x3_m = m; }
  SKRotationMatrix(int3 v1, int3 v2, int3 v3);

  SKRotationMatrix inverse();
  const SKRotationMatrix proper() const;
  SKRotationMatrix::RotationType type() const;
  int3 rotationAxis() const;
  std::vector<int3> orthogonalToAxisDirection(size_t rotationOrder);

  inline int determinant() { return this->int3x3_m.determinant(); }
  inline SKRotationMatrix operator-() const { return -this->int3x3_m; }
  inline bool operator==(const SKRotationMatrix& b) const { return this->int3x3_m == b.int3x3_m; }

  struct hashFunction
  {
    std::size_t operator()(const SKRotationMatrix& k) const
    {
      std::size_t h = 0;
      hash_combine(h, k.int3x3_m.m11);
      hash_combine(h, k.int3x3_m.m12);
      hash_combine(h, k.int3x3_m.m13);
      hash_combine(h, k.int3x3_m.m21);
      hash_combine(h, k.int3x3_m.m22);
      hash_combine(h, k.int3x3_m.m23);
      hash_combine(h, k.int3x3_m.m31);
      hash_combine(h, k.int3x3_m.m32);
      hash_combine(h, k.int3x3_m.m33);
      return h;
    }
  };
  friend std::ostream& operator<<(std::ostream& os, const SKRotationMatrix& setting);

  static SKRotationMatrix zero;
  static SKRotationMatrix identity;
  static SKRotationMatrix inversionIdentity;

  // rotations for principle axes
  static SKRotationMatrix r_2_100;
  static SKRotationMatrix r_2i_100;
  static SKRotationMatrix r_3_100;
  static SKRotationMatrix r_3i_100;
  static SKRotationMatrix r_4_100;
  static SKRotationMatrix r_4i_100;
  static SKRotationMatrix r_6_100;
  static SKRotationMatrix r_6i_100;

  static SKRotationMatrix r_2_010;
  static SKRotationMatrix r_2i_010;
  static SKRotationMatrix r_3_010;
  static SKRotationMatrix r_3i_010;
  static SKRotationMatrix r_4_010;
  static SKRotationMatrix r_4i_010;
  static SKRotationMatrix r_6_010;
  static SKRotationMatrix r_6i_010;

  static SKRotationMatrix r_2_001;
  static SKRotationMatrix r_2i_001;
  static SKRotationMatrix r_3_001;
  static SKRotationMatrix r_3i_001;
  static SKRotationMatrix r_4_001;
  static SKRotationMatrix r_4i_001;
  static SKRotationMatrix r_6_001;
  static SKRotationMatrix r_6i_001;

  static SKRotationMatrix r_3_111;
  static SKRotationMatrix r_3i_111;

  static SKRotationMatrix r_2prime_100;
  static SKRotationMatrix r_2iprime_100;
  static SKRotationMatrix r_2doubleprime_100;
  static SKRotationMatrix r_2idoubleprime_100;

  static SKRotationMatrix r_2prime_010;
  static SKRotationMatrix r_2iprime_010;
  static SKRotationMatrix r_2doubleprime_010;
  static SKRotationMatrix r_2idoubleprime_010;

  static SKRotationMatrix r_2prime_001;
  static SKRotationMatrix r_2iprime_001;
  static SKRotationMatrix r_2doubleprime_001;
  static SKRotationMatrix r_2idoubleprime_001;

  static SKRotationMatrix monoclinicB1toA1;
  static SKRotationMatrix monoclinicB1toA2;
  static SKRotationMatrix monoclinicB1toA3;
  static SKRotationMatrix monoclinicB1toB2;
  static SKRotationMatrix monoclinicB1toB3;
  static SKRotationMatrix monoclinicB1toC1;
  static SKRotationMatrix monoclinicB1toC2;
  static SKRotationMatrix monoclinicB1toC3;

  static SKRotationMatrix orthorhombicCABtoABC;
  static SKRotationMatrix orthorhombicBCAtoABC;
  static SKRotationMatrix orthorhombicBAmCtoABC;
  static SKRotationMatrix orthorhombicAmCBtoABC;
  static SKRotationMatrix orthorhombicmCBAtoABC;

  static std::vector<std::tuple<SKRotationMatrix, int3, int3>> twoFoldSymmetryOperations;
  static std::vector<int3> allPossibleRotationAxes;
};

export inline SKRotationMatrix operator*(const SKRotationMatrix& a, const SKRotationMatrix& b)
{
  return SKRotationMatrix(a.int3x3_m * b.int3x3_m);
}

export inline int3 operator*(const SKRotationMatrix& a, const int3& b) { return a.int3x3_m * b; }

export inline double3 operator*(const SKRotationMatrix& a, const double3& b) { return a.int3x3_m * b; }

export inline SKRotationMatrix operator+(const SKRotationMatrix& a, const SKRotationMatrix& b)
{
  return a.int3x3_m + b.int3x3_m;
}

export inline double3x3 operator*(const double3x3& v1, const SKRotationMatrix& v2)
{
  const int3x3 value = v2.int3x3_m;
  return v1 * double3x3(value);
}

export inline SKRotationMatrix operator-(const SKRotationMatrix& v1, const SKRotationMatrix& v2)
{
  return SKRotationMatrix(v1.int3x3_m - v2.int3x3_m);
}

/*
namespace std
{
    template <> struct hash<SKRotationMatrix>
    {
        size_t operator()(const SKRotationMatrix& k) const
        {
            std::hash<int3x3> hash_fn;
            return hash_fn(k.int3x3_m);
        }
    };
}*/
