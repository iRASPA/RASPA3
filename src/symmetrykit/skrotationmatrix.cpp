module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <tuple>
#include <vector>
#endif

module skrotationmatrix;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <iostream>;
import <tuple>;
#endif

import int3;
import int3x3;
import double3;
import double3x3;

SKRotationMatrix::SKRotationMatrix() {}

// SKRotationMatrix::SKRotationMatrix(const SKTransformationMatrix &m)
//{
//     assert(abs(m.determinant()) == 1);
//     this->int3x3_m = m.transformation;
// }

SKRotationMatrix::SKRotationMatrix(int3 v1, int3 v2, int3 v3)
{
  this->int3x3_m.m11 = v1.x;
  this->int3x3_m.m21 = v1.y;
  this->int3x3_m.m31 = v1.z;
  this->int3x3_m.m12 = v2.x;
  this->int3x3_m.m22 = v2.y;
  this->int3x3_m.m32 = v2.z;
  this->int3x3_m.m13 = v3.x;
  this->int3x3_m.m23 = v3.y;
  this->int3x3_m.m33 = v3.z;
}

const SKRotationMatrix SKRotationMatrix::proper() const
{
  if (this->int3x3_m.determinant() == 1)
  {
    return *this;
  }
  else
  {
    return -(*this);
  }
}

SKRotationMatrix::RotationType SKRotationMatrix::type() const
{
  int determinant = this->int3x3_m.determinant();

  if (determinant == -1)
  {
    switch (this->int3x3_m.trace())
    {
      case -3:
        return RotationType::axis_1m;
      case -2:
        return RotationType::axis_6m;
      case -1:
        return RotationType::axis_4m;
      case 0:
        return RotationType::axis_3m;
      case 1:
        return RotationType::axis_2m;

      default:
        return RotationType::none;
    }
  }
  else
  {
    switch (this->int3x3_m.trace())
    {
      case -1:
        return RotationType::axis_2;
      case 0:
        return RotationType::axis_3;
      case 1:
        return RotationType::axis_4;
      case 2:
        return RotationType::axis_6;
      case 3:
        return RotationType::axis_1;
      default:
        return RotationType::none;
    }
  }
}

int3 SKRotationMatrix::rotationAxis() const
{
  // rotation axis is the eigenvector with eigenvalue lambda==1
  for (size_t i = 0; i < SKRotationMatrix::allPossibleRotationAxes.size(); i++)
  {
    int3 axis = SKRotationMatrix::allPossibleRotationAxes[i];
    if ((this->int3x3_m * axis) == axis)
    {
      return axis;
    }
  }

  return int3(0, 0, 0);
}

/// Computes a list of integer-vectors that are orthogonal to the rotation axis for a given rotation matrix
///
/// - parameter rotationOrder: the rotation order
///
/// - returns: a list of perpendicular eigenvectors
///
/// Note : Theorem TA4.1 in Boisen en Gibbs (1990) states that a vector x is in the plane perpendicular to the axis
/// direction 'e' of a proper rotation matrix 'W_p' with rotational order 'n' if and only if S.x=0 where S = W_p + W_p^2
/// + ... + W_p^n Ref: R.W. Grosse-Kunstleve, "Algorithms for deriving crystallographic space-group information", Acta
/// Cryst. A55, 383-395, 1999
///
/// The algorithm of Atsushi Togo is used: a search over all possible rotation axes.
std::vector<int3> SKRotationMatrix::orthogonalToAxisDirection(size_t rotationOrder)
{
  std::vector<int3> orthoAxes{};

  SKRotationMatrix properRotation = this->proper();
  SKRotationMatrix sumRot = SKRotationMatrix::identity;
  SKRotationMatrix rot = SKRotationMatrix::identity;

  for (size_t i = 0; i < rotationOrder - 1; i++)
  {
    rot = rot * properRotation;
    sumRot = sumRot + rot;
  }

  for (const int3& rotationAxes : SKRotationMatrix::allPossibleRotationAxes)
  {
    if (sumRot * rotationAxes == int3(0, 0, 0))
    {
      orthoAxes.push_back(rotationAxes);
    }
  }

  return orthoAxes;
}

std::ostream& operator<<(std::ostream& os, const SKRotationMatrix& m)
{
  os << "SKRotationMatrix: " << m.int3x3_m.m11 << '/' << m.int3x3_m.m12 << '/' << m.int3x3_m.m13 << '/'
     << m.int3x3_m.m21 << '/' << m.int3x3_m.m22 << '/' << m.int3x3_m.m23 << '/' << m.int3x3_m.m31 << '/'
     << m.int3x3_m.m32 << '/' << m.int3x3_m.m33 << '/';
  return os;
}

/// Inverse of the matrix if the determinant = 1 or -1, otherwise the contents of the resulting matrix are undefined.
SKRotationMatrix SKRotationMatrix::inverse()
{
  int determinant = int3x3_m.determinant();

  int3 c1 = int3(-this->int3x3_m[1][2] * int3x3_m[2][1] + this->int3x3_m[1][1] * this->int3x3_m[2][2],
                 this->int3x3_m[0][2] * this->int3x3_m[2][1] - this->int3x3_m[0][1] * this->int3x3_m[2][2],
                 -this->int3x3_m[0][2] * this->int3x3_m[1][1] + this->int3x3_m[0][1] * this->int3x3_m[1][2]);
  int3 c2 = int3(this->int3x3_m[1][2] * int3x3_m[2][0] - this->int3x3_m[1][0] * this->int3x3_m[2][2],
                 -this->int3x3_m[0][2] * this->int3x3_m[2][0] + this->int3x3_m[0][0] * this->int3x3_m[2][2],
                 this->int3x3_m[0][2] * this->int3x3_m[1][0] - this->int3x3_m[0][0] * this->int3x3_m[1][2]);
  int3 c3 = int3(-this->int3x3_m[1][1] * int3x3_m[2][0] + this->int3x3_m[1][0] * this->int3x3_m[2][1],
                 this->int3x3_m[0][1] * this->int3x3_m[2][0] - this->int3x3_m[0][0] * this->int3x3_m[2][1],
                 -this->int3x3_m[0][1] * this->int3x3_m[1][0] + this->int3x3_m[0][0] * this->int3x3_m[1][1]);

  switch (determinant)
  {
    case -1:
      return -SKRotationMatrix(c1, c2, c3);
    case 1:
      return SKRotationMatrix(c1, c2, c3);
    default:
      return SKRotationMatrix();
  }
}

SKRotationMatrix SKRotationMatrix::zero = SKRotationMatrix(int3(0, 0, 0), int3(0, 0, 0), int3(0, 0, 0));
SKRotationMatrix SKRotationMatrix::identity = SKRotationMatrix(int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::inversionIdentity = SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, -1));

// rotations for principle axes
SKRotationMatrix SKRotationMatrix::r_2_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, -1, 0), int3(0, 0, -1));
SKRotationMatrix SKRotationMatrix::r_2i_100 = r_2_100;
SKRotationMatrix SKRotationMatrix::r_3_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, 0, 1), int3(0, -1, -1));
SKRotationMatrix SKRotationMatrix::r_3i_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, -1, -1), int3(0, 1, 0));
SKRotationMatrix SKRotationMatrix::r_4_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, 0, 1), int3(0, -1, 0));
SKRotationMatrix SKRotationMatrix::r_4i_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, 0, -1), int3(0, 1, 0));
SKRotationMatrix SKRotationMatrix::r_6_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, 1, 1), int3(0, -1, 0));
SKRotationMatrix SKRotationMatrix::r_6i_100 = SKRotationMatrix(int3(1, 0, 0), int3(0, 0, -1), int3(0, 1, 1));

SKRotationMatrix SKRotationMatrix::r_2_010 = SKRotationMatrix(int3(-1, 0, 0), int3(0, 1, 0), int3(0, 0, -1));
SKRotationMatrix SKRotationMatrix::r_2i_010 = r_2_010;
SKRotationMatrix SKRotationMatrix::r_3_010 = SKRotationMatrix(int3(-1, 0, -1), int3(0, 1, 0), int3(1, 0, 0));
SKRotationMatrix SKRotationMatrix::r_3i_010 = SKRotationMatrix(int3(0, 0, 1), int3(0, 1, 0), int3(-1, 0, -1));
SKRotationMatrix SKRotationMatrix::r_4_010 = SKRotationMatrix(int3(0, 0, -1), int3(0, 1, 0), int3(1, 0, 0));
SKRotationMatrix SKRotationMatrix::r_4i_010 = SKRotationMatrix(int3(0, 0, 1), int3(0, 1, 0), int3(-1, 0, 0));
SKRotationMatrix SKRotationMatrix::r_6_010 = SKRotationMatrix(int3(0, 0, -1), int3(0, 1, 0), int3(1, 0, 1));
SKRotationMatrix SKRotationMatrix::r_6i_010 = SKRotationMatrix(int3(1, 0, 1), int3(0, 1, 0), int3(-1, 0, 0));

SKRotationMatrix SKRotationMatrix::r_2_001 = SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::r_2i_001 = r_2_001;
SKRotationMatrix SKRotationMatrix::r_3_001 = SKRotationMatrix(int3(0, 1, 0), int3(-1, -1, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::r_3i_001 = SKRotationMatrix(int3(-1, -1, 0), int3(1, 0, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::r_4_001 = SKRotationMatrix(int3(0, 1, 0), int3(-1, 0, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::r_4i_001 = SKRotationMatrix(int3(0, -1, 0), int3(1, 0, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::r_6_001 = SKRotationMatrix(int3(1, 1, 0), int3(-1, 0, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::r_6i_001 = SKRotationMatrix(int3(0, -1, 0), int3(1, 1, 0), int3(0, 0, 1));

SKRotationMatrix SKRotationMatrix::r_3_111 = SKRotationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0));
SKRotationMatrix SKRotationMatrix::r_3i_111 = SKRotationMatrix(int3(0, 0, 1), int3(1, 0, 0), int3(0, 1, 0));

SKRotationMatrix SKRotationMatrix::r_2prime_100 =
    SKRotationMatrix(int3(-1, 0, 0), int3(0, 0, -1), int3(0, -1, 0));  // b-c
SKRotationMatrix SKRotationMatrix::r_2iprime_100 = r_2prime_100;
SKRotationMatrix SKRotationMatrix::r_2doubleprime_100 =
    SKRotationMatrix(int3(-1, 0, 0), int3(0, 0, 1), int3(0, 1, 0));  // b+c
SKRotationMatrix SKRotationMatrix::r_2idoubleprime_100 = r_2doubleprime_100;

SKRotationMatrix SKRotationMatrix::r_2prime_010 =
    SKRotationMatrix(int3(0, 0, -1), int3(0, -1, 0), int3(-1, 0, 0));  // a-c
SKRotationMatrix SKRotationMatrix::r_2iprime_010 = r_2prime_010;
SKRotationMatrix SKRotationMatrix::r_2doubleprime_010 =
    SKRotationMatrix(int3(0, 0, 1), int3(0, -1, 0), int3(1, 0, 0));  // a+c
SKRotationMatrix SKRotationMatrix::r_2idoubleprime_010 = r_2doubleprime_010;

SKRotationMatrix SKRotationMatrix::r_2prime_001 =
    SKRotationMatrix(int3(0, -1, 0), int3(-1, 0, 0), int3(0, 0, -1));  // a-b
SKRotationMatrix SKRotationMatrix::r_2iprime_001 = r_2prime_001;
SKRotationMatrix SKRotationMatrix::r_2doubleprime_001 =
    SKRotationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(0, 0, -1));  // a+b
SKRotationMatrix SKRotationMatrix::r_2idoubleprime_001 = r_2doubleprime_001;

SKRotationMatrix SKRotationMatrix::monoclinicB1toA1 = SKRotationMatrix(int3(0, 0, 1), int3(1, 0, 0), int3(0, 1, 0));
SKRotationMatrix SKRotationMatrix::monoclinicB1toA2 = SKRotationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(0, -1, -1));
SKRotationMatrix SKRotationMatrix::monoclinicB1toA3 = SKRotationMatrix(int3(0, -1, -1), int3(1, 0, 0), int3(0, 0, 1));
SKRotationMatrix SKRotationMatrix::monoclinicB1toB2 = SKRotationMatrix(int3(0, 0, 1), int3(0, 1, 0), int3(-1, 0, -1));
SKRotationMatrix SKRotationMatrix::monoclinicB1toB3 = SKRotationMatrix(int3(-1, 0, -1), int3(0, 1, 0), int3(1, 0, 0));
SKRotationMatrix SKRotationMatrix::monoclinicB1toC1 = SKRotationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0));
SKRotationMatrix SKRotationMatrix::monoclinicB1toC2 = SKRotationMatrix(int3(1, 0, 0), int3(0, 0, 1), int3(-1, -1, 0));
SKRotationMatrix SKRotationMatrix::monoclinicB1toC3 = SKRotationMatrix(int3(-1, -1, 0), int3(0, 0, 1), int3(0, 1, 0));

SKRotationMatrix SKRotationMatrix::orthorhombicCABtoABC = SKRotationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0));
SKRotationMatrix SKRotationMatrix::orthorhombicBCAtoABC = SKRotationMatrix(int3(0, 0, 1), int3(1, 0, 0), int3(0, 1, 0));
SKRotationMatrix SKRotationMatrix::orthorhombicBAmCtoABC =
    SKRotationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(0, 0, -1));
SKRotationMatrix SKRotationMatrix::orthorhombicAmCBtoABC =
    SKRotationMatrix(int3(1, 0, 0), int3(0, 0, 1), int3(0, -1, 0));
SKRotationMatrix SKRotationMatrix::orthorhombicmCBAtoABC =
    SKRotationMatrix(int3(0, 0, 1), int3(0, 1, 0), int3(-1, 0, 0));

std::vector<std::tuple<SKRotationMatrix, int3, int3>> SKRotationMatrix::twoFoldSymmetryOperations = {
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(-1, 0, 1), int3(-1, 1, 0)), int3(-1, 1, 1), int3(0, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(-1, 1, -1), int3(0, 0, -1)), int3(1, -2, 1), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(-1, 1, 0), int3(0, 0, -1)), int3(-1, 2, 0), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(-1, 1, 1), int3(0, 0, -1)), int3(-1, 2, 1), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(-1, 0, -1), int3(1, -1, 0)), int3(1, -1, 1), int3(0, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(-1, -1, 1)), int3(-1, -1, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(-1, 0, 1)), int3(-1, 0, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(-1, 1, 1)), int3(-1, 1, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, -1, 1), int3(0, 0, -1), int3(0, -1, 0)), int3(0, -1, 1), int3(1, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, -1, -1), int3(0, 0, 1), int3(0, 1, 0)), int3(0, 1, 1), int3(-1, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, -1, 0), int3(0, 1, 0), int3(0, -1, -1)), int3(0, 1, 0), int3(1, -2, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, -1, 0), int3(0, 1, 0), int3(0, 0, -1)), int3(0, 1, 0), int3(-1, 2, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, -1, 0), int3(0, 1, 0), int3(0, 1, -1)), int3(0, 1, 0), int3(-1, 2, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, -1, 1)), int3(0, -1, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, -1), int3(0, -1, -1), int3(0, 0, 1)), int3(0, 0, 1), int3(-1, -1, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, -1), int3(0, -1, 0), int3(0, 0, 1)), int3(0, 0, 1), int3(-1, 0, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, -1), int3(0, -1, 1), int3(0, 0, 1)), int3(0, 0, 1), int3(-1, 1, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, -1), int3(0, 0, 1)), int3(0, 0, 1), int3(0, -1, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, 1)), int3(0, 0, 1), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 1), int3(0, 0, 1)), int3(0, 0, 1), int3(0, 1, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 1), int3(0, -1, -1), int3(0, 0, 1)), int3(0, 0, 1), int3(1, -1, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 1), int3(0, -1, 0), int3(0, 0, 1)), int3(0, 0, 1), int3(1, 0, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 1), int3(0, -1, 1), int3(0, 0, 1)), int3(0, 0, 1), int3(1, 1, 2)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 1, 1)), int3(0, 1, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 0, -1), int3(0, -1, 0)), int3(0, -1, 1), int3(0, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 0, 1), int3(0, 1, 0)), int3(0, 1, 1), int3(0, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 1, 0), int3(0, -1, -1)), int3(0, 1, 0), int3(0, -2, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 1, -1), int3(0, 0, -1)), int3(0, -2, 1), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 1, 0), int3(0, 0, -1)), int3(0, 1, 0), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 1, 1), int3(0, 0, -1)), int3(0, 2, 1), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, 1, 0), int3(0, 1, -1)), int3(0, 1, 0), int3(0, 2, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 1, -1), int3(0, 0, -1), int3(0, -1, 0)), int3(0, -1, 1), int3(-1, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 1, 1), int3(0, 0, 1), int3(0, 1, 0)), int3(0, 1, 1), int3(1, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 1, 0), int3(0, 1, 0), int3(0, -1, -1)), int3(0, 1, 0), int3(-1, -2, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 1, 0), int3(0, 1, 0), int3(0, 0, -1)), int3(0, 1, 0), int3(1, 2, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 1, 0), int3(0, 1, 0), int3(0, 1, -1)), int3(0, 1, 0), int3(1, 2, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(1, -1, 1)), int3(1, -1, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(1, 0, 1)), int3(1, 0, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(1, 1, 1)), int3(1, 1, 2), int3(0, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(1, 0, -1), int3(-1, -1, 0)), int3(-1, -1, 1), int3(0, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(1, 1, -1), int3(0, 0, -1)), int3(-1, -2, 1), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(1, 1, 0), int3(0, 0, -1)), int3(1, 2, 0), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(1, 1, 1), int3(0, 0, -1)), int3(1, 2, 1), int3(0, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(-1, 0, 0), int3(1, 0, 1), int3(1, 1, 0)), int3(1, 1, 1), int3(0, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, 0), int3(-1, 0, 0), int3(-1, 1, -1)), int3(-1, 1, 0), int3(-1, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 0, -1), int3(-1, -1, 1), int3(-1, 0, 0)), int3(-1, 0, 1), int3(-1, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, -1), int3(-1, 0, 1), int3(0, 0, -1)), int3(-1, 1, 1), int3(-1, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, 0), int3(-1, 0, 0), int3(0, 0, -1)), int3(-1, 1, 0), int3(-1, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, 1), int3(-1, 0, -1), int3(0, 0, -1)), int3(1, -1, 1), int3(-1, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, 0), int3(-1, 0, 0), int3(1, -1, -1)), int3(-1, 1, 0), int3(1, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 0, 1), int3(-1, -1, -1), int3(1, 0, 0)), int3(1, 0, 1), int3(1, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, -1), int3(0, -1, 0), int3(-1, 1, 0)), int3(-1, 1, 1), int3(-1, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 0, -1), int3(0, -1, 0), int3(-1, 0, 0)), int3(-1, 0, 1), int3(-1, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, -1), int3(0, -1, 0), int3(-1, -1, 0)), int3(-1, -1, 1), int3(-1, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, -1, 1), int3(0, -1, 0), int3(1, -1, 0)), int3(1, -1, 1), int3(1, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 0, 1), int3(0, -1, 0), int3(1, 0, 0)), int3(1, 0, 1), int3(1, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, 1), int3(0, -1, 0), int3(1, 1, 0)), int3(1, 1, 1), int3(1, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 0, -1), int3(1, -1, -1), int3(-1, 0, 0)), int3(-1, 0, 1), int3(-1, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(-1, -1, -1)), int3(1, 1, 0), int3(-1, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, -1), int3(1, 0, -1), int3(0, 0, -1)), int3(-1, -1, 1), int3(1, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(0, 0, -1)), int3(1, 1, 0), int3(1, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, 1), int3(1, 0, 1), int3(0, 0, -1)), int3(1, 1, 1), int3(1, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(0, 0, 1), int3(1, -1, 1), int3(1, 0, 0)), int3(1, 0, 1), int3(1, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(1, 1, -1)), int3(1, 1, 0), int3(1, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(-1, -1, 0), int3(-1, 0, -1)), int3(1, 0, 0), int3(-2, 1, 1)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(-1, -1, 0), int3(0, 0, -1)), int3(1, 0, 0), int3(-2, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(-1, -1, 0), int3(1, 0, -1)), int3(1, 0, 0), int3(2, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(0, -1, 0), int3(-1, 0, -1)), int3(1, 0, 0), int3(-2, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(1, -1, -1), int3(0, -1, 0), int3(0, 0, -1)), int3(-2, 1, 1), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, -1, 0), int3(0, -1, 0), int3(0, 0, -1)), int3(-2, 1, 0), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, -1, 1), int3(0, -1, 0), int3(0, 0, -1)), int3(2, -1, 1), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, -1), int3(0, -1, 0), int3(0, 0, -1)), int3(-2, 0, 1), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(0, -1, 0), int3(0, 0, -1)), int3(1, 0, 0), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 1), int3(0, -1, 0), int3(0, 0, -1)), int3(2, 0, 1), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 1, -1), int3(0, -1, 0), int3(0, 0, -1)), int3(-2, -1, 1), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 1, 0), int3(0, -1, 0), int3(0, 0, -1)), int3(2, 1, 0), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 1, 1), int3(0, -1, 0), int3(0, 0, -1)), int3(2, 1, 1), int3(1, 0, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(0, -1, 0), int3(1, 0, -1)), int3(1, 0, 0), int3(2, 0, 1)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(1, -1, 0), int3(-1, 0, -1)), int3(1, 0, 0), int3(-2, -1, 1)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(1, -1, 0), int3(0, 0, -1)), int3(1, 0, 0), int3(2, 1, 0)),
    std::make_tuple(SKRotationMatrix(int3(1, 0, 0), int3(1, -1, 0), int3(1, 0, -1)), int3(1, 0, 0), int3(2, 1, 1))};

// all possible rotation axes written in terms of integers
std::vector<int3> SKRotationMatrix::allPossibleRotationAxes = {
    int3(1, 0, 0),   int3(0, 1, 0),   int3(0, 0, 1),   int3(0, 1, 1),   int3(1, 0, 1),   int3(1, 1, 0),
    int3(0, 1, -1),  int3(-1, 0, 1),  int3(1, -1, 0),  int3(1, 1, 1),   int3(-1, 1, 1),  int3(1, -1, 1),
    int3(1, 1, -1),  int3(0, 1, 2),   int3(2, 0, 1),   int3(1, 2, 0),   int3(0, 2, 1),   int3(1, 0, 2),
    int3(2, 1, 0),   int3(0, -1, 2),  int3(2, 0, -1),  int3(-1, 2, 0),  int3(0, -2, 1),  int3(1, 0, -2),
    int3(-2, 1, 0),  int3(2, 1, 1),   int3(1, 2, 1),   int3(1, 1, 2),   int3(2, -1, -1), int3(-1, 2, -1),
    int3(-1, -1, 2), int3(2, 1, -1),  int3(-1, 2, 1),  int3(1, -1, 2),  int3(2, -1, 1),  int3(1, 2, -1),
    int3(-1, 1, 2),  int3(3, 1, 2),   int3(2, 3, 1),   int3(1, 2, 3),   int3(3, 2, 1),   int3(1, 3, 2),
    int3(2, 1, 3),   int3(3, -1, 2),  int3(2, 3, -1),  int3(-1, 2, 3),  int3(3, -2, 1),  int3(1, 3, -2),
    int3(-2, 1, 3),  int3(3, -1, -2), int3(-2, 3, -1), int3(-1, -2, 3), int3(3, -2, -1), int3(-1, 3, -2),
    int3(-2, -1, 3), int3(3, 1, -2),  int3(-2, 3, 1),  int3(1, -2, 3),  int3(3, 2, -1),  int3(-1, 3, 2),
    int3(2, -1, 3),  int3(1, 1, 3),   int3(-1, 1, 3),  int3(1, -1, 3),  int3(-1, -1, 3), int3(1, 3, 1),
    int3(-1, 3, 1),  int3(1, 3, -1),  int3(-1, 3, -1), int3(3, 1, 1),   int3(3, 1, -1),  int3(3, -1, 1),
    int3(3, -1, -1)};
