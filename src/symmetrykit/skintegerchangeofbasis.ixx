export module skintegerchangeofbasis;

import int3;
import int3x3;
import double3;
import double3x3;

import sktransformationmatrix;

export class SKIntegerChangeOfBasis
{
 public:
  // the identity (the matrices used to be left zero beside determinants of one)
  SKIntegerChangeOfBasis()
      : _changeOfBasis(int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1)),
        _inverseChangeOfBasis(int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1))
  {
  }
  SKIntegerChangeOfBasis(SKTransformationMatrix inversionTransformation);
  SKIntegerChangeOfBasis inverse();
  SKTransformationMatrix operator*(const SKTransformationMatrix& b) const;
  double3 operator*(const double3& b) const;
  int3 operator*(const int3& b) const;

 private:
  SKTransformationMatrix _changeOfBasis;
  SKTransformationMatrix _inverseChangeOfBasis;
  int _changeOfBasisDeterminant = 1;
  int _inverseChangeOfBasisDeterminant = 1;
};
