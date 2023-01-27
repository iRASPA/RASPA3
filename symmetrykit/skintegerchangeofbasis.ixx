export module skintegerchangeofbasis;

import int3;
import int3x3;
import double3;
import double3x3;

import sktransformationmatrix;

export class SKIntegerChangeOfBasis
{
public:
    SKIntegerChangeOfBasis() { ; }
    SKIntegerChangeOfBasis(SKTransformationMatrix inversionTransformation);
    SKIntegerChangeOfBasis inverse();
    SKTransformationMatrix operator * (const SKTransformationMatrix& b) const;
    double3 operator * (const double3& b) const;
    int3 operator * (const int3& b) const;
private:
    SKTransformationMatrix _changeOfBasis;
    SKTransformationMatrix _inverseChangeOfBasis;
    int _changeOfBasisDeterminant = 1;
    int _inverseChangeOfBasisDeterminant = 1;
};
