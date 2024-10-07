module;

module skintegerchangeofbasis;

import int3;
import int3x3;
import double3;
import double3x3;

import sktransformationmatrix;

SKIntegerChangeOfBasis::SKIntegerChangeOfBasis(SKTransformationMatrix inversionTransformation)
{
  _inverseChangeOfBasis = inversionTransformation;
  _changeOfBasis = inversionTransformation.adjugate();
  _changeOfBasisDeterminant = inversionTransformation.determinant();
  _inverseChangeOfBasisDeterminant = 1;
}

SKIntegerChangeOfBasis SKIntegerChangeOfBasis::inverse()
{
  SKIntegerChangeOfBasis basis{};
  basis._changeOfBasis = _inverseChangeOfBasis;
  basis._changeOfBasisDeterminant = _inverseChangeOfBasisDeterminant;
  basis._inverseChangeOfBasis = _changeOfBasis;
  basis._inverseChangeOfBasisDeterminant = _changeOfBasisDeterminant;
  return basis;
}

SKTransformationMatrix SKIntegerChangeOfBasis::operator*(const SKTransformationMatrix& right) const
{
  return this->_inverseChangeOfBasis * right * this->_changeOfBasis /
         int(this->_inverseChangeOfBasisDeterminant * this->_changeOfBasisDeterminant);
}

double3 SKIntegerChangeOfBasis::operator*(const double3& right) const
{
  return this->_inverseChangeOfBasis * right / double(this->_inverseChangeOfBasisDeterminant);
}

int3 SKIntegerChangeOfBasis::operator*(const int3& right) const
{
  return this->_inverseChangeOfBasis * right / this->_inverseChangeOfBasisDeterminant;
}
