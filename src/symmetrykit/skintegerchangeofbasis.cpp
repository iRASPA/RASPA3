module;

module skintegerchangeofbasis;

import std;

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
  // a singular transformation is no change of basis: every product with it divides by its determinant
  if (_changeOfBasisDeterminant == 0)
  {
    throw std::invalid_argument("SKIntegerChangeOfBasis: the transformation is singular (determinant zero)");
  }
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
         (this->_inverseChangeOfBasisDeterminant * this->_changeOfBasisDeterminant);
}

double3 SKIntegerChangeOfBasis::operator*(const double3& right) const
{
  return this->_inverseChangeOfBasis * right / double(this->_inverseChangeOfBasisDeterminant);
}

int3 SKIntegerChangeOfBasis::operator*(const int3& right) const
{
  return this->_inverseChangeOfBasis * right / this->_inverseChangeOfBasisDeterminant;
}
