module;

module generalized_hessian;

import std;

import double3x3;

GeneralizedHessian::GeneralizedHessian(std::size_t nDofs, std::size_t nStrain) { resize(nDofs, nStrain); }

void GeneralizedHessian::resize(std::size_t nDofs, std::size_t nStrain)
{
  _nDofs = nDofs;
  _nStrain = nStrain;
  _positionPosition.assign(nDofs * nDofs, 0.0);
  _positionStrain.assign(nDofs * nStrain, 0.0);
  _strainStrain.assign(nStrain * nStrain, 0.0);
  _strainGradient = double3x3{};
}

void GeneralizedHessian::zero()
{
  std::ranges::fill(_positionPosition, 0.0);
  std::ranges::fill(_positionStrain, 0.0);
  std::ranges::fill(_strainStrain, 0.0);
  _strainGradient = double3x3{};
}

void GeneralizedHessian::add(std::size_t i, std::size_t j, double value)
{
  if (i < _nDofs && j < _nDofs)
  {
    _positionPosition[i * _nDofs + j] += value;
  }
}

void GeneralizedHessian::addPositionStrain(std::size_t i, std::size_t strainIndex, double value)
{
  if (i < _nDofs && strainIndex < _nStrain)
  {
    _positionStrain[i * _nStrain + strainIndex] += value;
  }
}

void GeneralizedHessian::addStrainStrain(std::size_t i, std::size_t j, double value)
{
  if (i < _nStrain && j < _nStrain)
  {
    _strainStrain[i * _nStrain + j] += value;
  }
}
