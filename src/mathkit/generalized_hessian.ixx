module;

export module generalized_hessian;

import std;

import double3x3;
import entity_hessian_layout;

/**
 * Generalized Hessian for minimization: position-position, position-strain, and strain-strain blocks.
 */
export class GeneralizedHessian
{
 public:
  GeneralizedHessian() = default;

  GeneralizedHessian(std::size_t nDofs, std::size_t nStrain);

  std::size_t numDofs() const noexcept { return _nDofs; }
  std::size_t numStrain() const noexcept { return _nStrain; }

  void resize(std::size_t nDofs, std::size_t nStrain);

  void zero();

  void add(std::size_t i, std::size_t j, double value);
  void addPositionStrain(std::size_t i, std::size_t strainIndex, double value);
  void addStrainStrain(std::size_t i, std::size_t j, double value);

  std::span<double> positionPosition() noexcept { return _positionPosition; }
  std::span<const double> positionPosition() const noexcept { return _positionPosition; }

  std::span<double> positionStrain() noexcept { return _positionStrain; }
  std::span<const double> positionStrain() const noexcept { return _positionStrain; }

  std::span<double> strainStrain() noexcept { return _strainStrain; }
  std::span<const double> strainStrain() const noexcept { return _strainStrain; }

  double3x3 &strainGradient() noexcept { return _strainGradient; }
  const double3x3 &strainGradient() const noexcept { return _strainGradient; }

  double &operator()(std::size_t i, std::size_t j) { return _positionPosition[i * _nDofs + j]; }
  double operator()(std::size_t i, std::size_t j) const { return _positionPosition[i * _nDofs + j]; }

 private:
  std::size_t _nDofs{};
  std::size_t _nStrain{};
  std::vector<double> _positionPosition;
  std::vector<double> _positionStrain;
  std::vector<double> _strainStrain;
  double3x3 _strainGradient{};
};
