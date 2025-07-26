module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#endif

export module hessian_factor;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

export namespace Potentials
{
/**
 * \brief Represents the first and second derivative factors associated with an energy component.
 *
 * The HessianFactor struct encapsulates the energy, force scaling factor, and the derivative
 * of the potential energy with respect to lambda. It provides constructors for initializing
 * these values and overloaded operators for arithmetic operations. Scaling is linear,
 * first switching Lennard-Jones (LJ) interactions on in the range [0, 0.5], then
 * electrostatic interactions from [0.5, 1.0].
 */
struct HessianFactor
{
  double energy;                  ///< The energy component.
  double dUdlambda;               ///< The derivative of the potential energy with respect to lambda.
  double firstDerivativeFactor;   ///< The scaling factor for the first derivative.
  double secondDerivativeFactor;  ///< The scaling factor for the second derivative.

  /**
   * \brief Constructs a HessianFactor with specified energy, force factor, and dUdlambda.
   *
   * \param energy The energy value.
   * \param dUdlambda The derivative of the potential energy with respect to lambda.
   * \param firstDerivativeFactor The first-derivative scaling factor.
   * \param secondDerivativeFactor The second-derivative scaling factor.
   */
  HessianFactor(double energy, double dUdlambda, double firstDerivativeFactor, double secondDerivativeFactor)
      : energy(energy),
        dUdlambda(dUdlambda),
        firstDerivativeFactor(firstDerivativeFactor),
        secondDerivativeFactor(secondDerivativeFactor)
  {
  }

  bool operator==(HessianFactor const&) const = default;

  inline HessianFactor& operator+=(const HessianFactor& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    firstDerivativeFactor += b.firstDerivativeFactor;
    secondDerivativeFactor += b.secondDerivativeFactor;
    return *this;
  }

  inline HessianFactor& operator-=(const HessianFactor& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    firstDerivativeFactor -= b.firstDerivativeFactor;
    secondDerivativeFactor -= b.secondDerivativeFactor;
    return *this;
  }

  inline HessianFactor operator-() const
  {
    HessianFactor v(0.0, 0.0, 0.0, 0.0);
    v.energy = -energy;
    v.dUdlambda = -dUdlambda;
    v.firstDerivativeFactor = -firstDerivativeFactor;
    v.secondDerivativeFactor = -secondDerivativeFactor;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const HessianFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, HessianFactor& e);
};

inline HessianFactor operator+(const HessianFactor& a, const HessianFactor& b)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor + b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor + b.secondDerivativeFactor;

  return m;
}

inline HessianFactor operator-(const HessianFactor& a, const HessianFactor& b)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor - b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor - b.secondDerivativeFactor;

  return m;
}

inline HessianFactor operator*(const HessianFactor& a, const HessianFactor& b)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b.secondDerivativeFactor;

  return m;
}

inline HessianFactor operator*(const double& a, const HessianFactor& b)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;
  m.firstDerivativeFactor = a * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a * b.secondDerivativeFactor;

  return m;
}

inline HessianFactor operator*(const HessianFactor& a, const double& b)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b;

  return m;
}

inline HessianFactor operator/(const HessianFactor& a, const double& b)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;
  m.firstDerivativeFactor = a.firstDerivativeFactor / b;
  m.secondDerivativeFactor = a.secondDerivativeFactor / b;

  return m;
}

inline HessianFactor sqrt(const HessianFactor& a)
{
  HessianFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);
  m.firstDerivativeFactor = std::sqrt(a.firstDerivativeFactor);
  m.secondDerivativeFactor = std::sqrt(a.secondDerivativeFactor);

  return m;
}
}  // namespace Potentials
