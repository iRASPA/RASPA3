module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#endif

export module gradient_factor;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <fstream>;
#endif

import archive;

export namespace Potentials
{
/**
 * \brief Represents the force factors associated with an energy component.
 *
 * The GradientFactor struct encapsulates the energy, force scaling factor, and the derivative
 * of the potential energy with respect to lambda. It provides constructors for initializing
 * these values and overloaded operators for arithmetic operations. Scaling is linear,
 * first switching Lennard-Jones (LJ) interactions on in the range [0, 0.5], then
 * electrostatic interactions from [0.5, 1.0].
 */
struct GradientFactor
{
  double energy;          ///< The energy component.
  double dUdlambda;       ///< The derivative of the potential energy with respect to lambda.
  double gradientFactor;  ///< The scaling factor for forces.

  /**
   * \brief Constructs a GradientFactor with specified energy, force factor, and dUdlambda.
   *
   * \param energy The energy value.
   * \param gradientFactor The force scaling factor.
   * \param dUdlambda The derivative of the potential energy with respect to lambda.
   */
  GradientFactor(double energy, double dUdlambda, double gradientFactor)
      : energy(energy), dUdlambda(dUdlambda), gradientFactor(gradientFactor)
  {
  }

  bool operator==(GradientFactor const&) const = default;

  inline GradientFactor& operator+=(const GradientFactor& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    gradientFactor += b.gradientFactor;
    return *this;
  }

  inline GradientFactor& operator-=(const GradientFactor& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    gradientFactor -= b.gradientFactor;
    return *this;
  }

  inline GradientFactor operator-() const
  {
    GradientFactor v(0.0, 0.0, 0.0);
    v.dUdlambda = -dUdlambda;
    v.energy = -energy;
    v.gradientFactor = -gradientFactor;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const GradientFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, GradientFactor& e);
};

// scaling is linear and first switch LJ on in 0-0.5, then the electrostatics from 0.5 to 1.0
inline GradientFactor operator+(const GradientFactor& a, const GradientFactor& b)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.dUdlambda = a.dUdlambda + b.dUdlambda;
  m.energy = a.energy + b.energy;
  m.gradientFactor = a.gradientFactor + b.gradientFactor;

  return m;
}

inline GradientFactor operator-(const GradientFactor& a, const GradientFactor& b)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;
  m.gradientFactor = a.gradientFactor - b.gradientFactor;

  return m;
}

inline GradientFactor operator*(const GradientFactor& a, const GradientFactor& b)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;
  m.gradientFactor = a.gradientFactor * b.gradientFactor;

  return m;
}

inline GradientFactor operator*(const double& a, const GradientFactor& b)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;
  m.gradientFactor = a * b.gradientFactor;

  return m;
}

inline GradientFactor operator*(const GradientFactor& a, const double& b)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;
  m.gradientFactor = a.gradientFactor * b;

  return m;
}

inline GradientFactor operator/(const GradientFactor& a, const double& b)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;
  m.gradientFactor = a.gradientFactor / b;

  return m;
}

inline GradientFactor sqrt(const GradientFactor& a)
{
  GradientFactor m(0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);
  m.gradientFactor = std::sqrt(a.gradientFactor);

  return m;
}
}  // namespace Potentials
