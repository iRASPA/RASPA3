module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <fstream>
#endif

export module force_factor;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <fstream>;
#endif

import archive;

/**
 * \brief Represents the force factors associated with an energy component.
 *
 * The ForceFactor struct encapsulates the energy, force scaling factor, and the derivative
 * of the potential energy with respect to lambda. It provides constructors for initializing
 * these values and overloaded operators for arithmetic operations. Scaling is linear,
 * first switching Lennard-Jones (LJ) interactions on in the range [0, 0.5], then
 * electrostatic interactions from [0.5, 1.0].
 */
export struct ForceFactor
{
  double energy;       ///< The energy component.
  double forceFactor;  ///< The scaling factor for forces.
  double dUdlambda;    ///< The derivative of the potential energy with respect to lambda.

  /**
   * \brief Constructs a ForceFactor with specified energy, force factor, and dUdlambda.
   *
   * \param energy The energy value.
   * \param forceFactor The force scaling factor.
   * \param dUdlambda The derivative of the potential energy with respect to lambda.
   */
  ForceFactor(double energy, double forceFactor, double dUdlambda)
      : energy(energy), forceFactor(forceFactor), dUdlambda(dUdlambda)
  {
  }

  bool operator==(ForceFactor const&) const = default;

  inline ForceFactor& operator+=(const ForceFactor& b)
  {
    energy += b.energy;
    forceFactor += b.forceFactor;
    dUdlambda += b.dUdlambda;
    return *this;
  }

  inline ForceFactor& operator-=(const ForceFactor& b)
  {
    energy -= b.energy;
    forceFactor -= b.forceFactor;
    dUdlambda -= b.dUdlambda;
    return *this;
  }

  inline ForceFactor operator-() const
  {
    ForceFactor v(0.0, 0.0, 0.0);
    v.energy = -energy;
    v.forceFactor = -forceFactor;
    v.dUdlambda = -dUdlambda;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ForceFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ForceFactor& e);
};

// scaling is linear and first switch LJ on in 0-0.5, then the electrostatics from 0.5 to 1.0
export inline ForceFactor operator+(const ForceFactor& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.forceFactor = a.forceFactor + b.forceFactor;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;

  return m;
}

export inline ForceFactor operator-(const ForceFactor& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.forceFactor = a.forceFactor - b.forceFactor;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;

  return m;
}

export inline ForceFactor operator*(const ForceFactor& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.forceFactor = a.forceFactor * b.forceFactor;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;

  return m;
}

export inline ForceFactor operator*(const double& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.forceFactor = a * b.forceFactor;
  m.dUdlambda = a * b.dUdlambda;

  return m;
}

export inline ForceFactor operator*(const ForceFactor& a, const double& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.forceFactor = a.forceFactor * b;
  m.dUdlambda = a.dUdlambda * b;

  return m;
}

export inline ForceFactor operator/(const ForceFactor& a, const double& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.forceFactor = a.forceFactor / b;
  m.dUdlambda = a.dUdlambda / b;

  return m;
}

export inline ForceFactor sqrt(const ForceFactor& a)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.forceFactor = std::sqrt(a.forceFactor);
  m.dUdlambda = std::sqrt(a.dUdlambda);

  return m;
}
