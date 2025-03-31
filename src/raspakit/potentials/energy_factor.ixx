module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#endif

export module energy_factor;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <fstream>;
#endif

import archive;

export namespace Potentials
{
/**
 * \brief Represents an energy factor with energy and its derivative with respect to lambda.
 *
 * The EnergyFactor struct holds the energy value and its derivative with respect to the scaling parameter lambda.
 * It provides constructors and operator overloads for arithmetic operations.
 */
struct EnergyFactor
{
  double energy;     ///< The energy value.
  double dUdlambda;  ///< The derivative of energy with respect to lambda.

  /**
   * \brief Constructs an EnergyFactor with specified energy and derivative.
   *
   * \param energy The energy value.
   * \param dUdlambda The derivative of energy with respect to lambda.
   */
  EnergyFactor(double energy, double dUdlambda) : energy(energy), dUdlambda(dUdlambda) {}

  bool operator==(EnergyFactor const&) const = default;

  inline EnergyFactor& operator+=(const EnergyFactor& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    return *this;
  }

  inline EnergyFactor& operator-=(const EnergyFactor& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    return *this;
  }

  inline EnergyFactor operator-() const
  {
    EnergyFactor v(0.0, 0.0);
    v.energy = -energy;
    v.dUdlambda = -dUdlambda;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Potentials::EnergyFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Potentials::EnergyFactor& e);
};

inline EnergyFactor operator+(const EnergyFactor& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;

  return m;
}

inline EnergyFactor operator-(const EnergyFactor& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;

  return m;
}

inline EnergyFactor operator*(const EnergyFactor& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;

  return m;
}

inline EnergyFactor operator*(const double& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;

  return m;
}

inline EnergyFactor operator*(const EnergyFactor& a, const double& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;

  return m;
}

inline EnergyFactor operator/(const EnergyFactor& a, const double& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;

  return m;
}

inline EnergyFactor sqrt(const EnergyFactor& a)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);

  return m;
}
}  // namespace Potentials
