module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <fstream>
#endif

export module third_derivative_factor;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <fstream>;
#endif

import archive;

/**
 * \brief Represents the force factors associated with an energy component.
 *
 * The ThirdDerivativeFactor struct encapsulates the energy, force scaling factor, and the derivative
 * of the potential energy with respect to lambda. It provides constructors for initializing
 * these values and overloaded operators for arithmetic operations. Scaling is linear,
 * first switching Lennard-Jones (LJ) interactions on in the range [0, 0.5], then
 * electrostatic interactions from [0.5, 1.0].
 */
export struct ThirdDerivativeFactor
{
  double energy;                  ///< The energy component.
  double dUdlambda;               ///< The derivative of the potential energy with respect to lambda.
  double firstDerivativeFactor;   ///< The scaling factor for the first derivative.
  double secondDerivativeFactor;  ///< The scaling factor for the second derivative.
  double thirdDerivativeFactor;   ///< The scaling factor for the third derivative.

  /**
   * \brief Constructs a ThirdDerivativeFactor with specified energy, force factor, and dUdlambda.
   *
   * \param energy The energy value.
   * \param dUdlambda The derivative of the potential energy with respect to lambda.
   * \param firstDerivativeFactor The first-derivative scaling factor.
   * \param secondDerivativeFactor The second-derivative scaling factor.
   * \param thirdDerivativeFactor The third-derivative scaling factor.
   */
  ThirdDerivativeFactor(double energy, double dUdlambda, double firstDerivativeFactor, double secondDerivativeFactor, double thirdDerivativeFactor):
        energy(energy), 
        dUdlambda(dUdlambda), 
        firstDerivativeFactor(firstDerivativeFactor), 
        secondDerivativeFactor(secondDerivativeFactor),
        thirdDerivativeFactor(thirdDerivativeFactor)
  {
  }

  bool operator==(ThirdDerivativeFactor const&) const = default;

  inline ThirdDerivativeFactor& operator+=(const ThirdDerivativeFactor& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    firstDerivativeFactor += b.firstDerivativeFactor;
    secondDerivativeFactor += b.secondDerivativeFactor;
    thirdDerivativeFactor += b.thirdDerivativeFactor;
    return *this;
  }

  inline ThirdDerivativeFactor& operator-=(const ThirdDerivativeFactor& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    firstDerivativeFactor -= b.firstDerivativeFactor;
    secondDerivativeFactor -= b.secondDerivativeFactor;
    thirdDerivativeFactor -= b.thirdDerivativeFactor;
    return *this;
  }

  inline ThirdDerivativeFactor operator-() const
  {
    ThirdDerivativeFactor v(0.0, 0.0, 0.0, 0.0, 0.0);
    v.energy = -energy;
    v.dUdlambda = -dUdlambda;
    v.firstDerivativeFactor = -firstDerivativeFactor;
    v.secondDerivativeFactor = -secondDerivativeFactor;
    v.thirdDerivativeFactor = -thirdDerivativeFactor;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ThirdDerivativeFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ThirdDerivativeFactor& e);
};


export inline ThirdDerivativeFactor operator+(const ThirdDerivativeFactor& a, const ThirdDerivativeFactor& b)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor + b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor + b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor + b.thirdDerivativeFactor;

  return m;
}

export inline ThirdDerivativeFactor operator-(const ThirdDerivativeFactor& a, const ThirdDerivativeFactor& b)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor - b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor - b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor - b.thirdDerivativeFactor;

  return m;
}

export inline ThirdDerivativeFactor operator*(const ThirdDerivativeFactor& a, const ThirdDerivativeFactor& b)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b.thirdDerivativeFactor;

  return m;
}

export inline ThirdDerivativeFactor operator*(const double& a, const ThirdDerivativeFactor& b)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;
  m.firstDerivativeFactor = a * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a * b.thirdDerivativeFactor;

  return m;
}

export inline ThirdDerivativeFactor operator*(const ThirdDerivativeFactor& a, const double& b)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b;

  return m;
}

export inline ThirdDerivativeFactor operator/(const ThirdDerivativeFactor& a, const double& b)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;
  m.firstDerivativeFactor = a.firstDerivativeFactor / b;
  m.secondDerivativeFactor = a.secondDerivativeFactor / b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor / b;

  return m;
}

export inline ThirdDerivativeFactor sqrt(const ThirdDerivativeFactor& a)
{
  ThirdDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);
  m.firstDerivativeFactor = std::sqrt(a.firstDerivativeFactor);
  m.secondDerivativeFactor = std::sqrt(a.secondDerivativeFactor);
  m.thirdDerivativeFactor = std::sqrt(a.thirdDerivativeFactor);

  return m;
}
