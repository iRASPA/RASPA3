module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cmath>
#include <fstream>
#endif

export module sixth_derivative_factor;

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
export struct SixthDerivativeFactor
{
  double energy;                  ///< The energy component.
  double dUdlambda;               ///< The derivative of the potential energy with respect to lambda.
  double firstDerivativeFactor;   ///< The scaling factor for the first derivative.
  double secondDerivativeFactor;  ///< The scaling factor for the second derivative.
  double thirdDerivativeFactor;   ///< The scaling factor for the third derivative.
  double fourthDerivativeFactor;  ///< The scaling factor for the fourth derivative.
  double fifthDerivativeFactor;   ///< The scaling factor for the fifth derivative.
  double sixthDerivativeFactor;   ///< The scaling factor for the sixth derivative.

  /**
   * \brief Constructs a ThirdDerivativeFactor with specified energy, force factor, and dUdlambda.
   *
   * \param energy The energy value.
   * \param dUdlambda The derivative of the potential energy with respect to lambda.
   * \param firstDerivativeFactor The first-derivative scaling factor.
   * \param secondDerivativeFactor The second-derivative scaling factor.
   * \param thirdDerivativeFactor The third-derivative scaling factor.
   */
  SixthDerivativeFactor(double energy, double dUdlambda, double firstDerivativeFactor, double secondDerivativeFactor, double thirdDerivativeFactor,
                        double fourthDerivativeFactor, double fifthDerivativeFactor, double sixthDerivativeFactor):
        energy(energy), 
        dUdlambda(dUdlambda), 
        firstDerivativeFactor(firstDerivativeFactor), 
        secondDerivativeFactor(secondDerivativeFactor),
        thirdDerivativeFactor(thirdDerivativeFactor),
        fourthDerivativeFactor(fourthDerivativeFactor),
        fifthDerivativeFactor(fifthDerivativeFactor),
        sixthDerivativeFactor(sixthDerivativeFactor)
  {
  }

  bool operator==(SixthDerivativeFactor const&) const = default;

  inline SixthDerivativeFactor& operator+=(const SixthDerivativeFactor& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    firstDerivativeFactor += b.firstDerivativeFactor;
    secondDerivativeFactor += b.secondDerivativeFactor;
    thirdDerivativeFactor += b.thirdDerivativeFactor;
    fourthDerivativeFactor += b.fourthDerivativeFactor;
    fifthDerivativeFactor += b.fifthDerivativeFactor;
    sixthDerivativeFactor += b.sixthDerivativeFactor;

    return *this;
  }

  inline SixthDerivativeFactor& operator-=(const SixthDerivativeFactor& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    firstDerivativeFactor -= b.firstDerivativeFactor;
    secondDerivativeFactor -= b.secondDerivativeFactor;
    thirdDerivativeFactor -= b.thirdDerivativeFactor;
    fourthDerivativeFactor -= b.fourthDerivativeFactor;
    fifthDerivativeFactor -= b.fifthDerivativeFactor;
    sixthDerivativeFactor -= b.sixthDerivativeFactor;

    return *this;
  }

  inline SixthDerivativeFactor operator-() const
  {
    SixthDerivativeFactor v(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    v.energy = -energy;
    v.dUdlambda = -dUdlambda;
    v.firstDerivativeFactor = -firstDerivativeFactor;
    v.secondDerivativeFactor = -secondDerivativeFactor;
    v.thirdDerivativeFactor = -thirdDerivativeFactor;
    v.fourthDerivativeFactor = -fourthDerivativeFactor;
    v.fifthDerivativeFactor = -fifthDerivativeFactor;
    v.sixthDerivativeFactor = -sixthDerivativeFactor;

    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const SixthDerivativeFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, SixthDerivativeFactor& e);
};


export inline SixthDerivativeFactor operator+(const SixthDerivativeFactor& a, const SixthDerivativeFactor& b)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor + b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor + b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor + b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor + b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor + b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor + b.sixthDerivativeFactor;

  return m;
}

export inline SixthDerivativeFactor operator-(const SixthDerivativeFactor& a, const SixthDerivativeFactor& b)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor - b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor - b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor - b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor - b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor - b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor - b.sixthDerivativeFactor;

  return m;
}

export inline SixthDerivativeFactor operator*(const SixthDerivativeFactor& a, const SixthDerivativeFactor& b)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor * b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor * b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor * b.sixthDerivativeFactor;

  return m;
}

export inline SixthDerivativeFactor operator*(const double& a, const SixthDerivativeFactor& b)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;
  m.firstDerivativeFactor = a * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a * b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a * b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a * b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a * b.sixthDerivativeFactor;

  return m;
}

export inline SixthDerivativeFactor operator*(const SixthDerivativeFactor& a, const double& b)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor * b;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor * b;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor * b;

  return m;
}

export inline SixthDerivativeFactor operator/(const SixthDerivativeFactor& a, const double& b)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;
  m.firstDerivativeFactor = a.firstDerivativeFactor / b;
  m.secondDerivativeFactor = a.secondDerivativeFactor / b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor / b;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor / b;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor / b;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor / b;

  return m;
}

export inline SixthDerivativeFactor sqrt(const SixthDerivativeFactor& a)
{
  SixthDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);
  m.firstDerivativeFactor = std::sqrt(a.firstDerivativeFactor);
  m.secondDerivativeFactor = std::sqrt(a.secondDerivativeFactor);
  m.thirdDerivativeFactor = std::sqrt(a.thirdDerivativeFactor);
  m.fourthDerivativeFactor = std::sqrt(a.fourthDerivativeFactor);
  m.fifthDerivativeFactor = std::sqrt(a.fifthDerivativeFactor);
  m.sixthDerivativeFactor = std::sqrt(a.sixthDerivativeFactor);

  return m;
}
