module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#endif

export module triquintic_derivative_factor;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

export namespace Potentials
{
/**
 * \brief Represents the force factors associated with an energy component.
 *
 * The ThirdDerivativeFactor struct encapsulates the energy, force scaling factor, and the derivative
 * of the potential energy with respect to lambda. It provides constructors for initializing
 * these values and overloaded operators for arithmetic operations. Scaling is linear,
 * first switching Lennard-Jones (LJ) interactions on in the range [0, 0.5], then
 * electrostatic interactions from [0.5, 1.0].
 */
struct TriquinticDerivativeFactor
{
  double energy;                  ///< The energy component.
  double firstDerivativeFactor;   ///< The scaling factor for the first derivative.
  double secondDerivativeFactor;  ///< The scaling factor for the second derivative.
  double thirdDerivativeFactor;   ///< The scaling factor for the third derivative.
  double fourthDerivativeFactor;  ///< The scaling factor for the fourth derivative.
  double fifthDerivativeFactor;   ///< The scaling factor for the fifth derivative.
  double sixthDerivativeFactor;   ///< The scaling factor for the sixth derivative.

  /**
   * \brief Constructs a ThirdDerivativeFactor with specified energy, force factor.
   *
   * \param energy The energy value.
   * \param firstDerivativeFactor The first-derivative scaling factor.
   * \param secondDerivativeFactor The second-derivative scaling factor.
   * \param thirdDerivativeFactor The third-derivative scaling factor.
   */
  TriquinticDerivativeFactor(double energy, double firstDerivativeFactor, double secondDerivativeFactor,
                             double thirdDerivativeFactor, double fourthDerivativeFactor, double fifthDerivativeFactor,
                             double sixthDerivativeFactor)
      : energy(energy),
        firstDerivativeFactor(firstDerivativeFactor),
        secondDerivativeFactor(secondDerivativeFactor),
        thirdDerivativeFactor(thirdDerivativeFactor),
        fourthDerivativeFactor(fourthDerivativeFactor),
        fifthDerivativeFactor(fifthDerivativeFactor),
        sixthDerivativeFactor(sixthDerivativeFactor)
  {
  }

  TriquinticDerivativeFactor()
      : energy(0.0),
        firstDerivativeFactor(0.0),
        secondDerivativeFactor(0.0),
        thirdDerivativeFactor(0.0),
        fourthDerivativeFactor(0.0),
        fifthDerivativeFactor(0.0),
        sixthDerivativeFactor(0.0)
  {
  }

  bool operator==(TriquinticDerivativeFactor const&) const = default;

  inline TriquinticDerivativeFactor& operator+=(const TriquinticDerivativeFactor& b)
  {
    energy += b.energy;
    firstDerivativeFactor += b.firstDerivativeFactor;
    secondDerivativeFactor += b.secondDerivativeFactor;
    thirdDerivativeFactor += b.thirdDerivativeFactor;
    fourthDerivativeFactor += b.fourthDerivativeFactor;
    fifthDerivativeFactor += b.fifthDerivativeFactor;
    sixthDerivativeFactor += b.sixthDerivativeFactor;

    return *this;
  }

  inline TriquinticDerivativeFactor& operator-=(const TriquinticDerivativeFactor& b)
  {
    energy -= b.energy;
    firstDerivativeFactor -= b.firstDerivativeFactor;
    secondDerivativeFactor -= b.secondDerivativeFactor;
    thirdDerivativeFactor -= b.thirdDerivativeFactor;
    fourthDerivativeFactor -= b.fourthDerivativeFactor;
    fifthDerivativeFactor -= b.fifthDerivativeFactor;
    sixthDerivativeFactor -= b.sixthDerivativeFactor;

    return *this;
  }

  inline TriquinticDerivativeFactor operator-() const
  {
    TriquinticDerivativeFactor v(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    v.energy = -energy;
    v.firstDerivativeFactor = -firstDerivativeFactor;
    v.secondDerivativeFactor = -secondDerivativeFactor;
    v.thirdDerivativeFactor = -thirdDerivativeFactor;
    v.fourthDerivativeFactor = -fourthDerivativeFactor;
    v.fifthDerivativeFactor = -fifthDerivativeFactor;
    v.sixthDerivativeFactor = -sixthDerivativeFactor;

    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const TriquinticDerivativeFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, TriquinticDerivativeFactor& e);
};

inline TriquinticDerivativeFactor operator+(const TriquinticDerivativeFactor& a, const TriquinticDerivativeFactor& b)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.firstDerivativeFactor = a.firstDerivativeFactor + b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor + b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor + b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor + b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor + b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor + b.sixthDerivativeFactor;

  return m;
}

inline TriquinticDerivativeFactor operator-(const TriquinticDerivativeFactor& a, const TriquinticDerivativeFactor& b)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.firstDerivativeFactor = a.firstDerivativeFactor - b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor - b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor - b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor - b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor - b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor - b.sixthDerivativeFactor;

  return m;
}

inline TriquinticDerivativeFactor operator*(const TriquinticDerivativeFactor& a, const TriquinticDerivativeFactor& b)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor * b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor * b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor * b.sixthDerivativeFactor;

  return m;
}

inline TriquinticDerivativeFactor operator*(const double& a, const TriquinticDerivativeFactor& b)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.firstDerivativeFactor = a * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a * b.thirdDerivativeFactor;
  m.fourthDerivativeFactor = a * b.fourthDerivativeFactor;
  m.fifthDerivativeFactor = a * b.fifthDerivativeFactor;
  m.sixthDerivativeFactor = a * b.sixthDerivativeFactor;

  return m;
}

inline TriquinticDerivativeFactor operator*(const TriquinticDerivativeFactor& a, const double& b)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor * b;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor * b;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor * b;

  return m;
}

inline TriquinticDerivativeFactor operator/(const TriquinticDerivativeFactor& a, const double& b)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.firstDerivativeFactor = a.firstDerivativeFactor / b;
  m.secondDerivativeFactor = a.secondDerivativeFactor / b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor / b;
  m.fourthDerivativeFactor = a.fourthDerivativeFactor / b;
  m.fifthDerivativeFactor = a.fifthDerivativeFactor / b;
  m.sixthDerivativeFactor = a.sixthDerivativeFactor / b;

  return m;
}

inline TriquinticDerivativeFactor sqrt(const TriquinticDerivativeFactor& a)
{
  TriquinticDerivativeFactor m(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.firstDerivativeFactor = std::sqrt(a.firstDerivativeFactor);
  m.secondDerivativeFactor = std::sqrt(a.secondDerivativeFactor);
  m.thirdDerivativeFactor = std::sqrt(a.thirdDerivativeFactor);
  m.fourthDerivativeFactor = std::sqrt(a.fourthDerivativeFactor);
  m.fifthDerivativeFactor = std::sqrt(a.fifthDerivativeFactor);
  m.sixthDerivativeFactor = std::sqrt(a.sixthDerivativeFactor);

  return m;
}
}  // namespace Potentials
