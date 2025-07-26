module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#endif

export module tricubic_derivative_factor;

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
 * The TricubicDerivativeFactor struct encapsulates the energy, force scaling factor, and the derivative
 * of the potential energy with respect to lambda. It provides constructors for initializing
 * these values and overloaded operators for arithmetic operations. Scaling is linear,
 * first switching Lennard-Jones (LJ) interactions on in the range [0, 0.5], then
 * electrostatic interactions from [0.5, 1.0].
 */
struct TricubicDerivativeFactor
{
  double energy;                  ///< The energy component.
  double firstDerivativeFactor;   ///< The scaling factor for the first derivative.
  double secondDerivativeFactor;  ///< The scaling factor for the second derivative.
  double thirdDerivativeFactor;   ///< The scaling factor for the third derivative.

  /**
   * \brief Constructs a TricubicDerivativeFactor with specified energy, force factor.
   *
   * \param energy The energy value.
   * \param firstDerivativeFactor The first-derivative scaling factor.
   * \param secondDerivativeFactor The second-derivative scaling factor.
   * \param thirdDerivativeFactor The third-derivative scaling factor.
   */
  TricubicDerivativeFactor(double energy, double firstDerivativeFactor, double secondDerivativeFactor,
                           double thirdDerivativeFactor)
      : energy(energy),
        firstDerivativeFactor(firstDerivativeFactor),
        secondDerivativeFactor(secondDerivativeFactor),
        thirdDerivativeFactor(thirdDerivativeFactor)
  {
  }

  TricubicDerivativeFactor()
      : energy(0.0), firstDerivativeFactor(0.0), secondDerivativeFactor(0.0), thirdDerivativeFactor(0.0)
  {
  }

  bool operator==(TricubicDerivativeFactor const&) const = default;

  inline TricubicDerivativeFactor& operator+=(const TricubicDerivativeFactor& b)
  {
    energy += b.energy;
    firstDerivativeFactor += b.firstDerivativeFactor;
    secondDerivativeFactor += b.secondDerivativeFactor;
    thirdDerivativeFactor += b.thirdDerivativeFactor;
    return *this;
  }

  inline TricubicDerivativeFactor& operator-=(const TricubicDerivativeFactor& b)
  {
    energy -= b.energy;
    firstDerivativeFactor -= b.firstDerivativeFactor;
    secondDerivativeFactor -= b.secondDerivativeFactor;
    thirdDerivativeFactor -= b.thirdDerivativeFactor;
    return *this;
  }

  inline TricubicDerivativeFactor operator-() const
  {
    TricubicDerivativeFactor v(0.0, 0.0, 0.0, 0.0);
    v.energy = -energy;
    v.firstDerivativeFactor = -firstDerivativeFactor;
    v.secondDerivativeFactor = -secondDerivativeFactor;
    v.thirdDerivativeFactor = -thirdDerivativeFactor;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const TricubicDerivativeFactor& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, TricubicDerivativeFactor& e);
};

inline TricubicDerivativeFactor operator+(const TricubicDerivativeFactor& a, const TricubicDerivativeFactor& b)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.firstDerivativeFactor = a.firstDerivativeFactor + b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor + b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor + b.thirdDerivativeFactor;

  return m;
}

inline TricubicDerivativeFactor operator-(const TricubicDerivativeFactor& a, const TricubicDerivativeFactor& b)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.firstDerivativeFactor = a.firstDerivativeFactor - b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor - b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor - b.thirdDerivativeFactor;

  return m;
}

inline TricubicDerivativeFactor operator*(const TricubicDerivativeFactor& a, const TricubicDerivativeFactor& b)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b.thirdDerivativeFactor;

  return m;
}

inline TricubicDerivativeFactor operator*(const double& a, const TricubicDerivativeFactor& b)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.firstDerivativeFactor = a * b.firstDerivativeFactor;
  m.secondDerivativeFactor = a * b.secondDerivativeFactor;
  m.thirdDerivativeFactor = a * b.thirdDerivativeFactor;

  return m;
}

inline TricubicDerivativeFactor operator*(const TricubicDerivativeFactor& a, const double& b)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.firstDerivativeFactor = a.firstDerivativeFactor * b;
  m.secondDerivativeFactor = a.secondDerivativeFactor * b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor * b;

  return m;
}

inline TricubicDerivativeFactor operator/(const TricubicDerivativeFactor& a, const double& b)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.firstDerivativeFactor = a.firstDerivativeFactor / b;
  m.secondDerivativeFactor = a.secondDerivativeFactor / b;
  m.thirdDerivativeFactor = a.thirdDerivativeFactor / b;

  return m;
}

inline TricubicDerivativeFactor sqrt(const TricubicDerivativeFactor& a)
{
  TricubicDerivativeFactor m(0.0, 0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.firstDerivativeFactor = std::sqrt(a.firstDerivativeFactor);
  m.secondDerivativeFactor = std::sqrt(a.secondDerivativeFactor);
  m.thirdDerivativeFactor = std::sqrt(a.thirdDerivativeFactor);

  return m;
}
}  // namespace Potentials
