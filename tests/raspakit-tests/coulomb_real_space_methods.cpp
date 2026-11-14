#include <gtest/gtest.h>

import std;

import forcefield;
import potential_coulomb_real_space;

namespace
{
constexpr double derivativeStep = 1.0e-5;

ForceField makeForceField(ForceField::ChargeMethod method)
{
  ForceField forceField;
  forceField.chargeMethod = method;
  forceField.cutOffCoulomb = 14.0;
  forceField.EwaldAlpha = 0.18;
  forceField.modifiedShiftedForceBeta = 0.3;
  return forceField;
}
}  // namespace

TEST(coulomb_real_space_methods, charge_method_names_and_fourier_dispatch)
{
  EXPECT_EQ(ForceField::chargeMethodFromString("Wolf"), ForceField::ChargeMethod::Wolf);
  EXPECT_EQ(ForceField::chargeMethodFromString("DSF"), ForceField::ChargeMethod::DampedShiftedForce);
  EXPECT_EQ(ForceField::chargeMethodFromString("mDSF"), ForceField::ChargeMethod::ModifiedShiftedForce);
  EXPECT_EQ(ForceField::chargeMethodFromString("Zero-Dipole"), ForceField::ChargeMethod::ZeroDipole);

  ForceField forceField;
  forceField.chargeMethod = ForceField::ChargeMethod::Ewald;
  EXPECT_TRUE(forceField.usesEwaldFourier());
  for (const ForceField::ChargeMethod method :
       {ForceField::ChargeMethod::Wolf, ForceField::ChargeMethod::DampedShiftedForce,
        ForceField::ChargeMethod::ModifiedShiftedForce, ForceField::ChargeMethod::ZeroDipole})
  {
    forceField.chargeMethod = method;
    EXPECT_FALSE(forceField.usesEwaldFourier());
  }
}

TEST(coulomb_real_space_methods, potentials_and_forces_vanish_as_defined_at_cutoff)
{
  for (const ForceField::ChargeMethod method :
       {ForceField::ChargeMethod::Wolf, ForceField::ChargeMethod::DampedShiftedForce,
        ForceField::ChargeMethod::ModifiedShiftedForce, ForceField::ChargeMethod::ZeroDipole})
  {
    const ForceField forceField = makeForceField(method);
    const Potentials::CoulombRealSpaceFactors factors =
        Potentials::coulombRealSpaceFactors(forceField, forceField.cutOffCoulomb);
    EXPECT_NEAR(factors.potential, 0.0, 1.0e-14);
    if (method != ForceField::ChargeMethod::Wolf)
    {
      EXPECT_NEAR(factors.firstDerivativeFactor, 0.0, 1.0e-14);
    }
  }
}

TEST(coulomb_real_space_methods, analytical_derivatives_match_finite_differences)
{
  for (const ForceField::ChargeMethod method :
       {ForceField::ChargeMethod::Wolf, ForceField::ChargeMethod::DampedShiftedForce,
        ForceField::ChargeMethod::ModifiedShiftedForce, ForceField::ChargeMethod::ZeroDipole})
  {
    const ForceField forceField = makeForceField(method);
    const double r = 6.7;
    const auto factors = Potentials::coulombRealSpaceFactors(forceField, r);
    const auto minus = Potentials::coulombRealSpaceFactors(forceField, r - derivativeStep);
    const auto plus = Potentials::coulombRealSpaceFactors(forceField, r + derivativeStep);

    const double numericalFirstDerivative = (plus.potential - minus.potential) / (2.0 * derivativeStep);
    const double numericalSecondDerivative =
        (plus.potential - 2.0 * factors.potential + minus.potential) / (derivativeStep * derivativeStep);
    const double analyticalFirstDerivative = factors.firstDerivativeFactor * r;
    const double analyticalSecondDerivative =
        factors.secondDerivativeFactor * r * r + factors.firstDerivativeFactor;

    EXPECT_NEAR(analyticalFirstDerivative, numericalFirstDerivative, 1.0e-10);
    EXPECT_NEAR(analyticalSecondDerivative, numericalSecondDerivative, 2.0e-6);
  }
}

// The Ewald exclusion of a core and its shell (CoreShellSpring) at and near r = 0: the closed forms of erf(alpha s)/s
// and its derivatives cancel catastrophically for small alpha s and are 0/0 at s = 0. The factors are finite there,
// take the exact limits, agree with an independent (longer) series below the switch and with the closed forms above.
TEST(coulomb_real_space_methods, ewald_exclusion_factors_are_finite_and_accurate_at_small_r)
{
  const double alpha = 0.3;
  const double twoOverSqrtPi = 2.0 / std::sqrt(std::numbers::pi);
  const auto finite = [](double x)
  { return std::to_string(x).find("nan") == std::string::npos && std::abs(x) < 1e300; };

  // r = 0, full scaling: erf(alpha r)/r -> 2 alpha/sqrt(pi), phi'/r -> -(4/3) alpha^3/sqrt(pi),
  // (phi'' - phi'/r)/r^2 -> (8/5) alpha^5/sqrt(pi)
  const Potentials::EwaldExclusionFactors atZero = Potentials::ewaldExclusionFactors(alpha, 1.0, 0.0);
  EXPECT_NEAR(atZero.potential, twoOverSqrtPi * alpha, 1e-15);
  EXPECT_NEAR(atZero.firstDerivativeFactor, -(2.0 / 3.0) * twoOverSqrtPi * std::pow(alpha, 3), 1e-15);
  EXPECT_NEAR(atZero.secondDerivativeFactor, 0.8 * twoOverSqrtPi * std::pow(alpha, 5), 1e-15);

  // r = 0 with a scaled (offset) pair: finite, no direction
  const Potentials::EwaldExclusionFactors scaled = Potentials::ewaldExclusionFactors(alpha, 0.0, 0.0);
  EXPECT_TRUE(finite(scaled.potential) && finite(scaled.dUdlambda));
  EXPECT_EQ(scaled.firstDerivativeFactor, 0.0);
  EXPECT_EQ(scaled.secondDerivativeFactor, 0.0);

  // erf(x)/x = (2/sqrt(pi)) sum_k c_k x^2k, c_k = (-1)^k / (k! (2k + 1)), to twelve terms
  const auto seriesFactors = [&](double r)
  {
    const double x = alpha * r;
    double potential = 0.0, first = 0.0, curvature = 0.0, term = 1.0;
    for (int k = 0; k < 12; ++k)
    {
      const double c = term / (2.0 * k + 1.0);
      potential += c * std::pow(x, 2 * k);
      if (k >= 1) first += 2.0 * k * c * std::pow(x, 2 * k - 2);
      if (k >= 2) curvature += 4.0 * k * (k - 1) * c * std::pow(x, 2 * k - 4);
      term *= -1.0 / (k + 1.0);
    }
    return std::array<double, 3>{alpha * twoOverSqrtPi * potential, std::pow(alpha, 3) * twoOverSqrtPi * first,
                                 std::pow(alpha, 5) * twoOverSqrtPi * curvature};
  };
  for (const double r : {1.0e-6, 1.0e-3, 2.0e-2, 3.32e-2, 3.34e-2, 0.1})
  {
    const Potentials::EwaldExclusionFactors f = Potentials::ewaldExclusionFactors(alpha, 1.0, r);
    const std::array<double, 3> reference = seriesFactors(r);
    EXPECT_NEAR(f.potential, reference[0], 1e-14 * std::abs(reference[0])) << r;
    // above the switch (alpha r = 0.01) the closed forms keep about 1e-16 / (alpha r)^2 and 1e-16 / (alpha r)^4
    EXPECT_NEAR(f.firstDerivativeFactor, reference[1], 1e-11 * std::abs(reference[1])) << r;
    EXPECT_NEAR(f.secondDerivativeFactor, reference[2], 2e-7 * std::abs(reference[2])) << r;
  }
}
