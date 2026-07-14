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
