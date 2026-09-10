#include <gtest/gtest.h>

import std;

import double4;
import units;
import pseudo_atom;
import vdwparameters;
import forcefield;
import potential_pair_derivatives;
import potential_pair_vdw;

// Tests for Class-II / PCFF (CFF) van der Waals cross interactions and mixing rules:
// - Sixth-power (Waldman-Hagler) mixing rule evaluation for CFFEpsilonSigma
// - String-to-enum alias resolution (cff_eps_sigma, class2, pcff, cff96)
// - Analytical energy and derivative agreement for mixed pairs

namespace
{
double analyticU96(double r, double eps, double sigma)
{
  double t = sigma / r;
  double t3 = t * t * t;
  double t6 = t3 * t3;
  double t9 = t6 * t3;
  return eps * (2.0 * t9 - 3.0 * t6);
}

double analyticGradientFactor96(double r, double eps, double sigma)
{
  double t = sigma / r;
  double t3 = t * t * t;
  double t6 = t3 * t3;
  double t9 = t6 * t3;
  return (18.0 * eps / (r * r)) * (t6 - t9);
}
}  // namespace

TEST(cff_cross_mixing, StringToEnumAliases)
{
  EXPECT_EQ(VDWParameters::stringToEnum("cff-eps-sigma"), VDWParameters::Type::CFFEpsilonSigma);
  EXPECT_EQ(VDWParameters::stringToEnum("cffepssigma"), VDWParameters::Type::CFFEpsilonSigma);
  EXPECT_EQ(VDWParameters::stringToEnum("cff_eps_sigma"), VDWParameters::Type::CFFEpsilonSigma);
  EXPECT_EQ(VDWParameters::stringToEnum("class2"), VDWParameters::Type::CFFEpsilonSigma);
  EXPECT_EQ(VDWParameters::stringToEnum("class-2"), VDWParameters::Type::CFFEpsilonSigma);
  EXPECT_EQ(VDWParameters::stringToEnum("pcff"), VDWParameters::Type::CFFEpsilonSigma);

  EXPECT_EQ(VDWParameters::stringToEnum("cff-9-6"), VDWParameters::Type::CFF9_6);
  EXPECT_EQ(VDWParameters::stringToEnum("cff9-6"), VDWParameters::Type::CFF9_6);
  EXPECT_EQ(VDWParameters::stringToEnum("cff_9_6"), VDWParameters::Type::CFF9_6);
  EXPECT_EQ(VDWParameters::stringToEnum("cff96"), VDWParameters::Type::CFF9_6);
}

TEST(cff_cross_mixing, SixthPowerMixingRule_MatchesAnalyticFormula)
{
  PseudoAtom atomA("C_cff", false, 12.011, 0.0, 0.0, 6, false);
  PseudoAtom atomB("H_cff", false, 1.008, 0.0, 0.0, 1, false);

  double epsA = 120.0;  // Kelvin
  double sigA = 3.60;   // Angstrom
  double epsB = 45.0;   // Kelvin
  double sigB = 2.80;   // Angstrom

  VDWParameters vdwA(VDWParameters::Type::CFFEpsilonSigma, {epsA, sigA});
  VDWParameters vdwB(VDWParameters::Type::CFFEpsilonSigma, {epsB, sigB});

  ForceField ff({atomA, atomB}, {vdwA, vdwB}, ForceField::MixingRule::SixthPower, 12.0, 12.0, 12.0, false, false,
                false);

  // Expected Waldman-Hagler (sixth-power) mixed parameters
  double s6A = std::pow(sigA, 6.0);
  double s6B = std::pow(sigB, 6.0);
  double expectedSigma = std::pow(0.5 * (s6A + s6B), 1.0 / 6.0);
  double expectedEps = 2.0 * std::sqrt(epsA * epsB) * (std::pow(sigA, 3.0) * std::pow(sigB, 3.0)) / (s6A + s6B);

  VDWParameters cross = ff(0, 1);
  EXPECT_EQ(cross.type, VDWParameters::Type::CFFEpsilonSigma);
  EXPECT_NEAR(cross.parameters.x * Units::EnergyToKelvin, expectedEps, 1e-8);
  EXPECT_NEAR(cross.parameters.y, expectedSigma, 1e-8);

  // Test pair energy evaluation at various distances against analytic 9-6 kernel
  for (double r : {2.5, 3.0, expectedSigma, 4.0, 5.5, 8.0})
  {
    double rr = r * r;
    Potentials::PairDerivatives<0> val = Potentials::potentialVDW<0>(ff, 1.0, 1.0, rr, 0, 1);
    double energyK = val.energy * Units::EnergyToKelvin;
    double expectedK = analyticU96(r, expectedEps, expectedSigma);
    EXPECT_NEAR(energyK, expectedK, std::max(1e-7, 1e-7 * std::abs(expectedK))) << " at r = " << r;
  }

  // At r = sigma, potential minimum is exactly -eps
  double rrMin = expectedSigma * expectedSigma;
  Potentials::PairDerivatives<0> valMin = Potentials::potentialVDW<0>(ff, 1.0, 1.0, rrMin, 0, 1);
  EXPECT_NEAR(valMin.energy * Units::EnergyToKelvin, -expectedEps, 1e-8);
}

TEST(cff_cross_mixing, SpatialDerivativesMatchAnalytic)
{
  PseudoAtom atomA("C_cff", false, 12.011, 0.0, 0.0, 6, false);
  PseudoAtom atomB("H_cff", false, 1.008, 0.0, 0.0, 1, false);

  double epsA = 120.0;
  double sigA = 3.60;
  double epsB = 45.0;
  double sigB = 2.80;

  VDWParameters vdwA(VDWParameters::Type::CFFEpsilonSigma, {epsA, sigA});
  VDWParameters vdwB(VDWParameters::Type::CFFEpsilonSigma, {epsB, sigB});

  ForceField ff({atomA, atomB}, {vdwA, vdwB}, ForceField::MixingRule::SixthPower, 12.0, 12.0, 12.0, false, false,
                false);

  double s6A = std::pow(sigA, 6.0);
  double s6B = std::pow(sigB, 6.0);
  double expectedSigma = std::pow(0.5 * (s6A + s6B), 1.0 / 6.0);
  double expectedEps = 2.0 * std::sqrt(epsA * epsB) * (std::pow(sigA, 3.0) * std::pow(sigB, 3.0)) / (s6A + s6B);

  for (double r : {3.0, expectedSigma, 4.0, 5.5, 7.0})
  {
    double rr = r * r;
    Potentials::PairDerivatives<1> val = Potentials::potentialVDW<1>(ff, 1.0, 1.0, rr, 0, 1);
    double gradFactor = val.firstDerivativeFactor * Units::EnergyToKelvin;
    double expectedGradFactor = analyticGradientFactor96(r, expectedEps, expectedSigma);
    EXPECT_NEAR(gradFactor, expectedGradFactor, std::max(1e-6, 1e-6 * std::abs(expectedGradFactor)))
        << " at r = " << r;
  }
}
