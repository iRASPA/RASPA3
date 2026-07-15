#include <gtest/gtest.h>

import std;

import double4;
import units;
import pseudo_atom;
import vdwparameters;
import forcefield;
import potential_pair_derivatives;
import potential_pair_vdw;

// Tests for the van der Waals potentials ported from RASPA2:
// - value parity with the RASPA2 functional forms at full coupling (lambda = 1)
// - the shifted potential is zero at the cutoff distance
// - finiteness at r -> 0 for fractional lambda (soft-core scaling)
// - the analytic dU/dlambda matches a finite-difference derivative

namespace
{

struct PotentialTestCase
{
  VDWParameters::Type type;
  std::vector<double> values;  // raw force-field input (energies in Kelvin)
  std::vector<double> distances;
};

std::vector<PotentialTestCase> testCases()
{
  return {
      {VDWParameters::Type::LennardJones, {119.8, 3.405}, {3.0, 3.405, 3.8, 5.0}},
      {VDWParameters::Type::BuckingHam, {3.0e5, 3.5, 1.2e4}, {2.0, 2.45, 3.0, 5.0}},
      {VDWParameters::Type::Morse, {500.0, 1.5, 2.0}, {1.0, 2.0, 3.0, 5.0}},
      {VDWParameters::Type::FeynmannHibbs, {36.7, 2.958, 1.0}, {2.6, 2.958, 3.4, 5.0}},
      {VDWParameters::Type::MM3, {120.0, 3.5}, {1.0, 3.2, 3.9, 5.0}},
      {VDWParameters::Type::BornHugginsMeyer, {6.0e4, 3.15, 2.34, 1.0e6, 2.0e6}, {3.0, 4.0, 4.5, 6.0}},
      {VDWParameters::Type::LennardJonesShiftedForce, {119.8, 3.405}, {3.0, 3.405, 3.8, 5.0}},
      {VDWParameters::Type::LennardJonesSecondOrderTaylorShifted, {119.8, 3.405}, {3.0, 3.405, 3.8, 5.0}},
      {VDWParameters::Type::Potential12_6, {6.0e7, 2.5e4}, {3.2, 3.66, 4.0, 5.0}},
      {VDWParameters::Type::Potential12_6_2_0, {6.0e7, -2.5e4, -100.0, 0.0}, {3.2, 3.66, 4.0, 5.0}},
      {VDWParameters::Type::CFF9_6, {1.0e6, 2.0e4}, {3.2, 3.68, 4.0, 5.0}},
      {VDWParameters::Type::CFFEpsilonSigma, {120.0, 3.4}, {3.0, 3.4, 3.8, 5.0}},
      {VDWParameters::Type::MatsuokaClementiYoshimine, {2.0e5, 3.0, -1.0e3, 1.2}, {1.0, 2.0, 3.0, 5.0}},
      {VDWParameters::Type::Generic, {3.0e5, 3.2, 1.0e3, 1.0e4, 5.0e4, 1.0e5}, {2.0, 2.2, 3.0, 5.0}},
      {VDWParameters::Type::PellenqNicholson, {1.0e5, 3.28, 1.0e4, 5.0e4, 2.0e5}, {1.0, 1.2, 2.0, 5.0}},
      {VDWParameters::Type::HydratedIonWater, {1.0e7, 3.5, 1.0e3, 1.0e4, 1.0e6}, {3.0, 4.3, 5.0, 6.0}},
      {VDWParameters::Type::Mie, {6.0e7, 12.0, 2.5e4, 6.0}, {3.2, 3.66, 4.0, 5.0}},
      {VDWParameters::Type::WeeksChandlerAndersen, {119.8, 3.405}, {3.0, 3.405, 3.7, 3.8}},
      {VDWParameters::Type::RepulsiveHarmonic, {500.0, 2.5}, {1.0, 2.0, 3.0, 5.0}},
  };
}

ForceField makeForceField(const PotentialTestCase &testCase, bool shifted)
{
  PseudoAtom pseudoAtom("A", false, 12.0, 0.0, 0.0, 6, false);
  VDWParameters vdw(testCase.type, testCase.values);
  return ForceField({pseudoAtom}, {vdw}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, shifted, false,
                    false);
}

// Independent implementation of the Tang-Toennies damping factor 1 - exp(-b r) sum_{k=0}^{n} (b r)^k / k!
double dampingFactor(std::size_t n, double b, double r)
{
  double sum = 0.0;
  double factorial = 1.0;
  for (std::size_t k = 0; k <= n; ++k)
  {
    if (k > 0) factorial *= static_cast<double>(k);
    sum += std::pow(b * r, static_cast<double>(k)) / factorial;
  }
  return 1.0 - std::exp(-b * r) * sum;
}

// Reference energies in Kelvin, written directly from the RASPA2 functional forms in potentials.c.
double referenceEnergyInKelvin(const PotentialTestCase &testCase, double r)
{
  const std::vector<double> &p = testCase.values;
  switch (testCase.type)
  {
    case VDWParameters::Type::LennardJones:
    {
      double u6 = std::pow(p[1] / r, 6);
      return 4.0 * p[0] * (u6 * u6 - u6);
    }
    case VDWParameters::Type::BuckingHam:
      return p[0] * std::exp(-p[1] * r) - p[2] / std::pow(r, 6);
    case VDWParameters::Type::Morse:
    {
      double expTerm = std::exp(-p[1] * (r - p[2]));
      return p[0] * ((1.0 - expTerm) * (1.0 - expTerm) - 1.0);
    }
    case VDWParameters::Type::FeynmannHibbs:
    {
      double u6 = std::pow(p[1] / r, 6);
      double fh = Units::FeymannHibbsConversionFactor / (p[2] * 300.0);
      return 4.0 * p[0] * (u6 * u6 - u6) + fh * 4.0 * p[0] * (132.0 * u6 * u6 - 30.0 * u6) / (r * r);
    }
    case VDWParameters::Type::MM3:
    {
      double P = p[1] / r;
      if (P > 3.02) return p[0] * 192.270 * P * P;
      return p[0] * (1.84e5 * std::exp(-12.0 / P) - 2.25 * std::pow(P, 6));
    }
    case VDWParameters::Type::BornHugginsMeyer:
      return p[0] * std::exp(p[1] * (p[2] - r)) - p[3] / std::pow(r, 6) - p[4] / std::pow(r, 8);
    case VDWParameters::Type::LennardJonesShiftedForce:
    {
      double rc = 12.0;
      double u6 = std::pow(p[1] / r, 6);
      double c6 = std::pow(p[1] / rc, 6);
      return 4.0 * p[0] * (u6 * u6 - u6 - (c6 * c6 - c6) + (12.0 * c6 * c6 - 6.0 * c6) * (r - rc) / rc);
    }
    case VDWParameters::Type::LennardJonesSecondOrderTaylorShifted:
    {
      double rc = 12.0;
      double u6 = std::pow(p[1] / r, 6);
      double c6 = std::pow(p[1] / rc, 6);
      double displacement = (r - rc) / rc;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      double quadraticCoefficient = 156.0 * c6 * c6 - 42.0 * c6;
      return 4.0 * p[0] *
             (u6 * u6 - u6 - (c6 * c6 - c6) + linearCoefficient * displacement -
              0.5 * quadraticCoefficient * displacement * displacement);
    }
    case VDWParameters::Type::Potential12_6:
      return p[0] / std::pow(r, 12) - p[1] / std::pow(r, 6);
    case VDWParameters::Type::Potential12_6_2_0:
      return p[0] / std::pow(r, 12) + p[1] / std::pow(r, 6) + p[2] / (r * r) + p[3];
    case VDWParameters::Type::CFF9_6:
      return p[0] / std::pow(r, 9) - p[1] / std::pow(r, 6);
    case VDWParameters::Type::CFFEpsilonSigma:
      return p[0] * (2.0 * std::pow(p[1] / r, 9) - 3.0 * std::pow(p[1] / r, 6));
    case VDWParameters::Type::MatsuokaClementiYoshimine:
      return p[0] * std::exp(-p[1] * r) + p[2] * std::exp(-p[3] * r);
    case VDWParameters::Type::Generic:
      return p[0] * std::exp(-p[1] * r) - p[2] / std::pow(r, 4) - p[3] / std::pow(r, 6) - p[4] / std::pow(r, 8) -
             p[5] / std::pow(r, 10);
    case VDWParameters::Type::PellenqNicholson:
      return p[0] * std::exp(-p[1] * r) - dampingFactor(6, p[1], r) * p[2] / std::pow(r, 6) -
             dampingFactor(8, p[1], r) * p[3] / std::pow(r, 8) - dampingFactor(10, p[1], r) * p[4] / std::pow(r, 10);
    case VDWParameters::Type::HydratedIonWater:
      return p[0] * std::exp(-p[1] * r) - p[2] / std::pow(r, 4) - p[3] / std::pow(r, 6) - p[4] / std::pow(r, 12);
    case VDWParameters::Type::Mie:
      return p[0] / std::pow(r, p[1]) - p[2] / std::pow(r, p[3]);
    case VDWParameters::Type::WeeksChandlerAndersen:
    {
      if (r > std::pow(2.0, 1.0 / 6.0) * p[1]) return 0.0;
      double u6 = std::pow(p[1] / r, 6);
      return 4.0 * p[0] * (u6 * u6 - u6) + p[0];
    }
    case VDWParameters::Type::RepulsiveHarmonic:
    {
      if (r >= p[1]) return 0.0;
      double displacement = 1.0 - r / p[1];
      return 0.5 * p[0] * displacement * displacement;
    }
    default:
      return 0.0;
  }
}

}  // namespace

TEST(vdw_potentials, RASPA2_reference_values_at_full_coupling)
{
  for (const PotentialTestCase &testCase : testCases())
  {
    ForceField forceField = makeForceField(testCase, false);
    for (double r : testCase.distances)
    {
      Potentials::PairDerivatives<0> value =
          Potentials::potentialVDW<0>(forceField, 1.0, 1.0, r * r, 0, 0);
      double energyInKelvin = value.energy * Units::EnergyToKelvin;
      double reference = referenceEnergyInKelvin(testCase, r);
      double tolerance = std::max(1e-8, 1e-8 * std::abs(reference));
      EXPECT_NEAR(energyInKelvin, reference, tolerance)
          << VDWParameters::nameOfType(testCase.type) << " at r = " << r;
    }
  }
}

TEST(vdw_potentials, shifted_potential_is_zero_at_cutoff)
{
  for (const PotentialTestCase &testCase : testCases())
  {
    ForceField forceField = makeForceField(testCase, true);
    double rc = 12.0;
    Potentials::PairDerivatives<0> value =
        Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rc * rc, 0, 0);
    EXPECT_NEAR(value.energy * Units::EnergyToKelvin, 0.0, 1e-8) << VDWParameters::nameOfType(testCase.type);
  }
}

TEST(vdw_potentials, finite_at_zero_distance_for_fractional_lambda)
{
  for (const PotentialTestCase &testCase : testCases())
  {
    ForceField forceField = makeForceField(testCase, true);
    for (double scaling : {0.0, 0.2, 0.5, 0.8, 0.95})
    {
      for (double rr : {1e-10, 1e-4, 0.01, 1.0})
      {
        Potentials::PairDerivatives<0> value =
            Potentials::potentialVDW<0>(forceField, scaling, 1.0, rr, 0, 0);
        EXPECT_TRUE(std::isfinite(value.energy))
            << VDWParameters::nameOfType(testCase.type) << " at rr = " << rr << ", lambda = " << scaling;
        EXPECT_TRUE(std::isfinite(value.dUdlambda))
            << VDWParameters::nameOfType(testCase.type) << " at rr = " << rr << ", lambda = " << scaling;
      }
    }
  }
}

TEST(vdw_potentials, lambda_derivative_matches_finite_difference)
{
  for (const PotentialTestCase &testCase : testCases())
  {
    ForceField forceField = makeForceField(testCase, true);
    double scalingB = 0.9;
    for (double scalingA : {0.3, 0.7})
    {
      std::vector<double> distances = testCase.distances;
      distances.push_back(0.5);  // deep inside the repulsive core, where the soft-core matters
      for (double r : distances)
      {
        double rr = r * r;
        Potentials::PairDerivatives<0> value =
            Potentials::potentialVDW<0>(forceField, scalingA, scalingB, rr, 0, 0);

        double delta = 1e-6;
        Potentials::PairDerivatives<0> forward =
            Potentials::potentialVDW<0>(forceField, scalingA + 0.5 * delta, scalingB, rr, 0, 0);
        Potentials::PairDerivatives<0> backward =
            Potentials::potentialVDW<0>(forceField, scalingA - 0.5 * delta, scalingB, rr, 0, 0);
        double finiteDifference = (forward.energy - backward.energy) / delta;

        // dUdlambda holds the symmetric derivative factor X with dU/d(scalingA) = scalingB * X
        double tolerance = 1e-5 * std::max(1.0, std::abs(finiteDifference));
        EXPECT_NEAR(scalingB * value.dUdlambda, finiteDifference, tolerance)
            << VDWParameters::nameOfType(testCase.type) << " at r = " << r << ", lambdaA = " << scalingA;
      }
    }
  }
}

TEST(vdw_potentials, all_spatial_derivatives_match_finite_difference)
{
  for (const PotentialTestCase &testCase : testCases())
  {
    ForceField forceField = makeForceField(testCase, true);
    for (const auto [scalingA, scalingB] : {std::pair{1.0, 1.0}, std::pair{0.7, 0.9}})
    {
      for (double r : testCase.distances)
      {
        const double rr = r * r;
        const Potentials::PairDerivatives<0> energy =
            Potentials::potentialVDW<0>(forceField, scalingA, scalingB, rr, 0, 0);
        const Potentials::PairDerivatives<1> gradient =
            Potentials::potentialVDW<1>(forceField, scalingA, scalingB, rr, 0, 0);
        const Potentials::PairDerivatives<2> hessian =
            Potentials::potentialVDW<2>(forceField, scalingA, scalingB, rr, 0, 0);

        const double fieldScale = std::max({1.0, std::abs(energy.energy), std::abs(energy.dUdlambda)});
        EXPECT_NEAR(gradient.energy, energy.energy, 1e-12 * fieldScale);
        EXPECT_NEAR(hessian.energy, energy.energy, 1e-12 * fieldScale);
        EXPECT_NEAR(gradient.dUdlambda, energy.dUdlambda, 1e-12 * fieldScale);
        EXPECT_NEAR(hessian.dUdlambda, energy.dUdlambda, 1e-12 * fieldScale);
        EXPECT_NEAR(gradient.firstDerivativeFactor, hessian.firstDerivativeFactor,
                    1e-12 * std::max({1.0, std::abs(gradient.firstDerivativeFactor),
                                      std::abs(hessian.firstDerivativeFactor)}));

        const double step = 1e-4 * std::max(1.0, r);
        auto energyAt = [&](double distance)
        {
          return Potentials::potentialVDW<0>(forceField, scalingA, scalingB, distance * distance, 0, 0).energy;
        };
        const double minus2 = energyAt(r - 2.0 * step);
        const double minus1 = energyAt(r - step);
        const double center = energyAt(r);
        const double plus1 = energyAt(r + step);
        const double plus2 = energyAt(r + 2.0 * step);
        const double numericalFirst = (minus2 - 8.0 * minus1 + 8.0 * plus1 - plus2) / (12.0 * step);
        const double numericalSecond =
            (-plus2 + 16.0 * plus1 - 30.0 * center + 16.0 * minus1 - minus2) /
            (12.0 * step * step);
        const double analyticFirst = gradient.firstDerivativeFactor * r;
        const double analyticSecond =
            hessian.firstDerivativeFactor + hessian.secondDerivativeFactor * rr;

        EXPECT_NEAR(analyticFirst, numericalFirst,
                    1e-5 * std::max({1.0, std::abs(analyticFirst), std::abs(numericalFirst)}))
            << VDWParameters::nameOfType(testCase.type) << " at r = " << r
            << ", lambda = " << scalingA * scalingB;
        EXPECT_NEAR(analyticSecond, numericalSecond,
                    1e-5 * std::max({1.0, std::abs(analyticSecond), std::abs(numericalSecond)}))
            << VDWParameters::nameOfType(testCase.type) << " at r = " << r
            << ", lambda = " << scalingA * scalingB;
      }
    }
  }
}

TEST(vdw_potentials, wca_continuous_at_internal_cutoff)
{
  PotentialTestCase testCase{VDWParameters::Type::WeeksChandlerAndersen, {119.8, 3.405}, {}};
  ForceField forceField = makeForceField(testCase, false);
  double sigma = 3.405;
  double rrCut = std::cbrt(2.0) * sigma * sigma;
  for (double scaling : {0.2, 0.5, 0.8, 1.0})
  {
    Potentials::PairDerivatives<0> below =
        Potentials::potentialVDW<0>(forceField, scaling, 1.0, rrCut * (1.0 - 1e-12), 0, 0);
    EXPECT_NEAR(below.energy * Units::EnergyToKelvin, 0.0, 1e-6) << "lambda = " << scaling;

    Potentials::PairDerivatives<1> outsideGradient =
        Potentials::potentialVDW<1>(forceField, scaling, 1.0, rrCut * (1.0 + 1e-8), 0, 0);
    Potentials::PairDerivatives<2> outsideHessian =
        Potentials::potentialVDW<2>(forceField, scaling, 1.0, rrCut * (1.0 + 1e-8), 0, 0);
    EXPECT_EQ(outsideGradient.energy, 0.0);
    EXPECT_EQ(outsideGradient.firstDerivativeFactor, 0.0);
    EXPECT_EQ(outsideHessian.firstDerivativeFactor, 0.0);
    EXPECT_EQ(outsideHessian.secondDerivativeFactor, 0.0);
  }
}

TEST(vdw_potentials, second_order_taylor_shifted_spatial_derivatives_match_finite_difference)
{
  PotentialTestCase testCase{VDWParameters::Type::LennardJonesSecondOrderTaylorShifted, {119.8, 3.405}, {}};
  ForceField forceField = makeForceField(testCase, true);

  for (double r : {3.0, 3.405, 3.8, 5.0, 11.0})
  {
    constexpr double delta = 1e-5;
    const Potentials::PairDerivatives<1> gradient =
        Potentials::potentialVDW<1>(forceField, 1.0, 1.0, r * r, 0, 0);
    const Potentials::PairDerivatives<2> hessian =
        Potentials::potentialVDW<2>(forceField, 1.0, 1.0, r * r, 0, 0);

    const double plus = Potentials::potentialVDW<0>(
                            forceField, 1.0, 1.0, (r + delta) * (r + delta), 0, 0)
                            .energy;
    const double center = Potentials::potentialVDW<0>(forceField, 1.0, 1.0, r * r, 0, 0).energy;
    const double minus = Potentials::potentialVDW<0>(
                             forceField, 1.0, 1.0, (r - delta) * (r - delta), 0, 0)
                             .energy;
    const double numericalFirst = (plus - minus) / (2.0 * delta);
    const double numericalSecond = (plus - 2.0 * center + minus) / (delta * delta);

    EXPECT_NEAR(gradient.firstDerivativeFactor * r, numericalFirst,
                1e-5 * std::max(1.0, std::abs(numericalFirst)))
        << "r = " << r;
    EXPECT_NEAR(hessian.firstDerivativeFactor + hessian.secondDerivativeFactor * r * r,
                numericalSecond, 1e-5 * std::max(1.0, std::abs(numericalSecond)))
        << "r = " << r;
  }
}

TEST(vdw_potentials, second_order_taylor_shifted_energy_gradient_and_hessian_vanish_at_cutoff)
{
  PotentialTestCase testCase{VDWParameters::Type::LennardJonesSecondOrderTaylorShifted, {119.8, 3.405}, {}};
  ForceField forceField = makeForceField(testCase, true);
  constexpr double cutoff = 12.0;

  const Potentials::PairDerivatives<1> value =
      Potentials::potentialVDW<1>(forceField, 1.0, 1.0, cutoff * cutoff, 0, 0);
  const Potentials::PairDerivatives<2> hessian =
      Potentials::potentialVDW<2>(forceField, 1.0, 1.0, cutoff * cutoff, 0, 0);
  EXPECT_NEAR(value.energy, 0.0, 1e-14);
  EXPECT_NEAR(value.firstDerivativeFactor, 0.0, 1e-14);
  EXPECT_NEAR(hessian.secondDerivativeFactor, 0.0, 1e-14);
}
