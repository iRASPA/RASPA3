#include <gtest/gtest.h>

import std;

import double4;
import units;
import pseudo_atom;
import vdwparameters;
import forcefield;
import energy_factor;
import potential_energy_vdw;

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
      Potentials::EnergyFactor value =
          Potentials::potentialVDWEnergy(forceField, 1.0, 1.0, r * r, 0, 0);
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
    Potentials::EnergyFactor value =
        Potentials::potentialVDWEnergy(forceField, 1.0, 1.0, rc * rc, 0, 0);
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
        Potentials::EnergyFactor value =
            Potentials::potentialVDWEnergy(forceField, scaling, 1.0, rr, 0, 0);
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
        Potentials::EnergyFactor value =
            Potentials::potentialVDWEnergy(forceField, scalingA, scalingB, rr, 0, 0);

        double delta = 1e-6;
        Potentials::EnergyFactor forward =
            Potentials::potentialVDWEnergy(forceField, scalingA + 0.5 * delta, scalingB, rr, 0, 0);
        Potentials::EnergyFactor backward =
            Potentials::potentialVDWEnergy(forceField, scalingA - 0.5 * delta, scalingB, rr, 0, 0);
        double finiteDifference = (forward.energy - backward.energy) / delta;

        // dUdlambda holds the symmetric derivative factor X with dU/d(scalingA) = scalingB * X
        double tolerance = std::max(1e-4, 1e-6 * std::abs(finiteDifference));
        EXPECT_NEAR(scalingB * value.dUdlambda, finiteDifference, tolerance)
            << VDWParameters::nameOfType(testCase.type) << " at r = " << r << ", lambdaA = " << scalingA;
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
    Potentials::EnergyFactor below =
        Potentials::potentialVDWEnergy(forceField, scaling, 1.0, rrCut * (1.0 - 1e-12), 0, 0);
    EXPECT_NEAR(below.energy * Units::EnergyToKelvin, 0.0, 1e-6) << "lambda = " << scaling;
  }
}
