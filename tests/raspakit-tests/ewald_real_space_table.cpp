#include <gtest/gtest.h>

import std;

import forcefield;
import potential_pair_derivatives;
import potential_pair_coulomb;
import potential_ewald_real_space_table;

namespace
{
ForceField makeEwaldForceField(double alpha, double cutOff)
{
  ForceField forceField;
  forceField.chargeMethod = ForceField::ChargeMethod::Ewald;
  forceField.cutOffCoulomb = cutOff;
  forceField.EwaldAlpha = alpha;
  forceField.updateEwaldRealSpaceTable();
  return forceField;
}

ForceField withoutTable(ForceField forceField)
{
  forceField.useEwaldRealSpaceTable = false;
  forceField.updateEwaldRealSpaceTable();
  return forceField;
}

void expectRelativeNear(double value, double reference, double tolerance, const char* what, double r)
{
  EXPECT_NEAR(value, reference, tolerance * std::max(1.0e-300, std::abs(reference))) << what << " at r = " << r;
}
}  // namespace

TEST(ewald_real_space_table, table_is_built_for_alpha_and_cutoff_and_cleared_when_disabled)
{
  const ForceField forceField = makeEwaldForceField(0.265, 12.0);
  const EwaldRealSpaceTable& table = forceField.ewaldRealSpaceTable;
  EXPECT_EQ(table.alpha, 0.265);
  EXPECT_TRUE(table.spans(12.0));
  EXPECT_TRUE(table.matches(0.265, 12.0));
  EXPECT_FALSE(table.matches(0.3, 12.0));
  EXPECT_FALSE(table.matches(0.265, 13.0));
  EXPECT_FALSE(table.rangeCapped);
  // 16384 intervals for the 12 Angstrom cutoff at the default spacing (plus the two guard nodes)
  EXPECT_NEAR(static_cast<double>(table.value.size()), 16384.0 + 2.0, 40.0);

  const ForceField exactOnly = withoutTable(forceField);
  EXPECT_TRUE(exactOnly.ewaldRealSpaceTable.value.empty());
  EXPECT_FALSE(exactOnly.ewaldRealSpaceTable.covers(0.265, 4.0));
}

TEST(ewald_real_space_table, large_cutoff_caps_the_range_and_keeps_the_resolution)
{
  const ForceField forceField = makeEwaldForceField(0.01, 500.0);
  const EwaldRealSpaceTable& table = forceField.ewaldRealSpaceTable;
  EXPECT_TRUE(table.rangeCapped);
  EXPECT_FALSE(table.spans(500.0));
  EXPECT_TRUE(table.matches(0.01, 500.0));  // at the cap: no rebuild for a larger cutoff either
  EXPECT_TRUE(table.matches(0.01, 600.0));
  EXPECT_EQ(table.spacing, EwaldRealSpaceTable::defaultSpacing);
  EXPECT_EQ(table.value.size(), EwaldRealSpaceTable::defaultMaximumIntervals + 2);

  // inside the covered range the table is used, beyond it the closed form: both must agree with the exact value
  const ForceField exactOnly = withoutTable(forceField);
  for (const double r : {1.5, 2.0, 10.0, 40.0, 60.0, 250.0, 499.0})
  {
    const Potentials::PairDerivatives<1> tabulated =
        Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, -1.0);
    const Potentials::PairDerivatives<1> exact = Potentials::potentialCoulomb<1>(exactOnly, 1.0, 1.0, r, 1.0, -1.0);
    expectRelativeNear(tabulated.energy, exact.energy, 1.0e-9, "energy", r);
    expectRelativeNear(tabulated.firstDerivativeFactor, exact.firstDerivativeFactor, 1.0e-8, "gradient factor", r);
  }
}

TEST(ewald_real_space_table, interpolated_energy_and_gradient_match_closed_form)
{
  // Error envelope of the cubic Hermite interpolant in r^2 at the default spacing: the fourth derivative of
  // erfc(alpha r)/r with respect to r^2 scales as r^-9, so the error is largest at the 1 Angstrom lower bound
  // and negligible at intermolecular distances.
  struct Band
  {
    double rMax;
    double energyTolerance;
    double gradientTolerance;
  };
  const std::array<Band, 3> bands{Band{1.5, 1.0e-8, 3.0e-7}, Band{2.5, 1.0e-9, 1.0e-8}, Band{12.0, 1.0e-10, 1.0e-9}};

  for (const double alpha : {0.15, 0.265, 0.4})
  {
    const ForceField forceField = makeEwaldForceField(alpha, 12.0);
    const ForceField exactOnly = withoutTable(forceField);
    ASSERT_FALSE(forceField.ewaldRealSpaceTable.value.empty());

    // dense sweep over the covered range, including points that straddle nodes
    for (double r = 1.0; r < 12.0; r += 0.0137)
    {
      const Band& band = *std::ranges::find_if(bands, [r](const Band& b) { return r < b.rMax; });

      const Potentials::PairDerivatives<0> energyTabulated =
          Potentials::potentialCoulomb<0>(forceField, 1.0, 1.0, r, 0.8, -0.4);
      const Potentials::PairDerivatives<0> energyExact =
          Potentials::potentialCoulomb<0>(exactOnly, 1.0, 1.0, r, 0.8, -0.4);
      expectRelativeNear(energyTabulated.energy, energyExact.energy, band.energyTolerance, "energy", r);
      expectRelativeNear(energyTabulated.dUdlambda, energyExact.dUdlambda, band.gradientTolerance, "dUdlambda", r);

      const Potentials::PairDerivatives<1> gradientTabulated =
          Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 0.8, -0.4);
      const Potentials::PairDerivatives<1> gradientExact =
          Potentials::potentialCoulomb<1>(exactOnly, 1.0, 1.0, r, 0.8, -0.4);
      expectRelativeNear(gradientTabulated.energy, gradientExact.energy, band.energyTolerance, "energy (order 1)", r);
      expectRelativeNear(gradientTabulated.firstDerivativeFactor, gradientExact.firstDerivativeFactor,
                         band.gradientTolerance, "gradient factor", r);
      // orders 0 and 1 read the same table entry
      EXPECT_EQ(gradientTabulated.energy, energyTabulated.energy);
    }
  }
}

TEST(ewald_real_space_table, closed_form_is_used_below_one_angstrom_for_scaled_pairs_and_for_stale_alpha)
{
  ForceField forceField = makeEwaldForceField(0.265, 12.0);
  const ForceField exactOnly = withoutTable(forceField);

  // Closed-form results from two inlined instances may differ by an ulp (contraction), hence EXPECT_DOUBLE_EQ;
  // the interpolant would differ by orders of magnitude more.

  // below the table's lower bound the closed form is used
  for (const double r : {0.05, 0.3, 0.7, 0.999})
  {
    const Potentials::PairDerivatives<1> a = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
    const Potentials::PairDerivatives<1> b = Potentials::potentialCoulomb<1>(exactOnly, 1.0, 1.0, r, 1.0, 1.0);
    EXPECT_DOUBLE_EQ(a.energy, b.energy);
    EXPECT_DOUBLE_EQ(a.firstDerivativeFactor, b.firstDerivativeFactor);
  }

  // fractional pairs keep the offset form of the scaled Ewald term
  for (const double scaling : {0.0, 0.25, 0.9})
  {
    const Potentials::PairDerivatives<1> a = Potentials::potentialCoulomb<1>(forceField, scaling, 1.0, 3.0, 1.0, 1.0);
    const Potentials::PairDerivatives<1> b = Potentials::potentialCoulomb<1>(exactOnly, scaling, 1.0, 3.0, 1.0, 1.0);
    EXPECT_DOUBLE_EQ(a.energy, b.energy);
    EXPECT_DOUBLE_EQ(a.dUdlambda, b.dUdlambda);
    EXPECT_DOUBLE_EQ(a.firstDerivativeFactor, b.firstDerivativeFactor);
  }

  // Hessian orders are closed-form throughout
  {
    const Potentials::PairDerivatives<2> a = Potentials::potentialCoulomb<2>(forceField, 1.0, 1.0, 3.0, 1.0, 1.0);
    const Potentials::PairDerivatives<2> b = Potentials::potentialCoulomb<2>(exactOnly, 1.0, 1.0, 3.0, 1.0, 1.0);
    EXPECT_DOUBLE_EQ(a.energy, b.energy);
    EXPECT_DOUBLE_EQ(a.firstDerivativeFactor, b.firstDerivativeFactor);
    EXPECT_DOUBLE_EQ(a.secondDerivativeFactor, b.secondDerivativeFactor);
  }

  // a direct change of alpha without a rebuild makes the table stale: the closed form at the new alpha is used
  forceField.EwaldAlpha = 0.31;
  EXPECT_FALSE(forceField.ewaldRealSpaceTable.covers(forceField.EwaldAlpha, 9.0));
  {
    const ForceField exactAtNewAlpha = withoutTable(forceField);
    const Potentials::PairDerivatives<1> a = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, 3.0, 1.0, 1.0);
    const Potentials::PairDerivatives<1> b = Potentials::potentialCoulomb<1>(exactAtNewAlpha, 1.0, 1.0, 3.0, 1.0, 1.0);
    EXPECT_DOUBLE_EQ(a.energy, b.energy);
    EXPECT_DOUBLE_EQ(a.firstDerivativeFactor, b.firstDerivativeFactor);
  }
  // and the rebuild picks the new alpha up
  forceField.updateEwaldRealSpaceTable();
  EXPECT_TRUE(forceField.ewaldRealSpaceTable.covers(0.31, 9.0));
}

TEST(ewald_real_space_table, tabulated_gradient_is_the_derivative_of_the_tabulated_energy)
{
  // The Hermite interpolant is C1 and the force is its exact derivative: a central difference of the tabulated
  // energy reproduces the tabulated gradient factor to finite-difference accuracy, also across node boundaries.
  const ForceField forceField = makeEwaldForceField(0.265, 12.0);
  const double step = 1.0e-6;
  for (double r = 1.2; r < 11.5; r += 0.31)
  {
    const double energyPlus = Potentials::potentialCoulomb<0>(forceField, 1.0, 1.0, r + step, 1.0, 1.0).energy;
    const double energyMinus = Potentials::potentialCoulomb<0>(forceField, 1.0, 1.0, r - step, 1.0, 1.0).energy;
    const double numerical = (energyPlus - energyMinus) / (2.0 * step) / r;
    const double analytic = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0).firstDerivativeFactor;
    EXPECT_NEAR(numerical, analytic, 1.0e-6 * std::abs(analytic) + 1.0e-6) << "r = " << r;
  }
}
