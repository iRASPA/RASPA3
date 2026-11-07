#include <gtest/gtest.h>

import std;

import double3;
import units;
import randomnumbers;
import urey_bradley_potential;
import bond_bond_potential;
import bend_torsion_potential;

TEST(INTRAMOLECULAR_POTENTIAL_FIXES, urey_bradley_lengths_follow_the_boltzmann_distribution)
{
  // The 1-3 distance is distributed as r^2 exp(-beta u(r)); the sampled mean must match the exact mean of that
  // density on the support of the sampler.
  const double beta = 1.0 / (Units::KB * 300.0);
  struct Case
  {
    UreyBradleyPotential potential;
    double lo, hi;
  };
  const std::array<Case, 3> cases{
      Case{UreyBradleyPotential({0, 1}, UreyBradleyType::Harmonic, {2000.0, 2.0}), 0.0, 6.0},
      Case{UreyBradleyPotential({0, 1}, UreyBradleyType::Morse, {2000.0, 2.0, 1.8}), 0.0, 3.0},
      Case{UreyBradleyPotential({0, 1}, UreyBradleyType::Quartic, {1000.0, 1.5, -200.0, 500.0}), 0.0, 3.0}};

  for (const auto& [potential, lo, hi] : cases)
  {
    constexpr std::size_t numberOfPoints = 200000;
    double weight = 0.0;
    double weightedLength = 0.0;
    for (std::size_t i = 0; i < numberOfPoints; ++i)
    {
      const double r = lo + (hi - lo) * (static_cast<double>(i) + 0.5) / static_cast<double>(numberOfPoints);
      const double density =
          r * r * std::exp(-beta * potential.calculateEnergy(double3(0.0, 0.0, 0.0), double3(r, 0.0, 0.0)));
      weight += density;
      weightedLength += r * density;
    }
    const double expectedMean = weightedLength / weight;

    RandomNumber random(17);
    constexpr std::size_t numberOfSamples = 200000;
    double sampledMean = 0.0;
    for (std::size_t i = 0; i < numberOfSamples; ++i) sampledMean += potential.generateUreyBradleyLength(random, beta);
    sampledMean /= static_cast<double>(numberOfSamples);

    EXPECT_NEAR(sampledMean, expectedMean, 4e-3) << "type " << static_cast<std::size_t>(potential.type);
  }
}

TEST(INTRAMOLECULAR_POTENTIAL_FIXES, bond_bond_with_fewer_parameters_reads_only_the_given_ones)
{
  const BondBondPotential potential({0, 1, 2}, BondBondType::CVFF, std::vector<double>{500.0});
  EXPECT_DOUBLE_EQ(potential.parameters[0], 500.0 * Units::KelvinToEnergy);
  EXPECT_EQ(potential.parameters[1], 0.0);
  EXPECT_EQ(potential.parameters[2], 0.0);
}

TEST(INTRAMOLECULAR_POTENTIAL_FIXES, smoothed_three_cosine_bend_torsion_gradient_matches_finite_difference)
{
  // Both bend angles near 175 degrees, inside the smoothing range (above 170 degrees) of S(theta).
  const BendTorsionPotential potential({0, 1, 2, 3}, BendTorsionType::SmoothedThreeCosine, {300.0, 150.0, 400.0});
  const double angle = 5.0 * std::numbers::pi / 180.0;
  std::array<double3, 4> positions{
      double3(-1.5 * std::cos(angle), 1.5 * std::sin(angle), 0.0), double3(0.0, 0.0, 0.0), double3(1.5, 0.0, 0.0),
      double3(1.5 + 1.5 * std::cos(angle), 0.6 * 1.5 * std::sin(angle), 0.8 * 1.5 * std::sin(angle))};

  const auto [energy, gradient, strain] =
      potential.potentialEnergyGradientStrain(positions[0], positions[1], positions[2], positions[3]);
  ASSERT_NE(energy, 0.0);

  constexpr double delta = 1.0e-6;
  for (std::size_t atom = 0; atom < 4; ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      std::array<double3, 4> plus = positions;
      std::array<double3, 4> minus = positions;
      plus[atom][axis] += delta;
      minus[atom][axis] -= delta;
      const double numerical = (potential.calculateEnergy(plus[0], plus[1], plus[2], plus[3]) -
                                potential.calculateEnergy(minus[0], minus[1], minus[2], minus[3])) /
                               (2.0 * delta);
      EXPECT_NEAR(gradient[atom][axis], numerical, 1e-5 * std::max(1.0, std::abs(numerical)))
          << "atom " << atom << ", axis " << axis;
    }
  }
}
