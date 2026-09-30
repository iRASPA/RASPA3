#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import forcefield;
import component;
import mc_moves_probabilities;
import property_molecule_backbone;

// Tests for the backbone chain statistics: internal distances, bond-vector correlation, form
// factor and the derived chain descriptors on analytic conformations, plus the sampling and
// block-statistics path through a component.

namespace
{

// A linear six-bead chain (backbone = graph diameter = all beads).
constexpr std::string_view kLinearChainJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.17, 1.41, 0.0]],
      ["CH2", [3.71, 1.41, 0.0]],
      ["CH2", [4.34, 2.82, 0.0]],
      ["CH2", [5.88, 2.82, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [700.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
  ],
  "VanDerWaals" : "auto"
}
)";

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}}, {{0.0, 3.95}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

Component makeComponent(const ForceField &forceField, const std::string &name, std::string_view json)
{
  TemporaryFile file(name + ".json", json);
  return Component(Component::Type::Adsorbate, 0, forceField, name, file.stemPath().string(), 5, 21,
                   MCMoveProbabilities(), std::nullopt, false);
}

std::vector<Atom> atomsAt(const std::vector<double3> &positions)
{
  std::vector<Atom> atoms(positions.size());
  for (std::size_t i = 0; i < positions.size(); ++i) atoms[i].position = positions[i];
  return atoms;
}

std::vector<std::size_t> identityBackbone(std::size_t n)
{
  std::vector<std::size_t> backbone(n);
  std::iota(backbone.begin(), backbone.end(), 0uz);
  return backbone;
}

// Rod of n beads with spacing l along a skew direction.
std::vector<Atom> makeRod(std::size_t n, double l)
{
  double3 direction = double3(1.0, -2.0, 0.5).normalized();
  std::vector<double3> positions;
  for (std::size_t i = 0; i < n; ++i) positions.push_back(static_cast<double>(i) * l * direction);
  return atomsAt(positions);
}

}  // namespace

// Rod: every bond parallel, C(k) = 1, <r^2(k)> = (k l)^2, projection equals the full length.
TEST(MOLECULE_BACKBONE, rod_internal_distances_correlation_and_moments)
{
  constexpr std::size_t n = 7;
  constexpr double l = 1.3;
  std::vector<Atom> rod = makeRod(n, l);
  std::vector<std::size_t> backbone = identityBackbone(n);

  std::vector<double> r2(n - 1), corr(n - 1);
  PropertyMoleculeBackbone::accumulateInternalDistances(rod, backbone, r2);
  PropertyMoleculeBackbone::accumulateBondCorrelation(rod, backbone, corr);
  for (std::size_t k = 1; k < n; ++k)
  {
    EXPECT_NEAR(r2[k - 1], static_cast<double>(k * k) * l * l, 1e-10);
    EXPECT_NEAR(corr[k - 1], 1.0, 1e-12);
  }

  auto m = PropertyMoleculeBackbone::computeMoments(rod, backbone);
  EXPECT_NEAR(m[PropertyMoleculeBackbone::BondLength], l, 1e-12);
  EXPECT_NEAR(m[PropertyMoleculeBackbone::BondLengthSquared], l * l, 1e-12);
  EXPECT_NEAR(m[PropertyMoleculeBackbone::EndToEndSquared], 36.0 * l * l, 1e-10);
  EXPECT_NEAR(m[PropertyMoleculeBackbone::Projection], 6.0 * l, 1e-10);
  // Rg^2 of n equally spaced points: l^2 (n^2 - 1) / 12.
  EXPECT_NEAR(m[PropertyMoleculeBackbone::RadiusOfGyrationSquared], l * l * (n * n - 1.0) / 12.0, 1e-10);
}

// Planar zig-zag with 90-degree bends: successive bonds orthogonal, C(k) alternates 1, 0, -1, 0, ...
TEST(MOLECULE_BACKBONE, zigzag_bond_correlation)
{
  std::vector<Atom> zigzag = atomsAt({double3(0, 0, 0), double3(1, 0, 0), double3(1, 1, 0), double3(0, 1, 0),
                                      double3(0, 2, 0), double3(1, 2, 0)});
  // Bonds: +x, +y, -x, +y, +x  ->  C(1) = 0, C(2) = mean(-1, 1, -1) = -1/3, C(3) = mean(0, 0) = 0, C(4) = 1.
  std::vector<double> corr(5);
  PropertyMoleculeBackbone::accumulateBondCorrelation(zigzag, identityBackbone(6), corr);
  EXPECT_NEAR(corr[0], 1.0, 1e-12);
  EXPECT_NEAR(corr[1], 0.0, 1e-12);
  EXPECT_NEAR(corr[2], -1.0 / 3.0, 1e-12);
  EXPECT_NEAR(corr[3], 0.0, 1e-12);
  EXPECT_NEAR(corr[4], 1.0, 1e-12);

  // Projection onto the first bond (+x): R = (1, 2, 0) -> 1; onto the last bond (+x) -> 1.
  auto m = PropertyMoleculeBackbone::computeMoments(zigzag, identityBackbone(6));
  EXPECT_NEAR(m[PropertyMoleculeBackbone::Projection], 1.0, 1e-12);
  EXPECT_NEAR(m[PropertyMoleculeBackbone::EndToEndSquared], 5.0, 1e-12);
}

// Form factor: P(q -> 0) = 1 exactly, and the two-bead value is (1 + sinc(q d)) / 2.
TEST(MOLECULE_BACKBONE, form_factor_limits)
{
  std::vector<Atom> dumbbell = atomsAt({double3(0, 0, 0), double3(2.0, 0, 0)});
  std::vector<double> q{1e-9, 0.5, 2.0};
  std::vector<double> p(3);
  PropertyMoleculeBackbone::accumulateFormFactor(dumbbell, q, p);
  EXPECT_NEAR(p[0], 1.0, 1e-9);
  EXPECT_NEAR(p[1], 0.5 * (1.0 + std::sin(1.0) / 1.0), 1e-12);
  EXPECT_NEAR(p[2], 0.5 * (1.0 + std::sin(4.0) / 4.0), 1e-12);

  // Guinier regime for the rod: P(q) ~ 1 - q^2 Rg^2 / 3.
  constexpr std::size_t n = 20;
  std::vector<Atom> rod = makeRod(n, 1.0);
  double rg2 = (n * n - 1.0) / 12.0;
  std::vector<double> qSmall{0.01};
  std::vector<double> pSmall(1);
  PropertyMoleculeBackbone::accumulateFormFactor(rod, qSmall, pSmall);
  // The next term of the expansion is O(q^4 <r^4>) ~ 5e-7 here.
  EXPECT_NEAR(pSmall[0], 1.0 - 0.01 * 0.01 * rg2 / 3.0, 2e-6);
}

// The histogram route reproduces the direct evaluation: pair distances binned at 0.005 Angstrom and
// transformed with sin(q r_b)/(q r_b) at the bin centres agree with the per-pair sum to the
// discretization error ~ (q dr)^2 / 24, and the histogram grows on demand.
TEST(MOLECULE_BACKBONE, form_factor_from_histogram_matches_direct)
{
  // Two irregular 16-bead conformations (deterministic pseudo-random walk).
  std::vector<std::vector<Atom>> molecules{};
  std::uint64_t state = 12345;
  auto next = [&state]()
  {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11) / 9007199254740992.0;
  };
  for (std::size_t m = 0; m < 2; ++m)
  {
    std::vector<double3> positions{double3(0.0, 0.0, 0.0)};
    for (std::size_t i = 1; i < 16; ++i)
    {
      double3 step(next() - 0.5, next() - 0.5, next() - 0.5);
      positions.push_back(positions.back() + 1.54 * step.normalized());
    }
    molecules.push_back(atomsAt(positions));
  }

  std::vector<double> q{0.01, 0.3, 1.0, 2.5, 5.0};
  constexpr double binWidth = 0.005;

  std::vector<double> direct(q.size());
  std::vector<double> histogram(10);  // deliberately too short: must grow
  for (const std::vector<Atom> &molecule : molecules)
  {
    PropertyMoleculeBackbone::accumulateFormFactor(molecule, q, direct);
    PropertyMoleculeBackbone::accumulatePairDistances(molecule, binWidth, histogram);
  }
  for (double &v : direct) v /= 2.0;

  EXPECT_GT(histogram.size(), 10);
  double pairs = std::accumulate(histogram.begin(), histogram.end(), 0.0);
  EXPECT_EQ(pairs, 2.0 * 16 * 15 / 2);

  std::vector<double> fromHistogram = PropertyMoleculeBackbone::formFactorFromHistogram(histogram, binWidth, 16, 2.0, q);
  ASSERT_EQ(fromHistogram.size(), q.size());
  // Rigorous bound: a pair displaced by at most dr/2 within its bin changes sin(x)/x by at most
  // q dr/4 (|d/dx sin(x)/x| <= 1/2), and the pair sum carries weight (N-1)/N < 1.
  for (std::size_t iq = 0; iq < q.size(); ++iq)
  {
    double tolerance = 1e-9 + 0.25 * q[iq] * binWidth;
    EXPECT_NEAR(fromHistogram[iq], direct[iq], tolerance) << "q = " << q[iq];
  }
  EXPECT_NEAR(fromHistogram[0], 1.0, 1e-3);

  // No molecules: zero, not NaN.
  std::vector<double> empty = PropertyMoleculeBackbone::formFactorFromHistogram(histogram, binWidth, 16, 0.0, q);
  for (double v : empty) EXPECT_EQ(v, 0.0);
}

TEST(MOLECULE_BACKBONE, debye_function)
{
  EXPECT_NEAR(PropertyMoleculeBackbone::debyeFunction(0.0), 1.0, 1e-12);
  EXPECT_NEAR(PropertyMoleculeBackbone::debyeFunction(1e-8), 1.0, 1e-8);
  // x = 2: 2 (e^-2 - 1 + 2) / 4 = (1 + e^-2) / 2.
  EXPECT_NEAR(PropertyMoleculeBackbone::debyeFunction(2.0), 0.5 * (1.0 + std::exp(-2.0)), 1e-12);
  // Large x: -> 2 / x.
  EXPECT_NEAR(PropertyMoleculeBackbone::debyeFunction(1000.0), 2.0 * 999.0 / 1e6, 1e-12);
}

// Derived descriptors on synthetic averages: an exactly exponential C(k) recovers l_p, an exact
// power law recovers nu, and the undetermined cases return zero.
TEST(MOLECULE_BACKBONE, fits_recover_exponential_and_power_law)
{
  PropertyMoleculeBackbone::Averages a{};
  constexpr std::size_t numberOfBeads = 64;
  constexpr double l = 1.5;
  constexpr double lp = 6.0;
  constexpr double nu = 0.588;
  a.moments[PropertyMoleculeBackbone::BondLength] = l;
  a.bondCorrelation.resize(numberOfBeads - 1);
  for (std::size_t k = 0; k < a.bondCorrelation.size(); ++k)
  {
    a.bondCorrelation[k] = std::exp(-static_cast<double>(k) * l / lp);
  }
  a.internalDistanceSquared.resize(numberOfBeads - 1);
  for (std::size_t k = 1; k < numberOfBeads; ++k)
  {
    a.internalDistanceSquared[k - 1] = 2.3 * std::pow(static_cast<double>(k), 2.0 * nu);
  }
  a.moments[PropertyMoleculeBackbone::Projection] = 4.2;

  EXPECT_NEAR(PropertyMoleculeBackbone::persistenceLengthFromFit(a), lp, 1e-9);
  EXPECT_NEAR(PropertyMoleculeBackbone::floryExponent(a), nu, 1e-9);
  EXPECT_NEAR(PropertyMoleculeBackbone::persistenceLengthFromProjection(a), 4.2, 1e-12);

  // Too few points: a five-bead chain has no fit window; a correlation at the noise floor after one
  // step leaves fewer than three fit points.
  PropertyMoleculeBackbone::Averages b{};
  b.moments[PropertyMoleculeBackbone::BondLength] = l;
  b.bondCorrelation = {1.0, 0.02, 0.01, 0.001};
  b.internalDistanceSquared = {1.0, 2.0, 3.0, 4.0};
  EXPECT_EQ(PropertyMoleculeBackbone::persistenceLengthFromFit(b), 0.0);
  EXPECT_EQ(PropertyMoleculeBackbone::floryExponent(b), 0.0);
}

// Sampling through a component: backbone inference, normalization of the per-k averages over
// molecules and blocks, ratio statistics and error bars, and the wave-vector grid.
TEST(MOLECULE_BACKBONE, sampling_through_component)
{
  ForceField forceField = makeZeroForceField();
  std::vector<Component> components{};
  components.push_back(makeComponent(forceField, "linear-backbone", kLinearChainJson));

  constexpr std::size_t numberOfBlocks = 4;
  PropertyMoleculeBackbone property(numberOfBlocks, components, 5, 0.1, 10.0, 1, 1);

  ASSERT_TRUE(property.isSampled(0));
  EXPECT_EQ(property.numberOfBackboneBeads(0), 6);
  EXPECT_EQ(property.backbonePerComponent[0].front(), 0);
  EXPECT_EQ(property.backbonePerComponent[0].back(), 5);
  // Contour length: five harmonic bonds located on the 5 Angstrom / 1024-point grid.
  EXPECT_NEAR(property.contourLengthPerComponent[0], 5.0 * 1.54, 5.0 * 0.005);

  // Log-spaced wave vectors from 0.1 to 10 in five points: 0.1, 0.316, 1, 3.16, 10.
  ASSERT_EQ(property.waveVectors.size(), 5);
  EXPECT_NEAR(property.waveVectors[0], 0.1, 1e-12);
  EXPECT_NEAR(property.waveVectors[2], 1.0, 1e-12);
  EXPECT_NEAR(property.waveVectors[4], 10.0, 1e-12);

  // Two molecules per sample: a rod with spacing 1 and a rod with spacing 2 (both along x).
  std::vector<Atom> pair(12);
  for (std::size_t i = 0; i < 6; ++i)
  {
    pair[i].position = double3(static_cast<double>(i), 0.0, 0.0);
    pair[6 + i].position = double3(2.0 * static_cast<double>(i), 0.0, 0.0);
  }
  for (std::size_t block = 0; block != numberOfBlocks; ++block)
  {
    property.sample(components, {2}, std::span<const Atom>(pair), 0, block);
  }

  // <r^2(k)> averaged over the two rods: k^2 (1 + 4) / 2.
  for (std::size_t k = 1; k < 6; ++k)
  {
    auto [value, error] = property.statistics(
        0, [k](const PropertyMoleculeBackbone::Averages &a) { return a.internalDistanceSquared[k - 1]; });
    EXPECT_NEAR(value, 2.5 * static_cast<double>(k * k), 1e-10);
    EXPECT_EQ(error, 0.0);
  }
  auto [corr3, errorCorr3] =
      property.statistics(0, [](const PropertyMoleculeBackbone::Averages &a) { return a.bondCorrelation[3]; });
  EXPECT_NEAR(corr3, 1.0, 1e-12);

  // <l> = 1.5, <R^2> = (25 + 100) / 2 = 62.5, C_N = <R^2> / (5 <l>^2) = 62.5 / 11.25.
  auto [meanL, errorL] = property.statistics(
      0, [](const PropertyMoleculeBackbone::Averages &a) { return a.moments[PropertyMoleculeBackbone::BondLength]; });
  EXPECT_NEAR(meanL, 1.5, 1e-12);
  auto [cn, errorCn] = property.statistics(
      0, [](const PropertyMoleculeBackbone::Averages &a)
      {
        double l = a.moments[PropertyMoleculeBackbone::BondLength];
        return a.moments[PropertyMoleculeBackbone::EndToEndSquared] / (5.0 * l * l);
      });
  EXPECT_NEAR(cn, 62.5 / 11.25, 1e-12);
  EXPECT_EQ(errorCn, 0.0);

  // Form factor through the sampling path (pair-distance histogram): the average of the two rods'
  // direct evaluations, identical in every block so the error is zero.
  {
    std::vector<double> direct(property.waveVectors.size());
    PropertyMoleculeBackbone::accumulateFormFactor(std::span<const Atom>(pair).subspan(0, 6), property.waveVectors,
                                                   direct);
    PropertyMoleculeBackbone::accumulateFormFactor(std::span<const Atom>(pair).subspan(6, 6), property.waveVectors,
                                                   direct);
    for (std::size_t iq = 0; iq < property.waveVectors.size(); ++iq)
    {
      auto [value, error] =
          property.statistics(0, [iq](const PropertyMoleculeBackbone::Averages &a) { return a.formFactor[iq]; });
      double q = property.waveVectors[iq];
      EXPECT_NEAR(value, 0.5 * direct[iq], 1e-9 + 0.25 * q * property.pairDistanceBinWidth);
      EXPECT_EQ(error, 0.0);
    }
  }

  // Skewing one block with an extra molecule makes the block scatter, and so the error, positive.
  property.sample(components, {1}, std::span<const Atom>(pair).subspan(0, 6), 0, 0);
  auto [meanLSkewed, errorLSkewed] = property.statistics(
      0, [](const PropertyMoleculeBackbone::Averages &a) { return a.moments[PropertyMoleculeBackbone::BondLength]; });
  EXPECT_GT(errorLSkewed, 0.0);
  EXPECT_LT(meanLSkewed, 1.5);
}
