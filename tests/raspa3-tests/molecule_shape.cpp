#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import forcefield;
import component;
import molecule_property_settings;
import mc_moves_probabilities;
import property_molecule_shape;

// Tests for the gyration-tensor shape descriptors: analytic conformations (rod, square, tetrahedron)
// pin down Rg, the eigenvalues, asphericity, acylindricity, shape anisotropy and prolateness, and
// the Kirkwood hydrodynamic radius; the sampling path checks the moments, the ratio statistics and
// the histogram normalization / error bars.

namespace
{

// A linear four-bead chain.
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
      ["CH2", [3.71, 1.41, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3]
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

// A three-unit comb polymer with 'RepeatUnits' [backbone, side1, side2] per unit.
constexpr std::string_view kCombPolymerJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 1.54, 0.0]],
      ["CH2", [0.0, 2.17, 1.41]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [1.54, 1.54, 0.0]],
      ["CH2", [1.54, 2.17, 1.41]],
      ["CH2", [3.08, 0.0, 0.0]],
      ["CH2", [3.08, 1.54, 0.0]],
      ["CH2", [3.08, 2.17, 1.41]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [0, 3],
    [3, 4],
    [4, 5],
    [3, 6],
    [6, 7],
    [7, 8]
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
  "VanDerWaals" : "auto",
  "RepeatUnits" : [
    [0, 1, 2],
    [3, 4, 5],
    [6, 7, 8]
  ]
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

std::vector<double> uniformWeights(std::size_t n) { return std::vector<double>(n, 1.0 / static_cast<double>(n)); }

}  // namespace

// Four beads on a line at spacing 1: Rg^2 = variance of {0,1,2,3} = 1.25, l2 = l3 = 0, k^2 = 1, S = 2.
TEST(MOLECULE_SHAPE, rod_is_maximally_anisotropic_and_prolate)
{
  std::vector<Atom> rod = atomsAt({double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(2.0, 0.0, 0.0),
                                   double3(3.0, 0.0, 0.0)});
  auto d = PropertyMoleculeShape::computeDescriptors(rod, uniformWeights(4));

  EXPECT_NEAR(d.radiusOfGyrationSquared, 1.25, 1e-12);
  EXPECT_NEAR(d.eigenvalues.x, 1.25, 1e-12);
  EXPECT_NEAR(d.eigenvalues.y, 0.0, 1e-12);
  EXPECT_NEAR(d.eigenvalues.z, 0.0, 1e-12);
  EXPECT_NEAR(d.asphericity, 1.25, 1e-12);
  EXPECT_NEAR(d.acylindricity, 0.0, 1e-12);
  EXPECT_NEAR(d.shapeAnisotropy, 1.0, 1e-12);
  EXPECT_NEAR(d.prolateness, 2.0, 1e-12);
}

// The descriptors are rotation invariant: the same rod along a skew direction.
TEST(MOLECULE_SHAPE, rotation_invariance)
{
  double3 direction = double3(1.0, 2.0, -0.5).normalized();
  std::vector<double3> positions;
  for (int i = 0; i < 4; ++i) positions.push_back(static_cast<double>(i) * direction + double3(7.0, -3.0, 2.0));
  auto d = PropertyMoleculeShape::computeDescriptors(atomsAt(positions), uniformWeights(4));

  EXPECT_NEAR(d.radiusOfGyrationSquared, 1.25, 1e-10);
  EXPECT_NEAR(d.eigenvalues.x, 1.25, 1e-10);
  EXPECT_NEAR(d.eigenvalues.y, 0.0, 1e-10);
  EXPECT_NEAR(d.eigenvalues.z, 0.0, 1e-10);
  EXPECT_NEAR(d.shapeAnisotropy, 1.0, 1e-10);
  EXPECT_NEAR(d.prolateness, 2.0, 1e-10);

  // The lab-frame tensor of a rod is Rg^2 n n^T.
  EXPECT_NEAR(d.tensor.ax, 1.25 * direction.x * direction.x, 1e-10);
  EXPECT_NEAR(d.tensor.by, 1.25 * direction.y * direction.y, 1e-10);
  EXPECT_NEAR(d.tensor.cz, 1.25 * direction.z * direction.z, 1e-10);
  EXPECT_NEAR(d.tensor.ay, 1.25 * direction.x * direction.y, 1e-10);
  EXPECT_NEAR(d.tensor.az, 1.25 * direction.x * direction.z, 1e-10);
  EXPECT_NEAR(d.tensor.bz, 1.25 * direction.y * direction.z, 1e-10);
  EXPECT_NEAR(d.tensor.ay, d.tensor.bx, 1e-14);
}

// Four beads at the corners of a unit square in a tilted plane: l1 = l2 = 0.25, l3 = 0 (an oblate
// disk): b = 0.125, c = 0.25, k^2 = 1/4, S = -1/4.
TEST(MOLECULE_SHAPE, square_is_an_oblate_disk)
{
  // Orthonormal in-plane axes so the eigenvalues are unchanged by the tilt.
  double3 u = double3(1.0, 1.0, 0.0).normalized();
  double3 v = double3(-1.0, 1.0, 1.0).normalized();
  std::vector<Atom> square = atomsAt({-0.5 * u - 0.5 * v, 0.5 * u - 0.5 * v, 0.5 * u + 0.5 * v, -0.5 * u + 0.5 * v});
  auto d = PropertyMoleculeShape::computeDescriptors(square, uniformWeights(4));

  EXPECT_NEAR(d.radiusOfGyrationSquared, 0.5, 1e-12);
  EXPECT_NEAR(d.eigenvalues.x, 0.25, 1e-12);
  EXPECT_NEAR(d.eigenvalues.y, 0.25, 1e-12);
  EXPECT_NEAR(d.eigenvalues.z, 0.0, 1e-12);
  EXPECT_NEAR(d.asphericity, 0.125, 1e-12);
  EXPECT_NEAR(d.acylindricity, 0.25, 1e-12);
  EXPECT_NEAR(d.shapeAnisotropy, 0.25, 1e-12);
  EXPECT_NEAR(d.prolateness, -0.25, 1e-12);
}

// A regular tetrahedron is isotropic: all eigenvalues equal, b = c = k^2 = S = 0.
TEST(MOLECULE_SHAPE, tetrahedron_is_spherical)
{
  std::vector<Atom> tetrahedron = atomsAt(
      {double3(1.0, 1.0, 1.0), double3(1.0, -1.0, -1.0), double3(-1.0, 1.0, -1.0), double3(-1.0, -1.0, 1.0)});
  auto d = PropertyMoleculeShape::computeDescriptors(tetrahedron, uniformWeights(4));

  EXPECT_NEAR(d.radiusOfGyrationSquared, 3.0, 1e-12);
  EXPECT_NEAR(d.eigenvalues.x, 1.0, 1e-12);
  EXPECT_NEAR(d.eigenvalues.y, 1.0, 1e-12);
  EXPECT_NEAR(d.eigenvalues.z, 1.0, 1e-12);
  EXPECT_NEAR(d.asphericity, 0.0, 1e-12);
  EXPECT_NEAR(d.acylindricity, 0.0, 1e-12);
  EXPECT_NEAR(d.shapeAnisotropy, 0.0, 1e-12);
  EXPECT_NEAR(d.prolateness, 0.0, 1e-12);
}

// Mass weights move the center and reweight the tensor: two beads with masses 1 and 3 at distance 4
// have the center at 3 from the light bead and Rg^2 = (1/4) 9 + (3/4) 1 = 3.
TEST(MOLECULE_SHAPE, mass_weighting)
{
  std::vector<Atom> dumbbell = atomsAt({double3(0.0, 0.0, 0.0), double3(4.0, 0.0, 0.0)});
  auto uniform = PropertyMoleculeShape::computeDescriptors(dumbbell, uniformWeights(2));
  EXPECT_NEAR(uniform.radiusOfGyrationSquared, 4.0, 1e-12);

  std::vector<double> weights{0.25, 0.75};
  auto weighted = PropertyMoleculeShape::computeDescriptors(dumbbell, weights);
  EXPECT_NEAR(weighted.radiusOfGyrationSquared, 3.0, 1e-12);
  EXPECT_NEAR(weighted.shapeAnisotropy, 1.0, 1e-12);
}

// Kirkwood sum for the four-bead rod: pairs at distance 1 (x3), 2 (x2), 3 (x1); ordered double sum
// is twice that, divided by N^2 = 16.
TEST(MOLECULE_SHAPE, kirkwood_inverse_hydrodynamic_radius)
{
  std::vector<Atom> rod = atomsAt({double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(2.0, 0.0, 0.0),
                                   double3(3.0, 0.0, 0.0)});
  double expected = 2.0 * (3.0 / 1.0 + 2.0 / 2.0 + 1.0 / 3.0) / 16.0;
  EXPECT_NEAR(PropertyMoleculeShape::computeInverseHydrodynamicRadius(rod), expected, 1e-12);
}

// Sampling through a component: the moments, the derived ratio statistics, the histogram
// normalization and the block error bars.
TEST(MOLECULE_SHAPE, sampling_moments_ratios_and_histograms)
{
  ForceField forceField = makeZeroForceField();
  std::vector<Component> components{};
  components.push_back(makeComponent(forceField, "linear-shape", kLinearChainJson));

  constexpr std::size_t numberOfBlocks = 5;
  // The Rg histogram is sized by its bin width: ceil(5.0 / 0.1) = 50 bins; the k^2 and S histograms
  // keep the fixed 'numberOfBins'.
  constexpr std::size_t numberOfBins = 50;
  constexpr double rgRange = 5.0;
  constexpr double deltaRg = 0.1;
  components[0].moleculeShapeSettings = MoleculeShapeSettings{.sampleEvery = 1,
                                                              .writeEvery = 1,
                                                              .numberOfBins = 40,
                                                              .massWeighted = false,
                                                              .radiusOfGyrationRange = rgRange,
                                                              .radiusOfGyrationBinWidth = deltaRg};
  PropertyMoleculeShape property(numberOfBlocks, forceField, components);

  ASSERT_TRUE(property.isSampled(0));
  ASSERT_EQ(property.weightsPerComponent[0].size(), 4);
  ASSERT_NEAR(property.deltaRadiusOfGyrationPerComponent[0], deltaRg, 1e-12);
  ASSERT_EQ(property.numberOfRadiusOfGyrationBinsPerComponent[0], numberOfBins);
  ASSERT_EQ(property.radiusOfGyrationHistogram[0][0].size(), numberOfBins);
  ASSERT_EQ(property.shapeAnisotropyHistogram[0][0].size(), 40uz);
  ASSERT_TRUE(property.endToEndAtomsPerComponent[0].has_value());

  std::vector<Atom> rod = components[0].atoms;
  for (std::size_t i = 0; i < 4; ++i) rod[i].position = double3(static_cast<double>(i), 0.0, 0.0);
  // Rod: Rg^2 = 1.25, end-to-end (0, 3) squared = 9.
  std::vector<Atom> tetrahedron = components[0].atoms;
  tetrahedron[0].position = double3(1.0, 1.0, 1.0);
  tetrahedron[1].position = double3(1.0, -1.0, -1.0);
  tetrahedron[2].position = double3(-1.0, 1.0, -1.0);
  tetrahedron[3].position = double3(-1.0, -1.0, 1.0);
  // Tetrahedron: Rg^2 = 3, end-to-end (0, 3) squared = |(2, 2, 0)|^2 = 8.

  for (std::size_t block = 0; block != numberOfBlocks; ++block)
  {
    property.sample(components, {1}, std::span<const Atom>(rod), 0, block);
    property.sample(components, {1}, std::span<const Atom>(tetrahedron), 0, block);
  }

  auto [meanRg2, errorRg2] = property.momentStatistics(0, PropertyMoleculeShape::RadiusOfGyrationSquared);
  EXPECT_NEAR(meanRg2, 0.5 * (1.25 + 3.0), 1e-12);
  EXPECT_EQ(errorRg2, 0.0);

  // Total tensor: trace equals <Rg^2>; the rod (along x) contributes 1.25 to S_xx only, the
  // tetrahedron is isotropic with 1 on each diagonal entry and zero off-diagonal.
  auto [sxx, eSxx] = property.momentStatistics(0, PropertyMoleculeShape::TensorXX);
  auto [syy, eSyy] = property.momentStatistics(0, PropertyMoleculeShape::TensorYY);
  auto [szz, eSzz] = property.momentStatistics(0, PropertyMoleculeShape::TensorZZ);
  auto [sxy, eSxy] = property.momentStatistics(0, PropertyMoleculeShape::TensorXY);
  EXPECT_NEAR(sxx, 0.5 * (1.25 + 1.0), 1e-12);
  EXPECT_NEAR(syy, 0.5, 1e-12);
  EXPECT_NEAR(szz, 0.5, 1e-12);
  EXPECT_NEAR(sxy, 0.0, 1e-12);
  EXPECT_NEAR(sxx + syy + szz, meanRg2, 1e-12);

  auto [meanKappa2, errorKappa2] = property.momentStatistics(0, PropertyMoleculeShape::ShapeAnisotropy);
  EXPECT_NEAR(meanKappa2, 0.5 * (1.0 + 0.0), 1e-12);

  auto [meanS, errorS] = property.momentStatistics(0, PropertyMoleculeShape::Prolateness);
  EXPECT_NEAR(meanS, 0.5 * (2.0 + 0.0), 1e-12);

  auto [meanR2, errorR2] = property.momentStatistics(0, PropertyMoleculeShape::EndToEndSquared);
  EXPECT_NEAR(meanR2, 0.5 * (9.0 + 8.0), 1e-12);

  // Ratio of block-combined means.
  auto [ratio, errorRatio] = property.statistics(
      0, [](const PropertyMoleculeShape::Moments &m)
      { return m[PropertyMoleculeShape::EndToEndSquared] / m[PropertyMoleculeShape::RadiusOfGyrationSquared]; });
  EXPECT_NEAR(ratio, 8.5 / 2.125, 1e-12);
  EXPECT_EQ(errorRatio, 0.0);

  // Ensemble-form anisotropy: rod has I2 = 0, Rg^4 = 1.5625; tetrahedron I2 = 3, Rg^4 = 9.
  auto [ensembleKappa2, errorEnsembleKappa2] = property.statistics(
      0, [](const PropertyMoleculeShape::Moments &m)
      {
        return 1.0 - 3.0 * m[PropertyMoleculeShape::SecondInvariant] / m[PropertyMoleculeShape::RadiusOfGyrationFourth];
      });
  EXPECT_NEAR(ensembleKappa2, 1.0 - 3.0 * 1.5 / 5.28125, 1e-12);

  // Rg histogram: normalized, half the mass at each of the two Rg values.
  auto [values, average, error] =
      property.result(property.radiusOfGyrationHistogram, 0, property.deltaRadiusOfGyrationPerComponent[0], 0.0);
  double integral = 0.0;
  for (std::size_t bin = 0; bin != numberOfBins; ++bin) integral += average[bin] * deltaRg;
  EXPECT_NEAR(integral, 1.0, 1e-12);
  std::size_t binRod = static_cast<std::size_t>(std::sqrt(1.25) / deltaRg);
  std::size_t binTet = static_cast<std::size_t>(std::sqrt(3.0) / deltaRg);
  EXPECT_NEAR(average[binRod], 0.5 / deltaRg, 1e-12);
  EXPECT_NEAR(average[binTet], 0.5 / deltaRg, 1e-12);
  for (std::size_t bin = 0; bin != numberOfBins; ++bin) EXPECT_EQ(error[bin], 0.0);

  // Anisotropy histogram: k^2 = 1 sits on the upper edge of the [0, 1] range and is discarded from the
  // histogram (the moments keep it); k^2 = 0 lands in the first bin.
  const double deltaShapeAnisotropy = property.deltaShapeAnisotropyPerComponent[0];
  const double deltaProlateness = property.deltaProlatenessPerComponent[0];
  auto [kValues, kAverage, kError] = property.result(property.shapeAnisotropyHistogram, 0, deltaShapeAnisotropy, 0.0);
  EXPECT_NEAR(kAverage[0], 0.5 / deltaShapeAnisotropy, 1e-12);

  // Prolateness histogram: S = 0 and S = 2 (the upper edge, discarded).
  auto [sValues, sAverage, sError] =
      property.result(property.prolatenessHistogram, 0, deltaProlateness, property.prolatenessLowerLimit);
  std::size_t binZero = static_cast<std::size_t>((0.0 - property.prolatenessLowerLimit) / deltaProlateness);
  EXPECT_NEAR(sAverage[binZero], 0.5 / deltaProlateness, 1e-12);

  // Skew one block; the errors become positive.
  property.sample(components, {1}, std::span<const Atom>(rod), 0, 0);
  auto [meanRg2Skewed, errorRg2Skewed] = property.momentStatistics(0, PropertyMoleculeShape::RadiusOfGyrationSquared);
  EXPECT_GT(errorRg2Skewed, 0.0);
  auto [valuesSkewed, averageSkewed, errorSkewed] =
      property.result(property.radiusOfGyrationHistogram, 0, property.deltaRadiusOfGyrationPerComponent[0], 0.0);
  EXPECT_GT(errorSkewed[binRod], 0.0);
}

// Per-monomer descriptors for a chain with repeat units: each unit's own gyration tensor, reported
// per unit index and pooled over units. Units 0 and 2 are placed collinear (k^2 = 1, S = 2), unit 1
// as an equilateral triangle (l1 = l2 = a^2/6, l3 = 0: k^2 = 1/4, S = -1/4).
TEST(MOLECULE_SHAPE, per_monomer_descriptors)
{
  ForceField forceField = makeZeroForceField();
  std::vector<Component> components{};
  components.push_back(makeComponent(forceField, "comb-monomers", kCombPolymerJson));

  constexpr std::size_t numberOfBlocks = 3;
  components[0].moleculeShapeSettings = MoleculeShapeSettings{
      .sampleEvery = 1, .writeEvery = 1, .numberOfBins = 16, .massWeighted = false, .radiusOfGyrationRange = 10.0};
  PropertyMoleculeShape property(numberOfBlocks, forceField, components);
  ASSERT_EQ(property.numberOfUnits(0), 3);
  ASSERT_EQ(property.unitAtomsPerComponent[0][1], (std::vector<std::size_t>{3, 4, 5}));
  for (double w : property.unitWeightsPerComponent[0][0]) EXPECT_NEAR(w, 1.0 / 3.0, 1e-12);

  constexpr double d = 1.5;  // collinear spacing
  constexpr double a = 2.0;  // triangle side
  std::vector<Atom> molecule = components[0].atoms;
  molecule[0].position = double3(0.0, 0.0, 0.0);
  molecule[1].position = double3(d, 0.0, 0.0);
  molecule[2].position = double3(2.0 * d, 0.0, 0.0);
  molecule[3].position = double3(10.0, 0.0, 0.0);
  molecule[4].position = double3(10.0 + a, 0.0, 0.0);
  molecule[5].position = double3(10.0 + 0.5 * a, 0.5 * std::sqrt(3.0) * a, 0.0);
  molecule[6].position = double3(20.0, 0.0, 0.0);
  molecule[7].position = double3(20.0, d, 0.0);
  molecule[8].position = double3(20.0, 2.0 * d, 0.0);

  for (std::size_t block = 0; block != numberOfBlocks; ++block)
  {
    property.sample(components, {1}, std::span<const Atom>(molecule), 0, block);
  }

  auto rg2 = [](const PropertyMoleculeShape::Moments &m) { return m[PropertyMoleculeShape::RadiusOfGyrationSquared]; };
  auto kappa2 = [](const PropertyMoleculeShape::Moments &m) { return m[PropertyMoleculeShape::ShapeAnisotropy]; };
  auto prolateness = [](const PropertyMoleculeShape::Moments &m) { return m[PropertyMoleculeShape::Prolateness]; };
  auto lambda3 = [](const PropertyMoleculeShape::Moments &m) { return m[PropertyMoleculeShape::Lambda3]; };

  // Collinear three points at spacing d: Rg^2 = 2 d^2 / 3.
  EXPECT_NEAR(property.unitStatistics(0, 0, rg2).first, 2.0 * d * d / 3.0, 1e-12);
  EXPECT_NEAR(property.unitStatistics(0, 0, kappa2).first, 1.0, 1e-12);
  EXPECT_NEAR(property.unitStatistics(0, 0, prolateness).first, 2.0, 1e-12);
  EXPECT_NEAR(property.unitStatistics(0, 2, kappa2).first, 1.0, 1e-12);

  // Equilateral triangle: Rg^2 = a^2 / 3, oblate disk.
  EXPECT_NEAR(property.unitStatistics(0, 1, rg2).first, a * a / 3.0, 1e-12);
  EXPECT_NEAR(property.unitStatistics(0, 1, lambda3).first, 0.0, 1e-12);
  EXPECT_NEAR(property.unitStatistics(0, 1, kappa2).first, 0.25, 1e-12);
  EXPECT_NEAR(property.unitStatistics(0, 1, prolateness).first, -0.25, 1e-12);

  // Pooled over the three units: every unit is one sample.
  auto [pooledKappa2, pooledError] = property.unitStatistics(0, std::nullopt, kappa2);
  EXPECT_NEAR(pooledKappa2, (1.0 + 0.25 + 1.0) / 3.0, 1e-12);
  EXPECT_EQ(pooledError, 0.0);
  EXPECT_NEAR(property.unitStatistics(0, std::nullopt, rg2).first, (2.0 * d * d / 3.0 * 2.0 + a * a / 3.0) / 3.0,
              1e-12);

  // The whole-molecule descriptors are unaffected by the per-unit bookkeeping.
  auto whole = PropertyMoleculeShape::computeDescriptors(molecule, property.weightsPerComponent[0]);
  EXPECT_NEAR(property.momentStatistics(0, PropertyMoleculeShape::RadiusOfGyrationSquared).first,
              whole.radiusOfGyrationSquared, 1e-12);
}

// The default Rg range is derived from the contour length of the bond-graph diameter (3 x 1.54 for
// the linear chain); components without settings are not sampled.
TEST(MOLECULE_SHAPE, default_range_from_contour_length)
{
  ForceField forceField = makeZeroForceField();
  std::vector<Component> components{};
  components.push_back(makeComponent(forceField, "linear-range", kLinearChainJson));
  components.push_back(makeComponent(forceField, "linear-unsampled", kLinearChainJson));

  components[0].moleculeShapeSettings =
      MoleculeShapeSettings{.sampleEvery = 1, .writeEvery = 1, .numberOfBins = 64, .massWeighted = true};
  PropertyMoleculeShape property(3, forceField, components);
  ASSERT_TRUE(property.isSampled(0));
  EXPECT_FALSE(property.isSampled(1));
  EXPECT_TRUE(property.weightsPerComponent[1].empty());
  // The harmonic equilibrium length is located on a 5 Angstrom / 1024-point grid (0.005 resolution);
  // the range is then rounded up to a whole number of 0.1 Angstrom bins (default width).
  const double requested = 0.6 * 3.0 * 1.54;
  EXPECT_EQ(property.numberOfRadiusOfGyrationBinsPerComponent[0], static_cast<std::size_t>(std::ceil(requested / 0.1)));
  EXPECT_NEAR(property.radiusOfGyrationRangePerComponent[0],
              0.1 * static_cast<double>(property.numberOfRadiusOfGyrationBinsPerComponent[0]), 1e-12);
  EXPECT_GE(property.radiusOfGyrationRangePerComponent[0], requested - 0.6 * 3.0 * 0.005);
  EXPECT_LT(property.radiusOfGyrationRangePerComponent[0], requested + 0.1 + 0.6 * 3.0 * 0.005);

  // Mass weights of identical beads are uniform.
  for (double w : property.weightsPerComponent[0]) EXPECT_NEAR(w, 0.25, 1e-12);
}
