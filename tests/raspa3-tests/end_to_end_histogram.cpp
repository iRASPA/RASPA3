#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import forcefield;
import component;
import mc_moves_probabilities;
import property_molecule_properties;

// Tests for the end-to-end distance histogram: the inference of the chain ends (explicit
// 'EndToEndAtoms', 'RepeatUnits' backbone ends, or the bond-graph diameter) and the histogram /
// error-bar machinery of PropertyMoleculeProperties.

namespace
{

// A linear four-bead chain: the ends are the graph-diameter endpoints (0, 3).
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

// The same chain with an explicit override: 'EndToEndAtoms' always wins over inference.
constexpr std::string_view kExplicitEndsJson =
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
  "VanDerWaals" : "auto",
  "EndToEndAtoms" : [1, 2]
}
)";

// A three-unit comb polymer whose side chains (two beads) are longer than half the backbone
// (three beads): the topological diameter runs from side-chain tip to side-chain tip (atoms 2, 8),
// but the 'RepeatUnits' inference must pick the backbone ends (atoms 0, 6). Every unit is
// [backbone, side1, side2]; the backbone bead of unit k bonds to the backbone bead of unit k+1.
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

}  // namespace

TEST(END_TO_END_HISTOGRAM, linear_chain_ends_inferred_from_graph_diameter)
{
  ForceField forceField = makeZeroForceField();
  Component component = makeComponent(forceField, "linear-ends", kLinearChainJson);

  ASSERT_TRUE(component.endToEndAtoms.has_value());
  EXPECT_EQ(component.endToEndAtoms->at(0), 0);
  EXPECT_EQ(component.endToEndAtoms->at(1), 3);
}

TEST(END_TO_END_HISTOGRAM, explicit_end_to_end_atoms_win)
{
  ForceField forceField = makeZeroForceField();
  Component component = makeComponent(forceField, "explicit-ends", kExplicitEndsJson);

  ASSERT_TRUE(component.endToEndAtoms.has_value());
  EXPECT_EQ(component.endToEndAtoms->at(0), 1);
  EXPECT_EQ(component.endToEndAtoms->at(1), 2);
}

TEST(END_TO_END_HISTOGRAM, repeat_units_pick_backbone_ends_not_side_chain_tips)
{
  ForceField forceField = makeZeroForceField();
  Component component = makeComponent(forceField, "comb-ends", kCombPolymerJson);

  // The graph diameter would give the side-chain tips (2, 8); the backbone ends are (0, 6).
  ASSERT_TRUE(component.endToEndAtoms.has_value());
  EXPECT_EQ(component.endToEndAtoms->at(0), 0);
  EXPECT_EQ(component.endToEndAtoms->at(1), 6);
}

// The histogram machinery: sampled distances land in the right bins, the probability density
// integrates to one, and the per-bin 95% confidence errors are zero for identical blocks and
// positive when the blocks differ.
TEST(END_TO_END_HISTOGRAM, histogram_normalization_and_error_bars)
{
  ForceField forceField = makeZeroForceField();
  std::vector<Component> components{};
  components.push_back(makeComponent(forceField, "linear-hist", kLinearChainJson));

  constexpr std::size_t numberOfBlocks = 5;
  constexpr std::size_t numberOfBins = 50;
  constexpr double range = 5.0;
  constexpr double delta = range / static_cast<double>(numberOfBins);
  PropertyMoleculeProperties properties(numberOfBlocks, components, numberOfBins, 4.0, 1, 1, range);

  ASSERT_TRUE(properties.endToEndAtomsPerComponent[0].has_value());
  ASSERT_NEAR(properties.deltaEndToEndPerComponent[0], delta, 1e-12);

  // One molecule, end beads 0 and 3 (as inferred): place the ends at a chosen distance.
  std::vector<Atom> moleculeAtoms = components[0].atoms;
  auto sampleAtDistance = [&](double distance, std::size_t block)
  {
    moleculeAtoms[0].position = double3(0.0, 0.0, 0.0);
    moleculeAtoms[3].position = double3(distance, 0.0, 0.0);
    properties.sample(components, {1}, std::span<const Atom>(moleculeAtoms), 0, block);
  };

  // Identical blocks: two samples per block at r = 2.05 and r = 3.55.
  for (std::size_t block = 0; block != numberOfBlocks; ++block)
  {
    sampleAtDistance(2.05, block);
    sampleAtDistance(3.55, block);
  }

  auto [values, average, error] = properties.result(properties.endToEndHistogram, 0, 0, delta, 0.0);

  // Normalization: the density integrates to one.
  double integral = 0.0;
  for (std::size_t bin = 0; bin != numberOfBins; ++bin) integral += average[bin] * delta;
  EXPECT_NEAR(integral, 1.0, 1e-12);

  // Each distance carries half the probability mass in its bin.
  std::size_t binA = static_cast<std::size_t>(2.05 / delta);
  std::size_t binB = static_cast<std::size_t>(3.55 / delta);
  EXPECT_NEAR(average[binA], 0.5 / delta, 1e-12);
  EXPECT_NEAR(average[binB], 0.5 / delta, 1e-12);

  // Identical blocks: the block scatter, and so the confidence-interval error, is exactly zero.
  for (std::size_t bin = 0; bin != numberOfBins; ++bin) EXPECT_EQ(error[bin], 0.0);

  // Skew one block; the error in the affected bins must become positive.
  sampleAtDistance(2.05, 0);
  auto [valuesSkewed, averageSkewed, errorSkewed] =
      properties.result(properties.endToEndHistogram, 0, 0, delta, 0.0);
  EXPECT_GT(errorSkewed[binA], 0.0);
  EXPECT_GT(errorSkewed[binB], 0.0);
}
