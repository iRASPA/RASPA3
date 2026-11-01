#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import forcefield;

// Force-field JSON options for the CBMC internal samplers (ring-closure Monte-Carlo tuning knobs).

namespace
{
std::string forceFieldJson(std::string_view extraOptions)
{
  return std::format(R"({{
  "PseudoAtoms": [
    {{
      "name": "CH4", "framework": false, "print_to_output": true,
      "element": "C", "print_as": "C", "mass": 16.04246, "charge": 0.0
    }}
  ],
  "SelfInteractions": [
    {{
      "name": "CH4", "type": "lennard-jones", "parameters": [158.5, 3.72]
    }}
  ],
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "shifted",
  "TailCorrections": false,
  "CutOffVDW": 11.8{}{}
}})",
                     extraOptions.empty() ? "" : ",\n  ", extraOptions);
}

std::optional<ForceField> parse(std::string_view extraOptions)
{
  TemporaryDirectory dir;
  dir.write("force_field.json", forceFieldJson(extraOptions));
  return ForceField::readForceField(dir.path().string(), "force_field.json");
}
}  // namespace

TEST(forcefield_options, ring_mc_knobs_have_defaults)
{
  std::optional<ForceField> forceField = parse("");
  ASSERT_TRUE(forceField.has_value());
  EXPECT_EQ(forceField->settings.numberOfTrialMovesPerOpenBead, 150uz);
  EXPECT_DOUBLE_EQ(forceField->settings.cbmcRingCrankshaftProbability, 0.2);
  EXPECT_DOUBLE_EQ(forceField->settings.cbmcRingTiltProbability, 0.25);
}

TEST(forcefield_options, ring_mc_knobs_are_parsed)
{
  std::optional<ForceField> forceField = parse(
      R"("NumberOfTrialMovesPerOpenBead": 60,
  "CBMCRingCrankshaftProbability": 0.35,
  "CBMCRingTiltProbability": 0.1)");
  ASSERT_TRUE(forceField.has_value());
  EXPECT_EQ(forceField->settings.numberOfTrialMovesPerOpenBead, 60uz);
  EXPECT_DOUBLE_EQ(forceField->settings.cbmcRingCrankshaftProbability, 0.35);
  EXPECT_DOUBLE_EQ(forceField->settings.cbmcRingTiltProbability, 0.1);
}

TEST(forcefield_options, ring_mc_probabilities_outside_unit_interval_throw)
{
  EXPECT_THROW(parse(R"("CBMCRingCrankshaftProbability": 1.5)"), std::runtime_error);
  EXPECT_THROW(parse(R"("CBMCRingTiltProbability": -0.1)"), std::runtime_error);
}

// operator== must compare the pseudo-atoms of the two force fields (it used to compare iterators into the
// same container, so force fields differing only in their pseudo-atoms compared equal).
TEST(forcefield_options, equality_compares_the_pseudo_atoms)
{
  const ForceField a({{"CH4", false, 16.04, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                     ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  ForceField b = a;
  EXPECT_TRUE(a == b);

  b.pseudoAtoms[0].mass = 15.0;
  EXPECT_FALSE(a == b);

  ForceField c = a;
  c.pseudoAtoms[0].name = "CH3";
  EXPECT_FALSE(a == c);
}
