#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import forcefield;
import simulationbox;

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

// "CutOff": "auto" sets the cut-off to half the box at run time; the potential shifts and the tail corrections
// (and the switching constants) must follow it instead of staying at the value of the force-field file.
TEST(forcefield_options, automatic_cutoff_recomputes_shifts_and_tail_corrections)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  const double automaticCutOff = 10.0 - std::numeric_limits<double>::epsilon();

  ForceField automatic({{"CH4", false, 16.04, 0.0, 0.0, 6, false}, {"O", true, 16.0, 0.0, 0.0, 8, false}},
                       {{158.5, 3.72}, {53.0, 3.3}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                       true, true, false);
  automatic.cutOffFrameworkVDWAutomatic = true;
  automatic.cutOffMoleculeVDWAutomatic = true;
  const double shiftAt12 = automatic.data[0].shift;
  const double tailAt12 = automatic.data[0].tailCorrectionEnergy;

  automatic.initializeAutomaticCutOff(box);
  EXPECT_DOUBLE_EQ(automatic.cutOffMoleculeVDW, automaticCutOff);
  EXPECT_DOUBLE_EQ(automatic.cutOffFrameworkVDW, automaticCutOff);

  const ForceField reference({{"CH4", false, 16.04, 0.0, 0.0, 6, false}, {"O", true, 16.0, 0.0, 0.0, 8, false}},
                             {{158.5, 3.72}, {53.0, 3.3}}, ForceField::MixingRule::Lorentz_Berthelot,
                             automaticCutOff, automaticCutOff, 12.0, true, true, false);
  for (std::size_t pair = 0; pair < 4; ++pair)
  {
    EXPECT_DOUBLE_EQ(automatic.data[pair].shift, reference.data[pair].shift);
    EXPECT_DOUBLE_EQ(automatic.data[pair].tailCorrectionEnergy, reference.data[pair].tailCorrectionEnergy);
    EXPECT_DOUBLE_EQ(automatic.data[pair].tailCorrectionPressure, reference.data[pair].tailCorrectionPressure);
  }
  EXPECT_NE(automatic.data[0].shift, shiftAt12);
  EXPECT_NE(automatic.data[0].tailCorrectionEnergy, tailAt12);

  // explicit cut-offs are left alone
  ForceField explicitCutOff = reference;
  explicitCutOff.initializeAutomaticCutOff(SimulationBox(30.0, 30.0, 30.0));
  EXPECT_DOUBLE_EQ(explicitCutOff.cutOffMoleculeVDW, automaticCutOff);
  EXPECT_DOUBLE_EQ(explicitCutOff.data[0].shift, reference.data[0].shift);
}

// Switching a shifted force field to plain truncation has to drop the shifts it computed before.
TEST(forcefield_options, truncating_a_shifted_force_field_clears_the_shifts)
{
  std::optional<ForceField> forceField = parse("");
  ASSERT_TRUE(forceField.has_value());
  ASSERT_NE(forceField->data[0].shift, 0.0);

  forceField->setTruncationMethod(ForceField::TruncationMethod::Truncated);
  EXPECT_DOUBLE_EQ(forceField->data[0].shift, 0.0);
}

// An integer JSON number is a valid value for a floating-point option.
TEST(forcefield_options, integer_valued_floating_point_options_are_read)
{
  std::optional<ForceField> forceField = parse(R"("EwaldPrecision": 1)");
  ASSERT_TRUE(forceField.has_value());
  EXPECT_DOUBLE_EQ(forceField->EwaldPrecision, 1.0);
}
