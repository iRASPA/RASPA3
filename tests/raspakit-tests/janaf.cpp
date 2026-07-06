#include <gtest/gtest.h>

import std;

import units;
import janaf;

// Reference values from Table A5 of the PhD thesis of A. Rahbari (TU Delft, 2020), JANAF
// columns: ideal-gas partition functions q/V in [Å⁻³] for nitrogen, hydrogen, and ammonia,
// referenced to the dissociated atoms in their electronic ground state.
TEST(janaf, partition_functions_thesis_table_A5)
{
  const double tolerance = 0.1;  // in ln(q/V); table values are given to 3 significant digits

  EXPECT_NEAR(JANAF::logPartitionFunction("N2", 573.0), std::log(2.67e90), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("N2", 673.0), std::log(7.04e77), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("N2", 773.0), std::log(3.52e68), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("N2", 873.0), std::log(2.48e61), tolerance);

  EXPECT_NEAR(JANAF::logPartitionFunction("H2", 573.0), std::log(6.53e40), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("H2", 673.0), std::log(1.36e35), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("H2", 773.0), std::log(8.79e30), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("H2", 873.0), std::log(5.38e27), tolerance);

  EXPECT_NEAR(JANAF::logPartitionFunction("NH3", 573.0), std::log(1.46e110), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("NH3", 673.0), std::log(5.26e94), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("NH3", 773.0), std::log(2.06e83), tolerance);
  EXPECT_NEAR(JANAF::logPartitionFunction("NH3", 873.0), std::log(3.58e74), tolerance);
}

// Atomization energies D₀ at 0 K computed from the JANAF enthalpies of formation (Eq. A78),
// compared with the experimental values in Table A4 of the thesis (McQuarrie).
TEST(janaf, atomization_energies_thesis_table_A4)
{
  EXPECT_NEAR(JANAF::atomizationEnergy("N2"), 941.6, 0.5);
  EXPECT_NEAR(JANAF::atomizationEnergy("H2"), 432.1, 0.5);
  EXPECT_NEAR(JANAF::atomizationEnergy("NH3"), 1158.0, 0.5);
}

TEST(janaf, dissociated_atom_and_ground_state_references_are_consistent)
{
  const double temperature = 700.0;
  for (std::string_view name : JANAF::availableMolecules())
  {
    const double expected = JANAF::logPartitionFunctionGroundStateReference(name, temperature) +
                            JANAF::atomizationEnergy(name) * 1000.0 / (Units::MolarGasConstant * temperature);
    EXPECT_NEAR(JANAF::logPartitionFunction(name, temperature), expected, 1e-10);
  }
}

TEST(janaf, name_matching)
{
  EXPECT_TRUE(JANAF::contains("CO2"));
  EXPECT_TRUE(JANAF::contains("co2"));
  EXPECT_TRUE(JANAF::contains("carbon dioxide"));
  EXPECT_TRUE(JANAF::contains("Carbon-Dioxide"));
  EXPECT_TRUE(JANAF::contains("ammonia"));
  EXPECT_TRUE(JANAF::contains("ethylene"));
  EXPECT_FALSE(JANAF::contains("unobtainium"));

  EXPECT_EQ(JANAF::logPartitionFunction("water", 600.0), JANAF::logPartitionFunction("H2O", 600.0));

  EXPECT_THROW(std::ignore = JANAF::logPartitionFunction("unobtainium", 600.0), std::runtime_error);
  EXPECT_THROW(std::ignore = JANAF::logPartitionFunction("N2", 50.0), std::runtime_error);
  EXPECT_THROW(std::ignore = JANAF::logPartitionFunction("N2", 5000.0), std::runtime_error);
}
