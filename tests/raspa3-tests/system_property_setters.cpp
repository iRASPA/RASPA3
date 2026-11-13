#include <gtest/gtest.h>

import std;

import forcefield;
import simulationbox;
import system;
import property_conventional_rdf;
import property_rdf;

// The property setters of a System (the Python 'set_property_*'): the RDFs are sized for the blocks of the system
// and keep the range asked for, and None clears a property.

namespace
{
System makeArgonBox(std::size_t numberOfBlocks)
{
  const ForceField forceField({{"Ar", false, 39.948, 0.0, 0.0, 18, false}}, {{119.8, 3.4}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {}, {}, {}, numberOfBlocks);
}
}  // namespace

// Five blocks and a range of 12 were fixed in the setters: a simulation of more blocks wrote past the block sums.
TEST(SYSTEM_PROPERTY_SETTERS, rdfs_get_the_blocks_of_the_system_and_their_own_range)
{
  System system = makeArgonBox(8);

  system.setPropertyConventionalRDF(PropertyConventionalRadialDistributionFunction(64, 9.0, 1, 100));
  ASSERT_TRUE(system.propertyConventionalRadialDistributionFunction.has_value());
  EXPECT_EQ(system.propertyConventionalRadialDistributionFunction->histogram.numberOfBlocks, 8uz);
  EXPECT_EQ(system.propertyConventionalRadialDistributionFunction->range, 9.0);

  system.setPropertyRDF(PropertyRadialDistributionFunction(1, 1, 64, 10.0, 1, 100));
  ASSERT_TRUE(system.propertyRadialDistributionFunction.has_value());
  EXPECT_EQ(system.propertyRadialDistributionFunction->histogram.numberOfBlocks, 8uz);

  system.setPropertyConventionalRDF(std::nullopt);
  system.setPropertyRDF(std::nullopt);
  EXPECT_FALSE(system.propertyConventionalRadialDistributionFunction.has_value());
  EXPECT_FALSE(system.propertyRadialDistributionFunction.has_value());
}
