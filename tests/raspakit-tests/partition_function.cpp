#include <gtest/gtest.h>

import std;

import janaf;
import nasa_polynomials;
import partition_function;

// The PartitionFunction facade uses the NASA polynomials (Burcat) as the primary source and
// the JANAF tables as fallback.
TEST(partition_function, nasa_polynomials_is_the_primary_source)
{
  // species present in both sources must be served by the primary source
  for (std::string_view name : {"H2O", "CO2", "CH4", "N2"})
  {
    EXPECT_EQ(PartitionFunction::logPartitionFunction(name, 600.0),
              NASAPolynomials::logPartitionFunction(name, 600.0))
        << name;
  }

  // species only available in the primary source
  EXPECT_TRUE(PartitionFunction::contains("propene"));
  EXPECT_EQ(PartitionFunction::logPartitionFunction("propene", 450.0),
            NASAPolynomials::logPartitionFunction("propene", 450.0));
}

TEST(partition_function, janaf_fallback_below_nasa_temperature_range)
{
  // the NASA polynomials start at 200 K, the JANAF tables at 100 K
  EXPECT_THROW(std::ignore = NASAPolynomials::logPartitionFunction("N2", 150.0), std::runtime_error);
  EXPECT_EQ(PartitionFunction::logPartitionFunction("N2", 150.0), JANAF::logPartitionFunction("N2", 150.0));
}

TEST(partition_function, janaf_only_source)
{
  EXPECT_TRUE(PartitionFunction::contains("N2", PartitionFunction::Source::JANAF));
  EXPECT_FALSE(PartitionFunction::contains("propene", PartitionFunction::Source::JANAF));
  EXPECT_NEAR(PartitionFunction::logPartitionFunction("N2", 600.0, PartitionFunction::Source::JANAF),
              JANAF::logPartitionFunction("N2", 600.0), 1e-10);
  EXPECT_THROW(std::ignore = PartitionFunction::logPartitionFunction("propene", 450.0, PartitionFunction::Source::JANAF),
               std::runtime_error);
}

TEST(partition_function, unknown_species)
{
  EXPECT_FALSE(PartitionFunction::contains("unobtainium"));
  EXPECT_THROW(std::ignore = PartitionFunction::logPartitionFunction("unobtainium", 300.0), std::runtime_error);
}
