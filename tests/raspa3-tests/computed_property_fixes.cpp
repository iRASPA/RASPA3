#include <gtest/gtest.h>

import std;

import double3;
import atom;
import randomnumbers;
import simulationbox;
import property_block_average;
import property_conventional_rdf;
import property_msd;

// Regression tests: the blocks that take part in the statistics of a biased run, the error of the conventional
// radial distribution function where the number of molecules changes, the volume of its shells, and the
// intervals a property cannot be sampled with.

TEST(COMPUTED_PROPERTY_FIXES, a_complete_block_of_low_weight_takes_part_in_the_statistics)
{
  // five blocks of ten samples each; the third spends its time where the bias is high (weights of 1e-3), and
  // the last holds two samples only
  BlockAverage<double> average(5);
  for (std::size_t block = 0; block != 5; ++block)
  {
    const std::size_t numberOfSamples = block == 4 ? 2uz : 10uz;
    for (std::size_t i = 0; i != numberOfSamples; ++i) average.addSample(block, 1.0 + 0.1 * static_cast<double>(block), block == 2 ? 1.0e-3 : 1.0);
  }

  EXPECT_TRUE(average.filled(0));
  EXPECT_TRUE(average.filled(1));
  // by its weight (a hundredth of the first block's) it was left out
  EXPECT_TRUE(average.filled(2));
  EXPECT_TRUE(average.filled(3));
  // a fifth of the samples of the first block: the block of a run that stopped part-way through it
  EXPECT_FALSE(average.filled(4));

  // a block without weight has no average, whatever its number of samples
  BlockAverage<double> weightless(3);
  weightless.addSample(0, 1.0, 1.0);
  weightless.addSample(1, 1.0, 0.0);
  EXPECT_FALSE(weightless.filled(1));
  EXPECT_FALSE(weightless.filled(2));
}

TEST(COMPUTED_PROPERTY_FIXES, the_error_of_the_rdf_of_an_ideal_gas_does_not_follow_the_number_of_molecules)
{
  // an ideal gas whose number of molecules differs from block to block, as in an open ensemble: g(r) is one in
  // every block, and the spread of the number of pairs (a factor of nine here) is no error of it
  const std::size_t numberOfBlocks = 5;
  const std::size_t numberOfBins = 10;
  const SimulationBox box(20.0, 20.0, 20.0);
  PropertyConventionalRadialDistributionFunction rdf(numberOfBlocks, 1, numberOfBins, 8.0, 1, std::nullopt);

  RandomNumber random(12);
  for (std::size_t block = 0; block != numberOfBlocks; ++block)
  {
    const std::size_t numberOfMolecules = (block % 2 == 0) ? 20uz : 60uz;
    for (std::size_t sample = 0; sample != 400; ++sample)
    {
      std::vector<Atom> atoms;
      for (std::size_t i = 0; i != numberOfMolecules; ++i)
      {
        const double3 position{20.0 * random.uniform(), 20.0 * random.uniform(), 20.0 * random.uniform()};
        atoms.push_back(Atom(position, 0.0, 1.0, static_cast<std::uint32_t>(i), 0, 0, false, false));
      }
      rdf.sample(box, std::span<Atom>{}, atoms, 0, block, 1.0);
    }
  }

  const auto results = rdf.result();
  const auto& [r, g, error] = results[0][0];
  ASSERT_EQ(g.size(), numberOfBins);
  for (std::size_t bin = 2; bin != numberOfBins; ++bin)
  {
    EXPECT_NEAR(g[bin], 1.0, 0.05) << "bin " << bin;
    // (normalized by the pairs of the whole run the blocks lay at 0.2 and 1.8, an error of about one)
    EXPECT_LT(error[bin], 0.1) << "bin " << bin;
  }

  // the first bin is a ball, of the volume (4 pi / 3) dr^3: with the mid-point value pi dr^3 in its place g(r)
  // there came out a third too large. (It holds few pairs: all bins up to it are taken together.)
  double weighted = 0.0;
  double volume = 0.0;
  for (std::size_t bin = 0; bin != 3; ++bin)
  {
    const double b = static_cast<double>(bin);
    const double shell = b * b + b + 1.0 / 3.0;
    weighted += g[bin] * shell;
    volume += shell;
  }
  EXPECT_NEAR(weighted / volume, 1.0, 0.05);
}

TEST(COMPUTED_PROPERTY_FIXES, an_interval_of_zero_and_a_single_block_element_are_refused)
{
  // the order-n scheme of the mean-squared displacement never ends with one element per block, and divides by
  // zero with none
  EXPECT_THROW(PropertyMeanSquaredDisplacement(10, std::nullopt, 1), std::runtime_error);
  EXPECT_THROW(PropertyMeanSquaredDisplacement(10, std::nullopt, 0), std::runtime_error);
  EXPECT_THROW(PropertyMeanSquaredDisplacement(0, std::nullopt, 25), std::runtime_error);
  EXPECT_THROW(PropertyMeanSquaredDisplacement(10, 0uz, 25), std::runtime_error);
  EXPECT_NO_THROW(PropertyMeanSquaredDisplacement(10, 5000uz, 25));

  EXPECT_THROW(PropertyConventionalRadialDistributionFunction(5, 1, 10, 8.0, 0, std::nullopt), std::runtime_error);
  EXPECT_THROW(PropertyConventionalRadialDistributionFunction(5, 1, 10, 8.0, 1, 0uz), std::runtime_error);
}
