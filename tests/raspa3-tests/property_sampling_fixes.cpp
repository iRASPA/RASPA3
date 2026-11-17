#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import units;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;

// Regression tests of the property-sampling fixes: the lambda histograms count a sample with the weight of the
// other fractional molecules, and the tail correction between two components is booked like their pair energies.

// The histogram along one lambda and dU/dlambda per lambda are kept under the bias of that histogram, which the
// results take out. The biases of the other fractional molecules in the box have to go with the sample: it
// counts with the weight of the configuration over the factor of this histogram. Counted with one, the
// distribution was that of a system whose other fractional molecules are biased.
TEST(PROPERTY_SAMPLING, a_lambda_histogram_counts_a_sample_with_the_weight_of_the_other_fractional_molecules)
{
  PropertyLambdaProbabilityHistogram histogram(5, 3);
  histogram.biasFactor = {-2.0, 0.5, 1.0};
  histogram.setCurrentBin(0);
  const double own = histogram.weight();
  ASSERT_DOUBLE_EQ(own, std::exp(2.0));

  // alone in the box: the weight of the configuration is the factor of this histogram, and the sample counts one
  histogram.sampleHistogram(0, 1.0, 10.0, true, own);
  EXPECT_DOUBLE_EQ(histogram.bookKeepingLambda[0][0], 1.0);
  EXPECT_DOUBLE_EQ(histogram.bookKeepingDUdlambda[0][0].first, 10.0);
  EXPECT_DOUBLE_EQ(histogram.bookKeepingDUdlambda[0][0].second, 1.0);

  // beside another fractional molecule whose bias gives the configuration a factor of 1/4
  histogram.sampleHistogram(1, 1.0, 10.0, true, 0.25 * own);
  EXPECT_NEAR(histogram.bookKeepingLambda[1][0], 0.25, 1.0e-14);
  EXPECT_NEAR(histogram.bookKeepingDUdlambda[1][0].first, 2.5, 1.0e-13);
  EXPECT_NEAR(histogram.bookKeepingDUdlambda[1][0].second, 0.25, 1.0e-14);
  // the average dU/dlambda at that lambda is the weighted one
  EXPECT_NEAR(histogram.bookKeepingDUdlambda[1][0].first / histogram.bookKeepingDUdlambda[1][0].second, 10.0, 1.0e-12);

  // a histogram whose bias is not in the weight (the lambda of an (N, lambda) transition-matrix chain) takes the
  // weight as it is
  histogram.sampleHistogram(2, 1.0, 10.0, true, 0.25, false);
  EXPECT_DOUBLE_EQ(histogram.bookKeepingLambda[2][0], 0.25);

  // a box without the fractional molecule does not count
  histogram.sampleHistogram(3, 1.0, 10.0, false, own);
  EXPECT_DOUBLE_EQ(histogram.bookKeepingLambda[3][0], 0.0);
}

// The pair energies of two components are booked half on (A, B) and half on (B, A), and the report of the pair
// adds the two. The tail correction was booked whole on (A, B): the report doubled it.
TEST(PROPERTY_SAMPLING, the_tail_correction_between_two_components_is_booked_like_their_pair_energy)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, false, true, false);
  const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("CH4"));
  const auto methane = [&](std::string name, std::uint8_t componentId)
  {
    return Component(forceField, std::move(name), 190.564, 4599200, 0.01142,
                     {Atom({0, 0, 0}, 0.0, 1.0, 0, type, componentId, false, false)}, {}, {}, 5, 21, {}, std::nullopt,
                     false);
  };
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
                {methane("first", 0), methane("second", 1)}, {}, {7, 11}, 5);

  const RunningEnergy total = system.computeTotalEnergies();
  ASSERT_NE(total.tail, 0.0);
  EnergyStatus status = system.computeMolecularPressure().first;

  const double firstFirst = status.componentEnergy(0, 0).VanDerWaalsTailCorrection.energy;
  const double secondSecond = status.componentEnergy(1, 1).VanDerWaalsTailCorrection.energy;
  const double firstSecond = status.componentEnergy(0, 1).VanDerWaalsTailCorrection.energy;
  const double secondFirst = status.componentEnergy(1, 0).VanDerWaalsTailCorrection.energy;
  EXPECT_NEAR(firstSecond, secondFirst, 1.0e-12 * std::abs(firstSecond));
  EXPECT_NEAR(firstFirst + secondSecond + firstSecond + secondFirst, total.tail, 1.0e-10 * std::abs(total.tail));
  // identical species: the tail energies go as the products of the counts, 7 x 7 : 11 x 11 : 2 x 7 x 11
  EXPECT_NEAR((firstSecond + secondFirst) / firstFirst, 2.0 * 7.0 * 11.0 / (7.0 * 7.0), 1.0e-10);
  EXPECT_NEAR(secondSecond / firstFirst, 121.0 / 49.0, 1.0e-10);
}
