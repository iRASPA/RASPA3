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
import property_conventional_rdf;
import property_density_grid;
import int3;
import framework;

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

// The structural samplers take the weight of the configuration, as the thermodynamic averages do: under a biased
// walk a configuration counts with its weight (they counted every configuration once, so their results were those
// of the biased walk and disagreed with the reweighted loadings).
TEST(PROPERTY_SAMPLING, the_conventional_rdf_counts_a_configuration_with_its_weight)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  auto pairAt = [](double distance)
  {
    std::vector<Atom> atoms(2);
    atoms[0].position = double3(5.0, 5.0, 5.0);
    atoms[1].position = double3(5.0 + distance, 5.0, 5.0);
    atoms[1].moleculeId = 1;
    return atoms;
  };
  std::vector<Atom> near = pairAt(2.1);
  std::vector<Atom> far = pairAt(4.1);
  std::vector<Atom> noFramework{};

  PropertyConventionalRadialDistributionFunction rdf(1, 1, 40, 10.0, 1, std::nullopt);
  rdf.sample(box, noFramework, near, 0, 0, 3.0);
  rdf.sample(box, noFramework, far, 1, 0, 1.0);

  // the pair at 2.1 counts three times the pair at 4.1: 3/4 and 1/4 of the weighted pairs, over the shell volume
  const auto results = rdf.result();
  const std::vector<double>& g = std::get<1>(results[0][0]);
  const double deltaR = 0.25;
  auto expected = [&](std::size_t bin, double fraction)
  {
    const double shell = 4.0 * std::numbers::pi * deltaR * deltaR * deltaR * (static_cast<double>(bin) + 0.5) *
                         (static_cast<double>(bin) + 0.5);
    return fraction * box.volume / shell;
  };
  EXPECT_NEAR(g[8], expected(8, 0.75), 1e-9 * expected(8, 0.75));
  EXPECT_NEAR(g[16], expected(16, 0.25), 1e-9 * expected(16, 0.25));

  // a configuration with a fractional molecule only is no sample
  std::vector<Atom> fractionalOnly(1);
  fractionalOnly[0].isFractional = true;
  const double volumeSamples = rdf.volumeSampleCount;
  rdf.sample(box, noFramework, fractionalOnly, 2, 0, 1.0);
  EXPECT_DOUBLE_EQ(rdf.volumeSampleCount, volumeSamples);
}

TEST(PROPERTY_SAMPLING, the_density_grid_counts_a_configuration_with_its_weight)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  PropertyDensityGrid grid(0, 1, int3(4, 4, 4), 1, 1, {}, PropertyDensityGrid::Normalization::NumberDensity,
                           PropertyDensityGrid::Binning::Standard);
  std::vector<Atom> one(1);
  one[0].position = double3(1.0, 1.0, 1.0);
  std::vector<Atom> three(3);
  for (Atom& atom : three) atom.position = double3(11.0, 11.0, 11.0);

  grid.sample(std::nullopt, box, one, 0, 3.0);
  grid.sample(std::nullopt, box, three, 1, 1.0);

  // the weighted number of atoms over the weight: (3 * 1 + 1 * 3) / 4
  const double atoms = std::accumulate(grid.grid_cell.begin(), grid.grid_cell.end(), 0.0);
  EXPECT_DOUBLE_EQ(atoms / grid.weightSum, 1.5);
  EXPECT_DOUBLE_EQ(grid.volumeSum / grid.weightSum, box.volume);
}
