#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import crystal;
import pair_interactions;
import exact_pore_size_distribution;
import voronoi_pore_size_distribution;
import apollonius_pore_size_distribution;

// The hybrid pore-size distribution with a floor short of the blocking probe, through the two drivers.
//
// Atoms of a simple-cubic lattice overlap along the edges of the cell, so there are crevices between them that
// lose volume from the smallest diameters up, and a floor has something to flatten. The windows are wider than
// the probe, so nothing is sealed and all of the difference between the curves is the floor's.

namespace
{

constexpr double edge = 6.0;
constexpr double atomRadius = 3.3;
constexpr double probeRadius = 0.5;
constexpr double floorRadius = 0.25;
constexpr double floorDiameter = 2.0 * floorRadius;
constexpr double maximumDiameter = 6.0;
constexpr std::size_t bins = 60;

PairInteractions atomAndProbe()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"O", "He"};
  interactions.charges = {0.0, 0.0};
  const double cross = atomRadius + probeRadius;
  interactions.parameters = {PairParameters{2.0 * atomRadius, 100.0, 0.0}, PairParameters{cross, 100.0, 0.0},
                             PairParameters{cross, 100.0, 0.0}, PairParameters{2.0 * probeRadius, 100.0, 0.0}};
  return interactions;
}

Crystal overlappingLattice()
{
  Crystal framework;
  framework.name = "hybrid_floor_test";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{double3(0.0, 0.0, 0.0), 0.0, 0});
  framework.fractionalPositions.push_back(double3(0.0, 0.0, 0.0));
  return framework;
}

// Removes the reports the drivers leave behind, whichever way the test ends.
struct ReportCleanup
{
  ~ReportCleanup()
  {
    for (const char* diagram : {"voronoi", "apollonius"})
    {
      for (const char* kind : {"psd", "psd-peaks"})
      {
        std::error_code ignored;
        std::filesystem::remove(std::format("hybrid_floor_test.{}.{}.txt", diagram, kind), ignored);
      }
    }
  }
};

template <typename Driver>
PoreSizeDistributionCurve sweep(double floor, bool peaksOnly = false)
{
  Driver driver;
  driver.run(atomAndProbe(), overlappingLattice(), "He", maximumDiameter, bins, 1, floor, peaksOnly);
  return driver.curve;
}

template <typename Driver>
void expectFlatBelowTheFloorAndBareAboveIt(const char* diagram)
{
  ReportCleanup cleanup;
  PoreSizeDistributionCurve curve = sweep<Driver>(floorRadius);

  EXPECT_DOUBLE_EQ(curve.probeRadius, probeRadius);
  EXPECT_DOUBLE_EQ(curve.floorRadius, floorRadius);
  EXPECT_TRUE(hybridCurve(curve));

  // Normalised by what is left at the floor, which the crevices have already lost from.
  EXPECT_LT(curve.probeAccessibleVolume, 0.999 * curve.voidVolume);
  EXPECT_GT(curve.probeAccessibleVolume, 0.0);
  EXPECT_NEAR(curve.probeAccessibleIntegral + curve.probeAccessibleSingularWeight, 1.0, 1.0e-2);

  bool bareCurveFallsBelowTheFloor = false;
  std::size_t rowsBelow = 0;
  std::size_t rowsAbove = 0;
  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    if (point.poreVolume <= 0.0) continue;

    if (point.diameter < floorDiameter)
    {
      ++rowsBelow;
      EXPECT_DOUBLE_EQ(point.probeAccessibleCumulative, 1.0) << "d " << point.diameter;
      EXPECT_DOUBLE_EQ(point.probeAccessibleDistribution, 0.0) << "d " << point.diameter;
      EXPECT_DOUBLE_EQ(point.probeAccessiblePoreVolume, curve.probeAccessibleVolume) << "d " << point.diameter;
      if (point.cumulative < 1.0 - 1.0e-6 && point.distribution > 0.0) bareCurveFallsBelowTheFloor = true;
    }
    else
    {
      // Nothing is blocked in this sweep, so above the floor the two curves are the one volume over two norms.
      ++rowsAbove;
      EXPECT_NEAR(point.probeAccessiblePoreVolume, point.poreVolume, 1.0e-9) << "d " << point.diameter;
      EXPECT_NEAR(point.probeAccessibleCumulative, point.poreVolume / curve.probeAccessibleVolume, 1.0e-9)
          << "d " << point.diameter;
    }
  }
  EXPECT_GT(rowsBelow, 0u);
  EXPECT_GT(rowsAbove, 0u);
  EXPECT_TRUE(bareCurveFallsBelowTheFloor) << "the columns of the bare curve are to keep what the floor flattens";

  std::ifstream report(std::format("hybrid_floor_test.{}.psd.txt", diagram));
  std::string text((std::istreambuf_iterator<char>(report)), std::istreambuf_iterator<char>());
  EXPECT_FALSE(text.empty());
  EXPECT_EQ(text.find("no diameter floor"), std::string::npos);
  EXPECT_EQ(text.find("mirror"), std::string::npos);
}

// With no floor the bare curve is the primary one, by way of the void volume rather than a second measure of it.
template <typename Driver>
void expectNoFloorToMirrorTheBareCurve()
{
  ReportCleanup cleanup;
  PoreSizeDistributionCurve curve = sweep<Driver>(0.0);

  EXPECT_TRUE(hybridCurve(curve));
  EXPECT_DOUBLE_EQ(curve.probeAccessibleVolume, curve.voidVolume);
  EXPECT_DOUBLE_EQ(curve.probeAccessibleIntegral, curve.integral);
  EXPECT_DOUBLE_EQ(curve.probeAccessibleSingularWeight, curve.singularWeight);
  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    EXPECT_DOUBLE_EQ(point.probeAccessibleCumulative, point.cumulative);
    EXPECT_DOUBLE_EQ(point.probeAccessibleDistribution, point.distribution);
    EXPECT_DOUBLE_EQ(point.probeAccessiblePoreVolume, point.poreVolume);
  }
}

// The peaks route takes the same promotion: its weights are shares of the volume at the floor.
template <typename Driver>
void expectPeaksNormalisedAtTheFloor()
{
  ReportCleanup cleanup;
  PoreSizeDistributionCurve floored = sweep<Driver>(floorRadius, true);

  EXPECT_LT(floored.probeAccessibleVolume, 0.999 * floored.voidVolume);
  ASSERT_FALSE(floored.probeAccessibleSpikes.empty());
  ASSERT_FALSE(floored.spikes.empty());

  const double ratio = floored.voidVolume / floored.probeAccessibleVolume;
  EXPECT_NEAR(floored.probeAccessibleSpikes.back().diameter, floored.spikes.back().diameter, 1.0e-9);
  EXPECT_NEAR(floored.probeAccessibleSpikes.back().weight, floored.spikes.back().weight * ratio, 1.0e-6);
}

}  // namespace


TEST(hybrid_pore_size_floor, voronoi_is_flat_below_the_floor)
{
  expectFlatBelowTheFloorAndBareAboveIt<VoronoiPoreSizeDistribution>("voronoi");
}

TEST(hybrid_pore_size_floor, apollonius_is_flat_below_the_floor)
{
  expectFlatBelowTheFloorAndBareAboveIt<ApolloniusPoreSizeDistribution>("apollonius");
}

TEST(hybrid_pore_size_floor, voronoi_without_a_floor_mirrors_the_bare_curve)
{
  expectNoFloorToMirrorTheBareCurve<VoronoiPoreSizeDistribution>();
}

TEST(hybrid_pore_size_floor, apollonius_without_a_floor_mirrors_the_bare_curve)
{
  expectNoFloorToMirrorTheBareCurve<ApolloniusPoreSizeDistribution>();
}

TEST(hybrid_pore_size_floor, voronoi_peaks_are_shares_of_the_volume_at_the_floor)
{
  expectPeaksNormalisedAtTheFloor<VoronoiPoreSizeDistribution>();
}

TEST(hybrid_pore_size_floor, apollonius_peaks_are_shares_of_the_volume_at_the_floor)
{
  expectPeaksNormalisedAtTheFloor<ApolloniusPoreSizeDistribution>();
}
