#include <gtest/gtest.h>

import std;

import double3;
import simulationbox;
import voronoi_accessibility;
import exact_pore_size_distribution;

// The pore-size distribution as a curve: the sweep over probe sizes, the spikes at the sizes that carry
// volume of their own, and the two ways the whole of the void is accounted for.
//
// A simple-cubic lattice of one atom is the case to test this on, because the answer is a formula. The
// Apollonius cell of the atom is the cube, so the deepest point of the void is the corner of that cube and
// the largest sphere that fits anywhere is sqrt(3) a - 2 r across. The distribution therefore has to end
// exactly there, and since that sphere carries volume with it, it has to end with a spike and not with a
// slide to zero. Nothing else about the curve is known in closed form, but everything about it is
// constrained: the volume falls, the distribution does not go negative, and the continuous part together
// with the spikes has to account for the whole void, since every point of it has some pore size.

namespace
{

struct Lattice
{
  SimulationBox box{6.0, 6.0, 6.0};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
};


Lattice simpleCubic(double edge, double radius)
{
  Lattice lattice;
  lattice.box = SimulationBox(edge, edge, edge);
  lattice.fractionalPositions = {double3(0.0, 0.0, 0.0)};
  lattice.radii = {radius};
  return lattice;
}


PoreSizeDistributionCurve curveOf(const Lattice& lattice, double maximumDiameter, std::size_t bins)
{
  auto build = [&](double probeRadius)
  {
    return VoronoiAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probeRadius);
  };
  return exactPoreSizeDistribution(build, lattice.box.volume, maximumDiameter, bins, 1);
}

}  // namespace


TEST(exact_pore_size_distribution, a_simple_cubic_lattice_ends_at_its_corner)
{
  const double edge = 6.0;
  const double radius = 1.5;
  const double corner = std::sqrt(3.0) * edge - 2.0 * radius;  // the largest sphere that fits anywhere

  PoreSizeDistributionCurve curve = curveOf(simpleCubic(edge, radius), 9.0, 90);

  ASSERT_FALSE(curve.spikes.empty());
  const PoreSizeSpike& last = curve.spikes.back();

  // Cornered by bisection, so it is placed to its own bracket and not to the spacing of the rows.
  EXPECT_NEAR(last.diameter, corner, 1.0e-4) << "bracket " << last.bracket;
  EXPECT_LT(last.bracket, 1.0e-4);

  // The largest sphere carries volume of its own, and here it is most of the void: the corners of the cube
  // are where nearly all of the room is at this density.
  EXPECT_GT(last.weight, 0.1);

  // Nothing is left above it, the range having been asked past it.
  EXPECT_LT(curve.truncatedWeight, 1.0e-9);
  EXPECT_NEAR(curve.largestDiameter, corner, 1.0e-4);
}


TEST(exact_pore_size_distribution, the_continuous_part_and_the_spikes_account_for_the_void)
{
  const double edge = 6.0;
  const double radius = 1.5;

  PoreSizeDistributionCurve curve = curveOf(simpleCubic(edge, radius), 9.0, 90);

  // Every point of the void has some pore size, so the two have to make up the whole of it. They are
  // computed by routes with no term in common, the one from the derivative and the other from the volume
  // the derivative fails to account for, so this is a real comparison and not an identity. What it is out
  // by is the truncation error of the trapezium rule over the rows.
  EXPECT_NEAR(curve.integral + curve.singularWeight, 1.0, 0.02)
      << "continuous " << curve.integral << ", spikes " << curve.singularWeight;

  // Refining the rows has to bring it nearer, which is what says the residual is the quadrature and not a
  // disagreement between the volume and its derivative.
  PoreSizeDistributionCurve finer = curveOf(simpleCubic(edge, radius), 9.0, 360);
  EXPECT_LT(std::abs(finer.integral + finer.singularWeight - 1.0),
            std::abs(curve.integral + curve.singularWeight - 1.0));

  // The spikes do not move under refinement either, being cornered to their own brackets rather than
  // assigned to a row.
  ASSERT_FALSE(finer.spikes.empty());
  EXPECT_NEAR(finer.spikes.back().diameter, curve.spikes.back().diameter, 1.0e-4);
  EXPECT_NEAR(finer.spikes.back().weight, curve.spikes.back().weight, 1.0e-3);
}


TEST(exact_pore_size_distribution, the_curve_falls_and_the_distribution_does_not_go_negative)
{
  PoreSizeDistributionCurve curve = curveOf(simpleCubic(6.0, 1.5), 9.0, 90);

  ASSERT_GT(curve.points.size(), 10uz);
  EXPECT_NEAR(curve.points.front().cumulative, 1.0, 1.0e-3);

  double previous = std::numeric_limits<double>::max();
  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    EXPECT_LE(point.cumulative, previous + 1.0e-9) << "at d = " << point.diameter;
    EXPECT_GE(point.cumulative, -1.0e-9) << "at d = " << point.diameter;
    EXPECT_GE(point.distribution, -1.0e-9) << "at d = " << point.diameter;

    // The three shares of the distribution are a division of it and not a second account of it.
    EXPECT_NEAR(point.accessible + point.inaccessible + point.undecided, point.distribution, 1.0e-9);
    previous = point.cumulative;
  }
  EXPECT_NEAR(curve.points.back().cumulative, 0.0, 1.0e-9);
}


// A lattice whose atoms are far enough apart that a probe of the sizes swept never has to squeeze past
// two of them at once. The excluded surface is then convex everywhere, nothing of it is reentrant, and
// the distribution is zero however much room there is between the atoms: the whole of the void is opened
// by every probe that has been asked about, so there is no pore size below the size at which it stops.
TEST(exact_pore_size_distribution, a_convex_void_has_no_distribution_at_all)
{
  const double edge = 12.0;
  const double radius = 1.0;
  // The inflated spheres of two neighbours touch at a probe radius of edge/2 - radius = 5, so up to a
  // diameter of 4 the boundary is a set of disjoint spheres.
  PoreSizeDistributionCurve curve = curveOf(simpleCubic(edge, radius), 4.0, 40);

  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    EXPECT_NEAR(point.distribution, 0.0, 1.0e-9) << "at d = " << point.diameter;
    EXPECT_NEAR(point.cumulative, 1.0, 1.0e-9) << "at d = " << point.diameter;
    EXPECT_EQ(point.numberOfArcs, 0uz);
  }
  EXPECT_TRUE(curve.spikes.empty());
  EXPECT_NEAR(curve.integral, 0.0, 1.0e-9);

  // The whole of the void is still open at the end of the range, which is the range stopping short and
  // is reported as such rather than as a distribution that failed to integrate to one.
  EXPECT_NEAR(curve.truncatedWeight, 1.0, 1.0e-9);
}
