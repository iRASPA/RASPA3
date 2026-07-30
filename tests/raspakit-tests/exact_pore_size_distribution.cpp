#include <gtest/gtest.h>

import std;

import double3;
import simulationbox;
import pore_accessibility;
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


// Twelve atoms at the vertices of an icosahedron, overlapping over the edges and over the faces as well, so
// that they wall in a small pore nothing can enter: a cage in the middle of a cell that is otherwise open.
// This is a sodalite cage in a zeolite in miniature, and it is there to be left out.
Lattice cageInAnOpenCell(double edge, double shellRadius, double atomRadius)
{
  Lattice lattice;
  lattice.box = SimulationBox(edge, edge, edge);

  const double golden = 0.5 * (1.0 + std::sqrt(5.0));
  const double scale = shellRadius / (std::sqrt(1.0 + golden * golden) * edge);
  for (double first : {-1.0, 1.0})
  {
    for (double second : {-golden, golden})
    {
      lattice.fractionalPositions.push_back(double3(0.5, 0.5 + first * scale, 0.5 + second * scale));
      lattice.fractionalPositions.push_back(double3(0.5 + first * scale, 0.5 + second * scale, 0.5));
      lattice.fractionalPositions.push_back(double3(0.5 + second * scale, 0.5, 0.5 + first * scale));
    }
  }
  lattice.radii.assign(lattice.fractionalPositions.size(), atomRadius);
  return lattice;
}


PoreSizeDistributionCurve curveOf(const Lattice& lattice, double maximumDiameter, std::size_t bins,
                                  double probe = 0.0)
{
  auto build = [&](double probeRadius)
  {
    return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probeRadius);
  };
  return exactPoreSizeDistribution(build, lattice.box.volume, maximumDiameter, bins, 1, probe);
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


// The second curve, over the volume one fixed probe can reach rather than over the whole void. The void of
// this lattice is a single channel system that runs through the cell, so there is nothing for any probe to
// be shut out of and the two curves have to be the same curve, differing only in what they are divided by:
// the accessible one by the volume open at the probe's own diameter, the whole one by the void.
TEST(exact_pore_size_distribution, an_open_void_leaves_the_accessible_curve_the_whole_of_it)
{
  const double probe = 0.8;
  PoreSizeDistributionCurve curve = curveOf(simpleCubic(6.0, 1.5), 9.0, 90, probe);

  EXPECT_EQ(curve.probeRadius, probe);
  ASSERT_GT(curve.probeAccessibleVolume, 0.0);

  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    if (point.diameter < 2.0 * probe)
    {
      // Room for the whole probe is room for its centre and then some, so every point of the region has a
      // pore size of at least the probe's diameter: none of the region has left, and none of it can leave.
      EXPECT_NEAR(point.probeAccessibleCumulative, 1.0, 1.0e-9) << "at d = " << point.diameter;
      EXPECT_NEAR(point.probeAccessibleDistribution, 0.0, 1.0e-9) << "at d = " << point.diameter;
    }
    else
    {
      EXPECT_NEAR(point.probeAccessiblePoreVolume, point.poreVolume, 1.0e-6) << "at d = " << point.diameter;
    }

    // Whatever the framework, the region is inside the void and the one volume bounds the other.
    EXPECT_LE(point.probeAccessiblePoreVolume, point.poreVolume + 1.0e-6) << "at d = " << point.diameter;
  }

  // Both are the same volumes over the same diameters, so the same account has to close on them, and the
  // largest pore is the corner of the cube for either of them.
  EXPECT_NEAR(curve.probeAccessibleIntegral + curve.probeAccessibleSingularWeight, 1.0, 0.02)
      << "continuous " << curve.probeAccessibleIntegral << ", spikes " << curve.probeAccessibleSingularWeight;
  ASSERT_FALSE(curve.probeAccessibleSpikes.empty());
  EXPECT_NEAR(curve.probeAccessibleLargestDiameter, curve.largestDiameter, 1.0e-4);
  EXPECT_LT(curve.probeAccessibleTruncatedWeight, 1.0e-9);
}


// A probe larger than anything the void has room for. The pore volume at its diameter is nothing, so it
// reaches nothing, and the accessible curve is not drawn at all rather than drawn over nothing.
TEST(exact_pore_size_distribution, a_probe_that_does_not_fit_has_no_accessible_curve)
{
  // The largest sphere the void of this lattice holds is sqrt(3) * 6 - 3 = 7.39 across, so a probe of
  // radius 5 has nowhere in it to be.
  PoreSizeDistributionCurve curve = curveOf(simpleCubic(6.0, 1.5), 9.0, 30, 5.0);

  EXPECT_EQ(curve.probeAccessibleVolume, 0.0);
  EXPECT_TRUE(curve.probeAccessibleSpikes.empty());
  EXPECT_NEAR(curve.probeAccessibleIntegral, 0.0, 1.0e-12);
  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    EXPECT_EQ(point.probeAccessiblePoreVolume, 0.0) << "at d = " << point.diameter;
    EXPECT_EQ(point.probeAccessibleCumulative, 0.0) << "at d = " << point.diameter;
    EXPECT_EQ(point.probeAccessibleDistribution, 0.0) << "at d = " << point.diameter;
  }

  // The whole void is unaffected by which probe was named: it is every pore the framework has.
  EXPECT_NEAR(curve.integral + curve.singularWeight, 1.0, 0.05);
}


// A cage the probe fits inside but cannot get inside, which is the case the accessible curve is drawn for. In
// experiment the molecule never enters such a pore, so it is not in the isotherm and it is not to be in this
// curve either, however much room it has: what the curve is over is the channel system, pockets removed.
//
// The cage here holds a sphere 1.8 across, the probe is 1.0 across, and the shell around it is closed. So the
// whole void has a family of pores at 1.8 that the accessible region does not have, and above 1.8 the cage is
// empty and the two agree again. That is a two-sided statement about where the difference is, which a curve
// that merely looked plausible would not satisfy.
TEST(exact_pore_size_distribution, a_cage_the_probe_cannot_enter_is_no_part_of_the_accessible_curve)
{
  const double shellRadius = 2.5;
  const double atomRadius = 1.6;
  const double cageDiameter = 2.0 * (shellRadius - atomRadius);  // the largest sphere the cage holds
  const double probe = 0.5;

  // The range runs past the largest sphere the open part of the cell holds, so that what the curve fails to
  // account for is its own error and not a sweep that stopped early.
  Lattice cage = cageInAnOpenCell(10.0, shellRadius, atomRadius);
  PoreSizeDistributionCurve curve = curveOf(cage, 14.0, 140, probe);

  ASSERT_GT(curve.probeAccessibleVolume, 0.0);
  EXPECT_LT(curve.probeAccessibleVolume, curve.voidVolume);

  double previous = std::numeric_limits<double>::max();
  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    // The volume the probe reaches at one diameter it reaches at every smaller one, so this falls. Nothing
    // else in the curve is as easily broken or as plainly wrong when it is: a cage's own volume put on the
    // wrong side of the wall for one row of the sweep and the right side for the next is a curve that rises.
    EXPECT_LE(point.probeAccessiblePoreVolume, previous + 1.0e-9) << "at d = " << point.diameter;
    previous = point.probeAccessiblePoreVolume;

    if (point.diameter > cageDiameter + 0.02)
    {
      // The cage is empty above the sphere it holds, so there is nothing left for the two to differ by.
      EXPECT_NEAR(point.probeAccessiblePoreVolume, point.poreVolume, 1.0e-6) << "at d = " << point.diameter;
    }
    else if (point.diameter > 2.0 * probe && point.diameter < cageDiameter - 0.02)
    {
      // Below it the cage is a pore of the void and no pore of the region, and the difference is its own
      // volume: little against the cell, but the whole of what the cage has.
      EXPECT_LT(point.probeAccessiblePoreVolume, point.poreVolume - 0.5) << "at d = " << point.diameter;
    }
  }

  // The cage empties all at once, at the size of the sphere it holds, so the whole void has a spike there and
  // the accessible curve has none: that family of pores is not among the ones the probe can reach.
  bool wholeHasIt = false;
  for (const PoreSizeSpike& spike : curve.spikes)
  {
    if (std::abs(spike.diameter - cageDiameter) < 0.05) wholeHasIt = true;
  }
  EXPECT_TRUE(wholeHasIt) << "the void's own spikes are the pores of the framework, the cage among them";
  for (const PoreSizeSpike& spike : curve.probeAccessibleSpikes)
  {
    EXPECT_GT(std::abs(spike.diameter - cageDiameter), 0.05) << "the cage is in the accessible curve after all";
  }

  // And what is left is a curve in its own right: the pores the probe reaches account for the whole of the
  // volume it reaches, by the same comparison of a volume against its own derivative.
  EXPECT_NEAR(curve.probeAccessibleIntegral + curve.probeAccessibleSingularWeight, 1.0, 0.02)
      << "continuous " << curve.probeAccessibleIntegral << ", spikes " << curve.probeAccessibleSingularWeight;
  EXPECT_LT(curve.probeAccessibleTruncatedWeight, 1.0e-9);
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
