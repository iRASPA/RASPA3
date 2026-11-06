#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import pore_accessibility;
import voronoi_network;
import exact_pore_size_distribution;
import exact_void_split;
import voronoi_blocking_spheres;

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
  UnitCell box{6.0, 6.0, 6.0};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
};


Lattice simpleCubic(double edge, double radius)
{
  Lattice lattice;
  lattice.box = UnitCell(edge, edge, edge);
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
  lattice.box = UnitCell(edge, edge, edge);

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
  return exactPoreSizeDistribution(build, lattice.box.volume, maximumDiameter, bins, 1, probe, probe);
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


TEST(exact_pore_size_distribution, peaks_only_finds_the_same_corner_spike)
{
  const double edge = 6.0;
  const double radius = 1.5;
  const double corner = std::sqrt(3.0) * edge - 2.0 * radius;
  const Lattice lattice = simpleCubic(edge, radius);

  auto build = [&](double probeRadius)
  { return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probeRadius); };

  PoreSizeDistributionCurve peaks = exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.0);

  ASSERT_FALSE(peaks.spikes.empty());
  EXPECT_NEAR(peaks.spikes.back().diameter, corner, 1.0e-4) << "bracket " << peaks.spikes.back().bracket;
  EXPECT_GT(peaks.spikes.back().weight, 0.1);
  EXPECT_TRUE(peaks.points.empty());
  EXPECT_LT(peaks.numberOfEvaluations, 40u);
}


// Refining the peaks is asking for the diameter to be measured rather than read off the network, so each spike has
// to come back with the diameter the bisection cornered and not the one it started from, in a bracket as narrow as
// that bisection got. And a room is one spike however the bisection fell: a cliff straddled by two of its intervals
// leaves an excess in both, and those are one room's weight.
TEST(exact_pore_size_distribution, refined_peaks_report_what_the_refinement_found_once_per_room)
{
  const std::vector<std::pair<std::string, Lattice>> cases = {
      {"simple cubic", simpleCubic(6.0, 1.5)},
      {"overlapping simple cubic", simpleCubic(6.0, 3.3)},
      {"cage in an open cell", cageInAnOpenCell(14.0, 3.0, 1.8)},
  };

  for (const auto& [name, lattice] : cases)
  {
    auto build = [&](double probeRadius)
    { return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probeRadius); };

    PoreSizeDistributionCurve coarse = exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.0, 0);
    PoreSizeDistributionCurve refined = exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.0, 12);

    ASSERT_FALSE(coarse.spikes.empty()) << name;
    ASSERT_EQ(refined.spikes.size(), coarse.spikes.size()) << name;
    ASSERT_EQ(refined.probeAccessibleSpikes.size(), coarse.probeAccessibleSpikes.size()) << name;

    for (std::size_t s = 0; s < refined.spikes.size(); ++s)
    {
      const PoreSizeSpike& before = coarse.spikes[s];
      const PoreSizeSpike& after = refined.spikes[s];

      EXPECT_LT(after.bracket, 1.0e-2 * before.bracket) << name << ", spike " << s;
      EXPECT_NE(after.diameter, before.diameter) << name << ", spike " << s;
      EXPECT_NEAR(after.diameter, before.diameter, before.bracket) << name << ", spike " << s;
      EXPECT_NEAR(after.weight, before.weight, 1.0e-3) << name << ", spike " << s;
      if (s > 0)
      {
        EXPECT_GT(after.diameter, refined.spikes[s - 1].diameter) << name << ", spike " << s;
      }
    }
    EXPECT_NEAR(refined.singularWeight, coarse.singularWeight, 1.0e-3) << name;
    EXPECT_NEAR(refined.probeAccessibleSingularWeight, coarse.probeAccessibleSingularWeight, 1.0e-3) << name;
  }
}


// Those lattices are symmetric enough that the bisection tends to land clear of the cliff. A disordered
// arrangement does not oblige: on these two, the cliffs of several rooms are straddled at twelve refinements and
// would each come back as a pair of spikes. Every room's weight has to come back together, as the one spike the
// unrefined bracket finds there.
TEST(exact_pore_size_distribution, a_room_whose_cliff_the_bisection_straddles_is_one_spike)
{
  for (unsigned seed : {35u, 38u})
  {
    std::mt19937 engine(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    Lattice lattice;
    lattice.box = UnitCell(10.0, 11.0, 12.0);
    const std::size_t atoms = 6 + seed % 20;
    for (std::size_t atom = 0; atom < atoms; ++atom)
    {
      lattice.fractionalPositions.push_back(double3(unit(engine), unit(engine), unit(engine)));
      lattice.radii.push_back(1.2 + unit(engine));
    }
    auto build = [&](double probeRadius)
    { return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probeRadius); };

    PoreSizeDistributionCurve coarse = exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.0, 0);
    PoreSizeDistributionCurve refined = exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.0, 12);
    ASSERT_FALSE(refined.spikes.empty()) << "seed " << seed;

    for (std::size_t s = 0; s < refined.spikes.size(); ++s)
    {
      const PoreSizeSpike& spike = refined.spikes[s];

      // The candidates are at least 1e-2 Å apart and each bracket reaches a quarter of the way to the next.
      if (s > 0) EXPECT_GT(spike.diameter - refined.spikes[s - 1].diameter, 5.0e-3) << "seed " << seed;

      auto nearest = std::ranges::min_element(coarse.spikes, {},
                                              [&](const PoreSizeSpike& other)
                                              { return std::abs(other.diameter - spike.diameter); });
      ASSERT_NE(nearest, coarse.spikes.end());
      EXPECT_NEAR(nearest->diameter, spike.diameter, 1.0e-2) << "seed " << seed << ", spike " << s;
      EXPECT_NEAR(nearest->weight, spike.weight, 1.0e-3) << "seed " << seed << ", spike " << s;
    }
  }
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


namespace
{

// A sealed cage of `count` spheres on a shell about `centre`, placed on the Fibonacci spiral.
std::vector<double3> fibonacciShell(const double3& centre, double shellRadius, std::size_t count)
{
  std::vector<double3> positions;
  const double golden = std::numbers::pi * (1.0 + std::sqrt(5.0));
  for (std::size_t i = 0; i < count; ++i)
  {
    double fraction = (static_cast<double>(i) + 0.5) / static_cast<double>(count);
    double cosinePolar = 1.0 - 2.0 * fraction;
    double sinePolar = std::sqrt(std::max(0.0, 1.0 - cosinePolar * cosinePolar));
    double azimuth = golden * (static_cast<double>(i) + 0.5);
    positions.push_back(centre + double3(shellRadius * sinePolar * std::cos(azimuth),
                                         shellRadius * sinePolar * std::sin(azimuth), shellRadius * cosinePolar));
  }
  return positions;
}


Lattice latticeOf(double edge, const std::vector<double3>& cartesian, double atomRadius)
{
  Lattice lattice;
  lattice.box = UnitCell(edge, edge, edge);
  for (const double3& position : cartesian)
  {
    lattice.fractionalPositions.push_back(double3::fract(lattice.box.inverseCell * position));
  }
  lattice.radii.assign(cartesian.size(), atomRadius);
  return lattice;
}


// The hybrid route as the drivers take it: the pockets sealed to `probe` are measured on that probe's network,
// their blocking spheres are put among the bare atoms, and the bare curve of that blocked structure is swept.
struct Hybrid
{
  ExactVoidSplit split;
  PoreSizeDistributionCurve blocked;
};

Hybrid hybridOf(const Lattice& lattice, double probe, double maximumDiameter, std::size_t bins)
{
  PoreAccessibility he = PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probe);
  ExactVoidSplit split = exactVoidSplitByComponents(he, lattice.box.volume, 1);
  std::vector<BlockingSphere> balls =
      hybridBlockedBalls(he, split, lattice.fractionalPositions, lattice.radii, probe, 1000);

  std::vector<double3> positions;
  std::vector<double> radii;
  for (const BlockingSphere& ball : balls)
  {
    positions.push_back(ball.centerFractional);
    radii.push_back(ball.radius);
  }
  auto build = [&](double inflation) { return PoreAccessibility::create(lattice.box, positions, radii, inflation); };
  return Hybrid{split, exactPoreSizeDistribution(build, lattice.box.volume, maximumDiameter, bins, 1, 0.0, 0.0)};
}

}  // namespace


// A blocking sphere is drawn where the probe's centre can be, and the hybrid puts it among the bare atoms. The
// cavity of a sealed cage reaches a probe radius beyond the wall the probe's centre sees, so a sphere of the
// covering radius alone leaves a shell of the cavity open, and at vanishing diameter that shell is counted as
// void of the open network: the normalisation of the whole curve is then wrong by it, and the shell shows up as
// pores of its own thickness. Carried into the bare structure the sphere grows by the probe radius, and the
// blocked void is the bare void less the whole cavity.
TEST(exact_pore_size_distribution, the_hybrid_blocks_the_whole_cavity_of_a_sealed_cage)
{
  const double edge = 14.0;
  const double probe = 0.5;
  Lattice cage = latticeOf(edge, fibonacciShell(double3(7.0, 7.0, 7.0), 3.0, 32), 1.8);

  PoreAccessibility bare = PoreAccessibility::create(cage.box, cage.fractionalPositions, cage.radii, 0.0);
  ExactVoidSplit bareSplit = exactVoidSplitByComponents(bare, cage.box.volume, 1);
  ASSERT_TRUE(bareSplit.reliable) << bareSplit.rejection;
  ASSERT_EQ(bareSplit.pockets.size(), 1uz);
  const double cavity = bareSplit.inaccessibleVolume;

  Hybrid hybrid = hybridOf(cage, probe, 4.0, 40);
  ASSERT_TRUE(hybrid.split.reliable) << hybrid.split.rejection;
  ASSERT_EQ(hybrid.split.pockets.size(), 1uz);
  ASSERT_TRUE(hybrid.split.pockets.front().coversPocket());

  // The ball of the covering radius plus the probe radius holds every point of the cavity a helium atom could
  // occupy; what it can miss is the crevice between three cage atoms that not even helium reaches, so the
  // blocked void may exceed the bare void less the cavity by a little, and never fall short of it.
  const double expected = bareSplit.voidVolume - cavity;
  EXPECT_GT(hybrid.blocked.voidVolume, expected - 1.0e-6);
  EXPECT_LT(hybrid.blocked.voidVolume, expected + 0.05 * cavity);
}


// The other limit goes the other way. Where the channel is nearer the pocket's centre than its own far wall,
// the sphere is capped by the channel radius, measured to where the probe's centre can get; a molecule's
// surface comes a probe radius nearer than that, so among the bare atoms the cap has to come down by one.
// Left as it was, the sphere inflated by each row's radius reaches a probe radius into the channel, and the
// blocked curve loses channel volume that the division of the same void by its surfaces keeps. Above the
// probe's diameter the two routes measure the same thing and have to agree.
TEST(exact_pore_size_distribution, the_hybrid_leaves_the_channel_round_a_capped_pocket_alone)
{
  const double edge = 24.0;
  const double probe = 0.5;
  const double shellRadius = 3.0;
  const double atomRadius = 1.0;

  // Three cages in a row, each shell cut where it enters a neighbour, so that the three cavities are one
  // elongated pocket: its far walls lie farther from its centre than the channel outside its middle does.
  std::vector<double3> centres = {double3(8.0, 12.0, 12.0), double3(12.0, 12.0, 12.0), double3(16.0, 12.0, 12.0)};
  std::vector<double3> positions;
  for (std::size_t cage = 0; cage < centres.size(); ++cage)
  {
    for (const double3& position : fibonacciShell(centres[cage], shellRadius, 64))
    {
      bool insideAnother = false;
      for (std::size_t other = 0; other < centres.size(); ++other)
      {
        if (other != cage && (position - centres[other]).length() < shellRadius - 1.0e-9) insideAnother = true;
      }
      if (!insideAnother) positions.push_back(position);
    }
  }
  Lattice chain = latticeOf(edge, positions, atomRadius);

  Hybrid hybrid = hybridOf(chain, probe, 4.4, 44);
  ASSERT_TRUE(hybrid.split.reliable) << hybrid.split.rejection;
  ASSERT_EQ(hybrid.split.pockets.size(), 1uz);
  const PocketGeometry& pocket = hybrid.split.pockets.front();
  ASSERT_LT(pocket.channelRadius, pocket.coveringRadius) << "the test wants a pocket its channel caps";

  // Above the diameter at which the capped sphere holds what is left of the pocket, the blocked curve and the
  // division by surfaces agree: the pocket is gone from both, and the channel is whole in both.
  PoreSizeDistributionCurve divided = curveOf(chain, 4.4, 44, probe);
  // The sphere of the channel radius less the probe radius, inflated by r, holds what is left of the pocket at r
  // once r reaches half the excess of the covering radius over the channel radius, plus the probe radius.
  const double covered = pocket.coveringRadius - pocket.channelRadius + 2.0 * probe;
  std::size_t compared = 0;
  for (std::size_t row = 0; row < hybrid.blocked.points.size(); ++row)
  {
    const PoreSizeDistributionPoint& blocked = hybrid.blocked.points[row];
    const PoreSizeDistributionPoint& reference = divided.points[row];
    ASSERT_NEAR(blocked.diameter, reference.diameter, 1.0e-12);
    if (blocked.diameter < covered + 0.1 || blocked.diameter > 2.0 * (shellRadius - atomRadius) - 0.1) continue;
    ++compared;
    EXPECT_NEAR(blocked.poreVolume, reference.probeAccessiblePoreVolume, 1.0e-4 * reference.probeAccessiblePoreVolume)
        << "at d = " << blocked.diameter;
  }
  EXPECT_GT(compared, 3uz) << "covering " << pocket.coveringRadius << ", channel " << pocket.channelRadius;
}


// Below the probe's diameter the surfaces cannot divide the void: a pocket whose bare windows are open is one
// surface with the channel outside it until the inflated atoms close them, so its volume would be kept until
// then and dropped after. The drivers take the hybrid through a blocked structure instead; asking this route
// for a floor below a positive probe is refused rather than answered wrongly.
TEST(exact_pore_size_distribution, a_floor_below_the_probe_is_refused)
{
  Lattice lattice = simpleCubic(6.0, 1.5);
  auto build = [&](double inflation)
  { return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, inflation); };
  EXPECT_THROW(exactPoreSizeDistribution(build, lattice.box.volume, 6.0, 12, 1, 0.5, 0.0), std::invalid_argument);
  EXPECT_THROW(exactPoreSizeDistribution(build, lattice.box.volume, 6.0, 12, 1, 0.5, 0.25), std::invalid_argument);
  EXPECT_THROW(exactPoreSizePeaks(build, lattice.box.volume, 1, 0.5, 0.0), std::invalid_argument);
  EXPECT_NO_THROW(exactPoreSizeDistribution(build, lattice.box.volume, 6.0, 12, 1, 0.5, 0.5));
  EXPECT_NO_THROW(exactPoreSizeDistribution(build, lattice.box.volume, 6.0, 12, 1, 0.0, 0.0));
}


// The sweep stops one step past the largest sphere the network reports, which is an estimate: its nodes are
// walked up to the peaks of the clearance, and a walk that fell short would leave the estimate below the
// sphere that is really there. Stopped there, the rows beyond would be zeros that are not true, a cliff where
// the void is still whole, and the search would corner a spike at that row and never see the real one. The
// estimate is handed in short on purpose here, by two Ångström, and the curve has to come out as if it had
// not been: the same corner spike, nothing truncated, the void accounted for.
TEST(exact_pore_size_distribution, a_short_estimate_of_the_largest_sphere_does_not_end_the_curve)
{
  const double edge = 6.0;
  const double radius = 1.5;
  const double corner = std::sqrt(3.0) * edge - 2.0 * radius;
  Lattice lattice = simpleCubic(edge, radius);

  auto honest = [&](double inflation)
  { return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, inflation); };
  auto understated = [&](double inflation)
  {
    PoreAccessibility accessibility = honest(inflation);
    if (inflation == 0.0)
    {
      // The estimate is read off the bare network's nodes, and only there.
      for (VoronoiNode& node : accessibility.network.nodes)
      {
        node.maximalRadius = std::min(node.maximalRadius, 0.5 * (corner - 2.0));
      }
    }
    return accessibility;
  };

  PoreSizeDistributionCurve reference = exactPoreSizeDistribution(honest, lattice.box.volume, 9.0, 90, 1);
  PoreSizeDistributionCurve curve = exactPoreSizeDistribution(understated, lattice.box.volume, 9.0, 90, 1);

  ASSERT_FALSE(curve.spikes.empty());
  EXPECT_NEAR(curve.spikes.back().diameter, corner, 1.0e-4) << "the curve ends where the void does";
  EXPECT_NEAR(curve.largestDiameter, corner, 1.0e-4);
  EXPECT_LT(curve.truncatedWeight, 1.0e-9);
  EXPECT_NEAR(curve.integral + curve.singularWeight, 1.0, 0.02);
  for (const PoreSizeSpike& spike : curve.spikes)
  {
    EXPECT_GT(std::abs(spike.diameter - (corner - 2.0)), 0.5) << "no spike at the short estimate";
  }

  // And it is the same curve, row for row, as the one taken from the sound estimate.
  ASSERT_EQ(curve.spikes.size(), reference.spikes.size());
  EXPECT_NEAR(curve.spikes.back().weight, reference.spikes.back().weight, 1.0e-9);
  ASSERT_EQ(curve.points.size(), reference.points.size());
  for (std::size_t row = 0; row < curve.points.size(); ++row)
  {
    EXPECT_NEAR(curve.points[row].poreVolume, reference.points[row].poreVolume, 1.0e-9) << "row " << row;
  }
}


// A room only just wider than the floor. The peaks route samples once below and once above each room, five
// thousandths of an Ångström either way, and the primary curve is collected only from the floor diameter up: so
// a room within that distance above the floor had the lower of its two samples under the floor, and its spike
// was left out of the primary curve altogether, though the room is in it and here is nearly all of it. The floor
// is the point the primary curve turns at, and what a bracket reaching across it holds of that curve is taken
// from the floor.
TEST(exact_pore_size_distribution, peaks_only_keeps_a_room_that_is_only_just_wider_than_the_floor)
{
  const double edge = 6.0;
  const double radius = 1.5;
  const double corner = std::sqrt(3.0) * edge - 2.0 * radius;
  const Lattice lattice = simpleCubic(edge, radius);

  auto build = [&](double probeRadius)
  { return PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probeRadius); };

  // The whole void does not know of the floor, so its one spike is the same whatever the floor is.
  const PoreSizeDistributionCurve unfloored = exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.0);
  ASSERT_EQ(unfloored.spikes.size(), 1uz);

  // From a room well clear of the floor, where nothing was ever wrong, down to one half a thousandth above it.
  for (const double above : {0.02, 0.004, 0.002, 0.0005})
  {
    std::array<double, 2> weights{};
    for (const std::size_t refinements : {0uz, 12uz})
    {
      const std::string which = std::format("room {} Å above the floor, {} refinements", above, refinements);
      const PoreSizeDistributionCurve peaks =
          exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.5 * (corner - above), refinements);

      ASSERT_EQ(peaks.spikes.size(), 1uz) << which;
      EXPECT_NEAR(peaks.spikes.back().weight, unfloored.spikes.back().weight, 1.0e-6) << which;

      // What is left at the floor is the corner room and next to nothing else, so its spike is all but the
      // whole of the primary curve.
      EXPECT_GT(peaks.probeAccessibleVolume, 0.0) << which;
      ASSERT_EQ(peaks.probeAccessibleSpikes.size(), 1uz) << which;
      const PoreSizeSpike& spike = peaks.probeAccessibleSpikes.back();
      EXPECT_NEAR(spike.diameter, corner, 1.0e-3) << which;
      EXPECT_GT(spike.weight, 0.998) << which;
      EXPECT_LE(spike.weight, 1.0 + 1.0e-9) << which;
      EXPECT_DOUBLE_EQ(peaks.probeAccessibleSingularWeight, spike.weight) << which;
      EXPECT_NEAR(peaks.probeAccessibleLargestDiameter, corner, 1.0e-3) << which;
      weights[refinements == 0 ? 0 : 1] = spike.weight;
    }
    // Refining the bracket corners the same cliff, so the weight does not depend on it.
    EXPECT_NEAR(weights[0], weights[1], 5.0e-4) << "room " << above << " Å above the floor";
  }

  // The floor above the room instead: nothing is left at the floor, and nothing of the primary curve is the
  // room's, while the whole void still has it.
  const PoreSizeDistributionCurve below =
      exactPoreSizePeaks(build, lattice.box.volume, 1, 0.0, 0.5 * (corner + 0.002));
  EXPECT_EQ(below.probeAccessibleVolume, 0.0);
  EXPECT_TRUE(below.probeAccessibleSpikes.empty());
  ASSERT_EQ(below.spikes.size(), 1uz);
  EXPECT_NEAR(below.spikes.back().weight, unfloored.spikes.back().weight, 1.0e-6);
}
