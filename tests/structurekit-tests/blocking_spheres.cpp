#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import pore_accessibility;
import voronoi_blocking_spheres;

// The sampled blocking spheres, and the allowances made for the sampling.
//
// A sphere may reach no nearer the channel than the nearest accessible sample, less how far that sample can
// stand off the channel's edge, which is the spacing of the accessible samples. The samples are thrown over the
// whole cell, so that spacing is set by the cell's volume and the number thrown, and not by how much of the
// cell turned out to be accessible.

namespace
{

struct Lattice
{
  UnitCell box{20.0, 20.0, 20.0};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;

  void add(const double3& cartesian, double radius)
  {
    this->fractionalPositions.push_back(double3::fract(this->box.inverseCell * cartesian));
    this->radii.push_back(radius);
  }
};

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

const double3 cageCentre(10.0, 10.0, 10.0);
constexpr double atomRadius = 1.8;
constexpr double probe = 0.5;

// A cage in the middle of the cell, with the void round it as the channel. With `filledBeyond` set, the
// corners of the cell farther than that from the cage's centre are packed solid, which takes away a good
// part of the channel but none of it anywhere near the cage.
Lattice cageInACell(std::optional<double> filledBeyond)
{
  Lattice lattice;
  for (const double3& position : fibonacciShell(cageCentre, 6.0, 128)) lattice.add(position, atomRadius);
  if (filledBeyond.has_value())
  {
    for (int i = 0; i < 10; ++i)
      for (int j = 0; j < 10; ++j)
        for (int k = 0; k < 10; ++k)
        {
          double3 position(2.0 * i, 2.0 * j, 2.0 * k);
          double3 delta = lattice.box.applyPeriodicBoundaryConditions(position - cageCentre);
          if (delta.length() >= filledBeyond.value()) lattice.add(position, atomRadius);
        }
  }
  return lattice;
}

std::vector<BlockingSphere> spheresOf(const Lattice& lattice, std::size_t samples)
{
  PoreAccessibility accessibility =
      PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probe);
  return computeBlockingSpheres(accessibility, samples);
}

}  // namespace


TEST(blocking_spheres, the_channel_margin_is_the_spacing_of_the_samples_and_not_of_the_accessible_ones)
{
  // Few enough samples that the margin, and not the far wall of the cavity, is what bounds the sphere.
  const std::size_t samples = 300;

  std::vector<BlockingSphere> open = spheresOf(cageInACell(std::nullopt), samples);
  std::vector<BlockingSphere> cramped = spheresOf(cageInACell(12.0), samples);

  ASSERT_FALSE(open.empty());
  ASSERT_EQ(cramped.size(), open.size());

  // The same samples land in the cavity and the same ones land just outside it, so the spheres are the same.
  // Taking the spacing from the count of accessible samples alone instead would have widened the margin by the
  // cube root of the part of the cell the corners took away, and held the sphere short by that.
  for (std::size_t i = 0; i < open.size(); ++i)
  {
    EXPECT_NEAR(cramped[i].radius, open[i].radius, 1.0e-6) << "sphere " << i;
    EXPECT_NEAR((cramped[i].centerFractional - open[i].centerFractional).length(), 0.0, 1.0e-6) << "sphere " << i;
  }

  // And the sphere is bounded by the channel, not by the free ball at its centre or by its own far wall: it
  // reaches past the cavity's room, and stops short of the cage's outer surface.
  const double room = 6.0 - atomRadius - probe;
  EXPECT_GT(open.front().radius, room + 1.0);
  EXPECT_LT(open.front().radius, 6.0 + atomRadius + probe);
}


// The same spheres wanted among the bare atoms, where a pocket's cavity goes on a probe radius past the
// outermost position of the probe's centre. A sampled position then stands for the ball of that radius about
// it, and is dealt with when that ball is inside a sphere, not when the position is.

namespace
{

constexpr double thinAtomRadius = 1.2;
constexpr double wideProbe = 1.0;

// A cage the shape of a capsule along x: a tube of rings closed by two half shells. A sphere at its middle
// that reaches the ends would reach through the sides, so the channel outside holds it short and the ends are
// left for further spheres.
Lattice capsuleCage()
{
  Lattice lattice;
  lattice.box = UnitCell(26.0, 16.0, 16.0);
  const double3 middle(13.0, 8.0, 8.0);
  const double tube = 4.0;
  const double half = 5.0;
  for (int ring = -2; ring <= 2; ++ring)
  {
    for (int k = 0; k < 13; ++k)
    {
      const double angle = 2.0 * std::numbers::pi * (k + 0.5 * (ring & 1)) / 13.0;
      lattice.add(middle + double3(2.0 * ring, tube * std::cos(angle), tube * std::sin(angle)), thinAtomRadius);
    }
  }
  for (const double3& position : fibonacciShell(double3(0.0, 0.0, 0.0), tube, 60))
  {
    // The shell is wound about z, and the half of it on either side of x = 0 closes that end of the tube.
    lattice.add(middle + position + double3(position.x >= 0.0 ? half : -half, 0.0, 0.0), thinAtomRadius);
  }
  return lattice;
}

// A cage of few enough atoms that its wall is open between them and shut only to the probe: the channel
// outside comes within a probe's width of the cavity at every gap.
Lattice thinWalledCage()
{
  Lattice lattice;
  lattice.box = UnitCell(22.0, 22.0, 22.0);
  for (const double3& position : fibonacciShell(double3(11.0, 11.0, 11.0), 6.0, 54))
  {
    lattice.add(position, thinAtomRadius);
  }
  return lattice;
}

struct BareCover
{
  std::vector<BlockingSphere> spheres;
  double resolution{0.0};
  std::size_t pocketPositions{0};
  std::size_t roomyPositions{0};      // of those, the ones with more room than the samples are apart
  std::size_t standingOut{0};         // positions whose ball of the probe's radius is inside no one sphere
  std::size_t roomyStandingOut{0};
  double nearestCentres{std::numeric_limits<double>::max()};
  double nearestToAChannelPosition{std::numeric_limits<double>::max()};  // surface of a sphere to a position
};

BareCover bareCoverOf(const Lattice& lattice, std::size_t samples)
{
  PoreAccessibility accessibility =
      PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, wideProbe);

  BareCover cover;
  cover.spheres = computeBlockingSpheres(accessibility, samples, wideProbe);
  cover.resolution = std::cbrt(lattice.box.volume / static_cast<double>(samples));

  std::vector<double3> centres;
  for (const BlockingSphere& sphere : cover.spheres) centres.push_back(lattice.box.cell * sphere.centerFractional);
  auto apart = [&](const double3& a, const double3& b)
  { return lattice.box.applyPeriodicBoundaryConditions(a - b).length(); };

  for (std::size_t i = 0; i < centres.size(); ++i)
  {
    for (std::size_t j = i + 1; j < centres.size(); ++j)
    {
      cover.nearestCentres = std::min(cover.nearestCentres, apart(centres[i], centres[j]));
    }
  }

  RandomNumber random{std::optional<std::size_t>(99)};
  for (std::size_t s = 0; s < 200000; ++s)
  {
    const double3 position = lattice.box.cell * double3(random.uniform(), random.uniform(), random.uniform());
    const PointClassification classification = accessibility.classify(position);
    if (classification.inside || classification.resample) continue;

    bool held = false;
    for (std::size_t k = 0; k < centres.size(); ++k)
    {
      const double distance = apart(position, centres[k]);
      if (distance + wideProbe < cover.spheres[k].radius) held = true;
      if (classification.accessible)
      {
        cover.nearestToAChannelPosition =
            std::min(cover.nearestToAChannelPosition, distance - cover.spheres[k].radius);
      }
    }
    if (classification.accessible) continue;

    const bool roomy = accessibility.clearance(position) > cover.resolution;
    ++cover.pocketPositions;
    if (roomy) ++cover.roomyPositions;
    if (!held)
    {
      ++cover.standingOut;
      if (roomy) ++cover.roomyStandingOut;
    }
  }
  return cover;
}

}  // namespace

TEST(blocking_spheres, among_the_bare_atoms_a_position_is_covered_when_its_whole_ball_is)
{
  // The sphere at the middle of the capsule is held short by the channel beside it. The positions in the last
  // probe radius of it have their balls standing out of it, and used to be struck off as covered all the same,
  // with nothing centred beyond to take them in: one position in twelve of the pocket at the coarser sampling
  // here, and one in sixty at the finer.
  const BareCover coarse = bareCoverOf(capsuleCage(), 20000);
  ASSERT_GT(coarse.pocketPositions, 3000uz);
  EXPECT_GT(coarse.spheres.size(), 2uz);
  EXPECT_LT(100 * coarse.standingOut, coarse.pocketPositions);
  EXPECT_EQ(coarse.roomyStandingOut, 0uz);

  const BareCover fine = bareCoverOf(capsuleCage(), 300000);
  EXPECT_LT(500 * fine.standingOut, fine.pocketPositions);
  EXPECT_EQ(fine.roomyStandingOut, 0uz);

  // And no sphere comes within a probe radius of a position the probe can take from the channel, which is
  // where the bare cavity of the channel begins.
  EXPECT_GE(coarse.nearestToAChannelPosition, wideProbe);
  EXPECT_GE(fine.nearestToAChannelPosition, wideProbe);
}

TEST(blocking_spheres, among_the_bare_atoms_a_thin_wall_does_not_draw_a_sphere_at_every_sample)
{
  // With the channel a probe's width away at every gap in the wall, the sphere over the middle is held well
  // short and the rest is left to spheres of little more than the room at their own centres. Those are never
  // thinner than the probe, so a test of their radius against the spacing of the samples lets every one of
  // them through, one at each sample along the rim and as close together as the samples are. Held to how far
  // they reach among the positions instead, the ones that reach less than a sample's spacing are not written,
  // and no two centres are nearer than that.
  for (std::size_t samples : {20000uz, 300000uz})
  {
    const BareCover cover = bareCoverOf(thinWalledCage(), samples);
    ASSERT_GT(cover.pocketPositions, 3000uz) << samples << " samples";
    EXPECT_GT(cover.spheres.size(), 1uz) << samples << " samples";
    EXPECT_GE(cover.nearestCentres, cover.resolution) << samples << " samples";
    EXPECT_EQ(cover.roomyStandingOut, 0uz) << samples << " samples";
  }
}
