#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import sampled_structure;
import voronoi_network;
import voronoi_pore_diameters;
import pore_accessibility;

// A cell may be handed in by its three vectors in any order, and in one order out of two they are a left-handed
// set, with a determinant below zero. It is the same cell, holding the same atoms in the same places, and a
// volume is not less than nothing because of the order its edges were named in.
//
// The volume used to be the determinant as it came. The widths across the cell were then worked out negative
// wherever they were taken from it, and a search that walks outward shell by shell until the shells are too far
// off to matter, the distance to a shell being its number times such a width, never got any further off.

namespace
{

struct TwoHands
{
  UnitCell right;
  UnitCell left;                      // the same vectors, the first two changed over
  std::vector<double3> positions;     // Cartesian, the same in both
  std::vector<double3> inRight;       // and in fractions of each
  std::vector<double3> inLeft;
  std::vector<double> radii;
};

TwoHands scattered(std::size_t atoms, std::size_t seed)
{
  const double3 a(11.0, 0.0, 0.0);
  const double3 b(2.0, 10.0, 0.0);
  const double3 c(-1.0, 3.0, 9.0);

  TwoHands both{.right = UnitCell(double3x3(a, b, c)), .left = UnitCell(double3x3(b, a, c))};
  RandomNumber random{std::optional<std::size_t>(seed)};
  for (std::size_t i = 0; i < atoms; ++i)
  {
    const double3 fractional(random.uniform(), random.uniform(), random.uniform());
    both.inRight.push_back(fractional);
    both.inLeft.push_back(double3(fractional.y, fractional.x, fractional.z));
    both.positions.push_back(both.right.cell * fractional);
    both.radii.push_back(0.8 + 0.6 * random.uniform());
  }
  return both;
}

SampledStructure sampled(const UnitCell &cell, const TwoHands &both)
{
  SampledStructure structure;
  structure.name = "sampled_structure_test";
  structure.unitCell = cell;
  structure.positions = both.positions;
  structure.radii = both.radii;
  return structure;
}

}  // namespace

TEST(sampled_structure, a_cell_is_as_large_whichever_hand_its_vectors_are)
{
  const TwoHands both = scattered(1, 1);
  EXPECT_GT(both.right.volume, 0.0);
  EXPECT_LT(both.left.cell.determinant(), 0.0);
  EXPECT_NEAR(both.left.volume, both.right.volume, 1.0e-9);

  const double3 rightWidths = both.right.perpendicularWidths();
  const double3 leftWidths = both.left.perpendicularWidths();
  EXPECT_NEAR(leftWidths.x, rightWidths.y, 1.0e-12);
  EXPECT_NEAR(leftWidths.y, rightWidths.x, 1.0e-12);
  EXPECT_NEAR(leftWidths.z, rightWidths.z, 1.0e-12);
}

TEST(sampled_structure, the_room_in_a_left_handed_cell_is_the_room_in_the_cell)
{
  const TwoHands both = scattered(24, 2);
  const SampledStructure right = sampled(both.right, both);
  const SampledStructure left = sampled(both.left, both);

  EXPECT_NEAR(left.density(), right.density(), 1.0e-9 * std::abs(right.density()));

  RandomNumber random{std::optional<std::size_t>(5)};
  for (std::size_t s = 0; s < 2000; ++s)
  {
    const double3 point = both.right.cell * double3(3.0 * random.uniform() - 1.0, 3.0 * random.uniform() - 1.0,
                                                    3.0 * random.uniform() - 1.0);
    const double3 step(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);

    EXPECT_NEAR(left.clearance(point), right.clearance(point), 1.0e-9) << "point " << s;
    EXPECT_EQ(left.overlaps(point, both.positions.size()), right.overlaps(point, both.positions.size()))
        << "point " << s;

    const SegmentBottleneck narrowest = left.segmentBottleneck(point, 2.0 * step);
    const SegmentBottleneck expected = right.segmentBottleneck(point, 2.0 * step);
    EXPECT_NEAR(narrowest.radius, expected.radius, 1.0e-9) << "point " << s;
  }
}

TEST(sampled_structure, the_network_of_a_left_handed_cell_is_the_network_of_the_cell)
{
  const TwoHands both = scattered(24, 2);

  const VoronoiNetwork right = VoronoiNetwork::create(both.right, both.inRight, both.radii);
  const VoronoiNetwork left = VoronoiNetwork::create(both.left, both.inLeft, both.radii);
  EXPECT_EQ(left.nodes.size(), right.nodes.size());
  EXPECT_EQ(left.edges.size(), right.edges.size());

  const PoreDiameters mine = PoreDiameters::compute(left);
  const PoreDiameters theirs = PoreDiameters::compute(right);
  EXPECT_NEAR(mine.includedSphereDiameter, theirs.includedSphereDiameter, 1.0e-9);
  EXPECT_NEAR(mine.freeSphereDiameter, theirs.freeSphereDiameter, 1.0e-9);
  EXPECT_NEAR(mine.includedAlongFreePathDiameter, theirs.includedAlongFreePathDiameter, 1.0e-9);

  const PoreAccessibility rightHanded = PoreAccessibility::create(both.right, both.inRight, both.radii, 0.3);
  const PoreAccessibility leftHanded = PoreAccessibility::create(both.left, both.inLeft, both.radii, 0.3);
  RandomNumber random{std::optional<std::size_t>(9)};
  for (std::size_t s = 0; s < 2000; ++s)
  {
    const double3 point = both.right.cell * double3(random.uniform(), random.uniform(), random.uniform());
    EXPECT_NEAR(leftHanded.clearance(point), rightHanded.clearance(point), 1.0e-9) << "point " << s;

    const PointClassification here = leftHanded.classify(point);
    const PointClassification there = rightHanded.classify(point);
    EXPECT_EQ(here.inside, there.inside) << "point " << s;
    EXPECT_EQ(here.accessible, there.accessible) << "point " << s;
  }
}


// The bins a structure's queries look through are built the first time one of them is asked, and built again
// when the structure is seen to have changed since. Many threads may ask at once, and all of them have to get
// the bins of the structure as it now is: none may take bins that another thread is in the middle of
// rebuilding for bins that are ready.
//
// A structure is changed here between rounds, by one thread while no other is looking, and then asked about by
// several at the same moment. What each is told is held against the answer worked out by looking at every atom
// and every copy of it near enough to matter.
TEST(sampled_structure, queries_from_many_threads_after_a_change_all_see_the_structure_as_it_is)
{
  const UnitCell cell(9.0, 10.0, 11.0);
  RandomNumber random{std::optional<std::size_t>(21)};

  SampledStructure structure;
  structure.name = "sampled_structure_test";
  structure.unitCell = cell;
  for (std::size_t i = 0; i < 40; ++i)
  {
    structure.positions.push_back(cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
    structure.radii.push_back(1.0);
  }

  auto roomByLookingAtEveryAtom = [&](const double3 &point)
  {
    double room = std::numeric_limits<double>::max();
    for (std::size_t j = 0; j < structure.positions.size(); ++j)
    {
      for (int x = -1; x <= 1; ++x)
      {
        for (int y = -1; y <= 1; ++y)
        {
          for (int z = -1; z <= 1; ++z)
          {
            const double3 copy = structure.positions[j] + cell.cell * double3(x, y, z);
            room = std::min(room, (point - copy).length() - structure.radii[j]);
          }
        }
      }
    }
    return room;
  };

  constexpr std::size_t numberOfThreads = 8;
  constexpr std::size_t numberOfRounds = 1500;
  constexpr std::size_t pointsPerThread = 3;

  std::vector<double3> points(numberOfThreads * pointsPerThread);
  std::vector<double> expected(points.size());
  std::vector<double> answered(points.size());
  std::atomic<std::size_t> round{0};
  std::atomic<std::size_t> finished{0};
  std::atomic<bool> done{false};

  std::vector<std::thread> threads;
  for (std::size_t t = 0; t < numberOfThreads; ++t)
  {
    threads.emplace_back(
        [&, t]()
        {
          std::size_t seen = 0;
          for (;;)
          {
            // Spinning and not sleeping, so that all of them set off within the same few instructions.
            std::size_t now = round.load(std::memory_order_acquire);
            while (now == seen && !done.load(std::memory_order_acquire)) now = round.load(std::memory_order_acquire);
            if (now == seen) return;
            seen = now;
            for (std::size_t p = 0; p < pointsPerThread; ++p)
            {
              answered[t * pointsPerThread + p] = structure.clearance(points[t * pointsPerThread + p]);
            }
            finished.fetch_add(1, std::memory_order_release);
          }
        });
  }

  std::size_t wrong = 0;
  for (std::size_t r = 0; r < numberOfRounds; ++r)
  {
    // Changed while nobody is looking: an atom more or an atom fewer, somewhere else, and the largest radius
    // moved, so that the bins of the round before are of no use for this one.
    if (r % 2 == 0)
    {
      structure.positions.push_back(cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
      structure.radii.push_back(1.0);
    }
    else
    {
      structure.positions.erase(structure.positions.begin());
      structure.radii.erase(structure.radii.begin());
    }
    for (double3 &position : structure.positions)
    {
      position = cell.cell * double3(random.uniform(), random.uniform(), random.uniform());
    }
    structure.radii.back() = 1.0 + 0.5 * random.uniform();

    for (std::size_t p = 0; p < points.size(); ++p)
    {
      points[p] = cell.cell * double3(random.uniform(), random.uniform(), random.uniform());
      expected[p] = roomByLookingAtEveryAtom(points[p]);
    }

    finished.store(0, std::memory_order_release);
    round.fetch_add(1, std::memory_order_release);
    while (finished.load(std::memory_order_acquire) != numberOfThreads)
    {
    }

    for (std::size_t p = 0; p < points.size(); ++p)
    {
      if (std::abs(answered[p] - expected[p]) > 1.0e-9) ++wrong;
    }
  }

  done.store(true, std::memory_order_release);
  for (std::thread &thread : threads) thread.join();

  EXPECT_EQ(wrong, 0uz) << "of " << numberOfRounds * points.size() << " answers";
}
