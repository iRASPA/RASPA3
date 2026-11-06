#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import sampled_structure;
import sampling_backend;
import mc_backend;
import sampled_roadmap;
import voronoi_pore_diameters;

// The roadmap draws its points afresh on every build, so what is checked here is what has to hold of any roadmap,
// whichever points it drew.

namespace
{

SampledStructure scattered(const UnitCell &cell, std::size_t atoms, double radius, std::size_t seed)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  SampledStructure structure;
  structure.name = "sampled_roadmap_test";
  structure.unitCell = cell;
  for (std::size_t i = 0; i < atoms; ++i)
  {
    structure.positions.push_back(cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
    structure.radii.push_back(radius);
  }
  return structure;
}

}  // namespace


// A hop is stored both ways round, and each way round reads it in the frame of the node it leaves from: the far node's
// frame is a lattice translation away, and the direction of travel is the other way.
TEST(sampled_roadmap, each_hop_is_read_from_the_frame_of_the_node_it_leaves)
{
  const UnitCell cell(16.0, 17.0, 18.0);
  SampledStructure structure = scattered(cell, 30, 1.2, 3);
  SampledRoadmap roadmap = SampledRoadmap::build(structure, samplingBackendCPU(), 3000, 50);

  ASSERT_GT(roadmap.network.edges.size(), 0uz);
  ASSERT_EQ(roadmap.network.edges.size() % 2, 0uz);

  std::size_t acrossTheBoundary = 0;
  for (std::size_t i = 0; i + 1 < roadmap.network.edges.size(); i += 2)
  {
    const VoronoiEdge &forward = roadmap.network.edges[i];
    const VoronoiEdge &backward = roadmap.network.edges[i + 1];
    ASSERT_EQ(forward.from, backward.to);
    ASSERT_EQ(forward.to, backward.from);

    const double3 shift = cell.cell * double3(static_cast<double>(forward.delta.x), static_cast<double>(forward.delta.y),
                                              static_cast<double>(forward.delta.z));
    if (shift.length() > 0.0) ++acrossTheBoundary;

    EXPECT_NEAR((backward.bottleneckPosition - (forward.bottleneckPosition - shift)).length(), 0.0, 1.0e-12);
    EXPECT_NEAR((backward.bottleneckDirection + forward.bottleneckDirection).length(), 0.0, 1.0e-12);
    EXPECT_NEAR(forward.bottleneckDirection.length(), 1.0, 1.0e-3);
  }
  EXPECT_GT(acrossTheBoundary, 0uz);
}

// The bottleneck of a hop is where its radius was measured. So the room at that point is that radius, on every hop
// there is: the ones between samples, the ones bent, and the ones between pockets that go in as the straight line
// left them, of which there are many in a cell with many pockets.
TEST(sampled_roadmap, the_room_at_a_hops_bottleneck_is_its_radius)
{
  const UnitCell cell(16.0, 17.0, 18.0);
  SampledStructure structure = scattered(cell, 30, 1.2, 3);
  SampledRoadmap roadmap = SampledRoadmap::build(structure, samplingBackendCPU(), 3000, 50);

  ASSERT_GT(roadmap.numberOfPocketCentres, 3uz);
  std::size_t betweenPockets = 0;
  std::size_t wrong = 0;
  double worst = 0.0;
  for (const VoronoiEdge &edge : roadmap.network.edges)
  {
    const double difference = std::abs(structure.clearance(edge.bottleneckPosition) - edge.radius);
    if (edge.from >= roadmap.firstPocketNode && edge.to >= roadmap.firstPocketNode) ++betweenPockets;
    worst = std::max(worst, difference);
    if (difference > 1.0e-9) ++wrong;
  }
  EXPECT_GT(betweenPockets, 2uz * 16uz) << "more hops than the ones that get bent";
  EXPECT_EQ(wrong, 0uz) << "worst by " << worst << " Å, of " << roadmap.network.edges.size();
}

// In a cell narrower than the link distance a sample has neighbours among its own images, and a pair of samples can
// be in reach by more than one image. Every one of them is a hop, measured like any other, and each is there once.
TEST(sampled_roadmap, a_thin_cell_links_a_sample_to_its_own_images_and_to_every_image_in_reach)
{
  const UnitCell cell(2.0, 20.0, 20.0);
  SampledStructure structure;
  structure.name = "sampled_roadmap_test";
  structure.unitCell = cell;
  structure.positions.push_back(double3(1.0, 10.0, 10.0));
  structure.radii.push_back(0.2);

  SampledRoadmap roadmap = SampledRoadmap::build(structure, samplingBackendCPU(), 200, 20);
  ASSERT_GT(roadmap.linkDistance, 2.0);

  std::size_t ownImage = 0;
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> hopsOfPair;
  std::set<std::tuple<std::size_t, std::size_t, int, int, int>> seen;
  for (std::size_t i = 0; i < roadmap.network.edges.size(); i += 2)
  {
    const VoronoiEdge &edge = roadmap.network.edges[i];
    if (edge.from >= roadmap.firstPocketNode || edge.to >= roadmap.firstPocketNode) continue;

    const double3 reach = roadmap.network.nodes[edge.to].position - roadmap.network.nodes[edge.from].position +
                          cell.cell * double3(static_cast<double>(edge.delta.x), static_cast<double>(edge.delta.y),
                                              static_cast<double>(edge.delta.z));
    EXPECT_LT(reach.length(), roadmap.linkDistance);
    EXPECT_GT(reach.length(), 0.0);
    EXPECT_NEAR(edge.length, reach.length(), 1.0e-9);

    EXPECT_TRUE(seen.insert({edge.from, edge.to, edge.delta.x, edge.delta.y, edge.delta.z}).second);
    if (edge.from == edge.to)
    {
      ++ownImage;
      EXPECT_FALSE(seen.contains({edge.from, edge.to, -edge.delta.x, -edge.delta.y, -edge.delta.z}))
          << "one hop, stored twice";
    }
    ++hopsOfPair[{edge.from, edge.to}];
  }
  EXPECT_GT(ownImage, 0uz);

  std::size_t reachedTwice = 0;
  for (const auto &[pair, hops] : hopsOfPair)
    if (pair.first != pair.second && hops > 1) ++reachedTwice;
  EXPECT_GT(reachedTwice, 0uz);
}


// A pocket is tried against the nearest copy of every other pocket and against that copy's neighbours, and
// against its own next repeats, in the compact cell of the lattice and not in the cell as given. A cell may be
// given by vectors that lean along one another, and then the copies one cell either way along each of those
// vectors are not the near ones: in a cubic lattice described with its second vector leaning two lengths of the
// first, the repeat along y is two cells back along the first and one along the second, which a search of one
// cell either way never tries. The way from a pocket to its own repeat is the straight channel through the
// crystal, and the middle of it is the window: left to the hops between scattered points, the window is found
// only as nearly as the points happen to fall.
//
// Simple cubic, one atom to a cell: one pocket at the body centre, a window in each face, and a free sphere
// that is the window's inscribed circle exactly, sqrt(2) a - 2 r, in whatever cell the lattice is given.
TEST(sampled_roadmap, a_pocket_reaches_its_own_next_repeat_in_a_cell_given_by_leaning_vectors)
{
  const double edge = 8.0;
  const double radius = 2.5;
  const double window = edge / std::sqrt(2.0) - radius;

  const double3 a(edge, 0.0, 0.0);
  const double3 b(0.0, edge, 0.0);
  const double3 c(0.0, 0.0, edge);
  const std::vector<std::pair<std::string, UnitCell>> cells{
      {"compact", UnitCell(double3x3(a, b, c))},
      {"leaning", UnitCell(double3x3(a, b + 2.0 * a, c + 2.0 * (b + 2.0 * a)))},
  };

  for (const auto &[name, cell] : cells)
  {
    SampledStructure structure;
    structure.name = "sampled_roadmap_test";
    structure.unitCell = cell;
    structure.positions.push_back(double3(0.0, 0.0, 0.0));
    structure.radii.push_back(radius);

    // The points are drawn afresh on every build, so what is asked is what holds of any of them. The roomiest
    // pocket is the one at the body centre and is the first of the pocket nodes; a walk that stalled on an edge
    // of the cube now and then leaves a second, which is nothing to do with this.
    const SampledRoadmap roadmap = SampledRoadmap::build(structure, samplingBackendCPU(), 2000, 2000);
    ASSERT_GE(roadmap.numberOfPocketCentres, 1uz) << name;
    const std::size_t pocket = roadmap.firstPocketNode;
    EXPECT_NEAR(roadmap.network.nodes[pocket].radius, 0.5 * std::sqrt(3.0) * edge - radius, 1.0e-2) << name;

    // That pocket's hops to its own repeats, by which way each runs.
    std::array<bool, 3> alongAxis{false, false, false};
    for (const VoronoiEdge &hop : roadmap.network.edges)
    {
      if (hop.from != pocket || hop.to != pocket) continue;
      const double3 reach = cell.cell * double3(static_cast<double>(hop.delta.x), static_cast<double>(hop.delta.y),
                                                static_cast<double>(hop.delta.z));
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const double3 repeat = (axis == 0) ? a : ((axis == 1) ? b : c);
        if ((reach - repeat).length() < 1.0e-9 || (reach + repeat).length() < 1.0e-9)
        {
          alongAxis[axis] = true;

          // No way through a face is wider than the face's window, and the line from the pocket to its repeat
          // misses the middle of the window only by as much as the walk uphill left the pocket off the body
          // centre.
          EXPECT_LE(hop.radius, window + 1.0e-9) << name << ", along axis " << axis;
          EXPECT_GT(hop.radius, window - 5.0e-2) << name << ", along axis " << axis;
        }
      }
    }
    EXPECT_TRUE(alongAxis[0]) << name << ": no hop from the pocket to its repeat along x";
    EXPECT_TRUE(alongAxis[1]) << name << ": no hop from the pocket to its repeat along y";
    EXPECT_TRUE(alongAxis[2]) << name << ": no hop from the pocket to its repeat along z";

    const PoreDiameters diameters = PoreDiameters::compute(roadmap.network);
    EXPECT_LE(diameters.freeSphereDiameter, 2.0 * window + 1.0e-9) << name;
    EXPECT_GT(diameters.freeSphereDiameter, 2.0 * window - 1.0e-1) << name;
  }
}

// Which of the hops between pockets are tried again with the line allowed to bend: the best few of every
// pocket, counting a hop for both of the pockets it joins. Counted only for the roomier of its two ends, a hop
// was bent or not by how it stood among that pocket's hops alone, and a small pocket beside a roomy one that has
// sixteen wider ways out got none of its own bent, though the one way it has may be the one that needs it.
TEST(sampled_roadmap, the_hops_bent_are_the_best_of_each_pocket_whichever_end_it_is)
{
  // Pocket 0 is roomy and has twenty wide ways to pockets 1 to 20; pocket 21 is small, and its one way out is
  // a narrow one to pocket 0, the narrowest of all of them.
  std::vector<PocketHop> hops;
  for (std::size_t other = 1; other <= 20; ++other)
  {
    hops.push_back(PocketHop{.from = 0, .to = other, .radius = 3.0 - 0.01 * static_cast<double>(other)});
  }
  hops.push_back(PocketHop{.from = 0, .to = 21, .radius = -0.2});
  const std::size_t theSmallPocketsOnlyWay = hops.size() - 1;

  const std::vector<std::uint8_t> bent = hopsToBend(22, hops, 16);
  ASSERT_EQ(bent.size(), hops.size());

  // Every one of them is the best its far end has, so every one is bent: the twenty for pockets 1 to 20, and the
  // narrow one for pocket 21.
  for (std::size_t h = 0; h < hops.size(); ++h) EXPECT_EQ(bent[h], 1) << "hop " << h;
  EXPECT_EQ(bent[theSmallPocketsOnlyWay], 1);

  // A pocket with more ways than are bent has only its best ones bent on its own account.
  std::vector<PocketHop> spokes;
  for (std::size_t n = 0; n < 40; ++n)
  {
    // forty ways from pocket 0 to its own repeats and to pocket 1, each wider than the last
    spokes.push_back(PocketHop{.from = 0, .to = n % 2, .radius = 0.1 * static_cast<double>(n)});
  }
  const std::vector<std::uint8_t> chosen = hopsToBend(2, spokes, 16);
  // The sixteen widest for pocket 0 are the last sixteen. Pocket 1 has the twenty odd ones, and its sixteen
  // widest of those reach down to the ninth: so from the ninth on the odd ones are bent, and from the
  // twenty-fourth on all of them.
  for (std::size_t n = 0; n < 40; ++n)
  {
    const bool forTheRoomyEnd = n >= 24;
    const bool forTheOtherEnd = (n % 2 == 1) && n >= 9;
    EXPECT_EQ(chosen[n] != 0, forTheRoomyEnd || forTheOtherEnd) << "way " << n;
  }

  // And none at all where none are asked for, or there are none.
  EXPECT_TRUE(hopsToBend(3, {}, 16).empty());
  for (std::uint8_t one : hopsToBend(22, hops, 0)) EXPECT_EQ(one, 0);
}
