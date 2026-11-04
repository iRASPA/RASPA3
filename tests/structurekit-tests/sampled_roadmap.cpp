#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import sampled_structure;
import sampling_backend;
import mc_backend;
import sampled_roadmap;

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
