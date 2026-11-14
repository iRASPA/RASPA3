#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import voronoi_network;
import pore_diameters;

// The widest percolating path and what hangs off it. The network here is written by hand, a few nodes and
// edges, so that which edge closes the loop and which nodes are reachable at its radius is a matter of reading
// the picture.

namespace
{

VoronoiNetwork tinyNetwork()
{
  VoronoiNetwork network;
  network.unitCell = UnitCell(10.0, 10.0, 10.0);
  auto node = [&](double3 position, double radius)
  {
    VoronoiNode n;
    n.position = position;
    n.fractional = position * 0.1;
    n.radius = radius;
    n.maximalRadius = radius;
    n.maximalPosition = position;
    network.nodes.push_back(n);
  };
  auto edge = [&](std::size_t from, std::size_t to, int3 delta, double radius)
  {
    VoronoiEdge forward;
    forward.from = from;
    forward.to = to;
    forward.delta = delta;
    forward.radius = radius;
    forward.length = 1.0;
    VoronoiEdge backward = forward;
    backward.from = to;
    backward.to = from;
    backward.delta = int3(0, 0, 0) - delta;
    network.edges.push_back(forward);
    network.edges.push_back(backward);
  };
  // Node 0 is a channel along x closing on its own image; node 1 is a wider room hanging off it through a
  // neck exactly as wide as the channel's own bottleneck.
  node(double3(1.0, 5.0, 5.0), 1.0);
  node(double3(5.0, 5.0, 5.0), 3.0);
  edge(0, 0, int3(1, 0, 0), 2.0);  // the loop, listed first so that it is met first among the ties
  edge(0, 1, int3(0, 0, 0), 2.0);  // the tied neck to the room
  return network;
}

}  // namespace


// Two edges of one radius, one closing the loop and one leading to a wider room: the room is reachable
// without passing anything narrower than the bottleneck, so it is on the path and sets Dif. Which of the two
// edges the sort happens to put first cannot be allowed to decide that.
TEST(pore_diameters, a_room_behind_a_neck_as_wide_as_the_bottleneck_is_on_the_path)
{
  VoronoiNetwork network = tinyNetwork();
  std::vector<std::size_t> all = {0, 1};

  PercolatingPath path = widestPercolatingPath(network, all);
  ASSERT_TRUE(path.percolates);
  EXPECT_NEAR(path.radius, 2.0, 1.0e-12);
  EXPECT_EQ(path.componentNodes.size(), 2uz) << "the room joins at the bottleneck's own radius";

  PoreDiameters diameters = PoreDiameters::compute(network);
  EXPECT_NEAR(diameters.freeSphereDiameter, 4.0, 1.0e-12);
  EXPECT_NEAR(diameters.includedAlongFreePathDiameter, 6.0, 1.0e-12) << "Dif is the room's, not the channel's";
  EXPECT_NEAR(diameters.includedSphereDiameter, 6.0, 1.0e-12);
}


// A neck narrower than the bottleneck, by however little, keeps the room off the path.
TEST(pore_diameters, a_room_behind_a_narrower_neck_is_not_on_the_path)
{
  VoronoiNetwork network = tinyNetwork();
  network.edges[2].radius = 2.0 - 1.0e-9;
  network.edges[3].radius = 2.0 - 1.0e-9;
  std::vector<std::size_t> all = {0, 1};

  PercolatingPath path = widestPercolatingPath(network, all);
  ASSERT_TRUE(path.percolates);
  EXPECT_EQ(path.componentNodes.size(), 1uz);
  EXPECT_NEAR(PoreDiameters::compute(network).includedAlongFreePathDiameter, 2.0, 1.0e-12);
}


// Two channels that never meet, with the same bottleneck. A sphere of that width goes down either, so both are the
// free path, and Dif is the larger of what they hold. It used to be the room in whichever of them closed its loop
// first, which is a question about the order of the edges and not about the crystal.
TEST(pore_diameters, disconnected_channels_tied_at_the_bottleneck_are_all_on_the_path)
{
  auto network = [](bool roomInFirst)
  {
    VoronoiNetwork tied = tinyNetwork();
    tied.nodes.clear();
    tied.edges.clear();
    auto node = [&](double3 position, double radius)
    {
      VoronoiNode n;
      n.position = position;
      n.fractional = position * 0.1;
      n.radius = radius;
      n.maximalRadius = radius;
      n.maximalPosition = position;
      tied.nodes.push_back(n);
    };
    auto edge = [&](std::size_t from, std::size_t to, int3 delta, double radius)
    {
      VoronoiEdge forward;
      forward.from = from;
      forward.to = to;
      forward.delta = delta;
      forward.radius = radius;
      forward.length = 1.0;
      VoronoiEdge backward = forward;
      backward.from = to;
      backward.to = from;
      backward.delta = int3(0, 0, 0) - delta;
      tied.edges.push_back(forward);
      tied.edges.push_back(backward);
    };
    node(double3(1.0, 2.0, 5.0), roomInFirst ? 4.0 : 2.5);
    node(double3(1.0, 8.0, 5.0), roomInFirst ? 2.5 : 4.0);
    node(double3(5.0, 5.0, 5.0), 9.0);  // a pocket wider than either, behind a neck narrower than the bottleneck
    edge(0, 0, int3(1, 0, 0), 2.0);
    edge(1, 1, int3(1, 0, 0), 2.0);
    edge(0, 2, int3(0, 0, 0), 1.0);
    return tied;
  };

  for (bool roomInFirst : {true, false})
  {
    VoronoiNetwork tied = network(roomInFirst);
    std::vector<std::size_t> all = {0, 1, 2};

    PercolatingPath path = widestPercolatingPath(tied, all);
    ASSERT_TRUE(path.percolates);
    EXPECT_NEAR(path.radius, 2.0, 1.0e-12);
    EXPECT_EQ(path.componentNodes, (std::vector<std::size_t>{0, 1})) << "room in first: " << roomInFirst;

    PoreDiameters diameters = PoreDiameters::compute(tied);
    EXPECT_NEAR(diameters.freeSphereDiameter, 4.0, 1.0e-12);
    EXPECT_NEAR(diameters.includedAlongFreePathDiameter, 8.0, 1.0e-12) << "room in first: " << roomInFirst;
    EXPECT_NEAR(diameters.includedSphereDiameter, 18.0, 1.0e-12);
  }
}

// A channel that only reaches its bottleneck's width later is not on the path, and neither is a piece as wide as the
// bottleneck that does not run through: it takes percolating at that width to be on the free path.
TEST(pore_diameters, a_piece_tied_at_the_bottleneck_that_does_not_run_through_is_not_on_the_path)
{
  VoronoiNetwork network = tinyNetwork();
  VoronoiNode pocket = network.nodes[1];
  pocket.maximalRadius = 7.0;
  network.nodes.push_back(pocket);
  network.nodes.push_back(pocket);
  VoronoiEdge across = network.edges[0];
  across.from = 2;
  across.to = 3;
  across.delta = int3(0, 0, 0);
  across.radius = 2.0;
  network.edges.push_back(across);
  std::swap(across.from, across.to);
  network.edges.push_back(across);

  std::vector<std::size_t> all = {0, 1, 2, 3};
  PercolatingPath path = widestPercolatingPath(network, all);
  ASSERT_TRUE(path.percolates);
  EXPECT_EQ(path.componentNodes, (std::vector<std::size_t>{0, 1}));
}

// A closed cage: its node is joined to its own images only through edges buried in the atoms (negative clearance,
// the network keeps them). No sphere gets through, so the structure does not percolate: Df and Dif are nought, and
// not a negative Df and the cage's own diameter.
TEST(pore_diameters, a_cage_joined_only_through_buried_edges_does_not_percolate)
{
  VoronoiNetwork network;
  network.unitCell = UnitCell(10.0, 10.0, 10.0);
  VoronoiNode cage;
  cage.position = double3(5.0, 5.0, 5.0);
  cage.fractional = double3(0.5, 0.5, 0.5);
  cage.radius = 3.0;
  cage.maximalRadius = 3.0;
  cage.maximalPosition = cage.position;
  network.nodes.push_back(cage);
  for (const int3 delta : {int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1)})
  {
    network.edges.push_back(VoronoiEdge{0, 0, delta, -0.4, 10.0});
    network.edges.push_back(VoronoiEdge{0, 0, int3(0, 0, 0) - delta, -0.4, 10.0});
  }

  const std::vector<std::size_t> nodes{0};
  EXPECT_FALSE(widestPercolatingPath(network, nodes).percolates);
  const PoreDiameters diameters = PoreDiameters::compute(network);
  EXPECT_NEAR(diameters.includedSphereDiameter, 6.0, 1e-12);
  EXPECT_EQ(diameters.freeSphereDiameter, 0.0);
  EXPECT_EQ(diameters.includedAlongFreePathDiameter, 0.0);
}
