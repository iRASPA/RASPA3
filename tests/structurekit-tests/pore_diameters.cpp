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
