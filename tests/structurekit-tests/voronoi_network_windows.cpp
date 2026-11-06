#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import apollonius_network;
import channel_analysis;
import voronoi_pore_diameters;
import pore_window;

// An edge of the radical network is a straight segment, and its radius is the least clearance anywhere along
// it. That least is taken somewhere, so the edge has a narrowest point and a direction there like any other,
// and the windows that are measured across the bottleneck of a path can be measured on this network as on the
// one built from the Apollonius diagram. They used not to be: the builder kept the radius and not the place,
// and every window asked of it came back unmeasured.

namespace
{

struct Structure
{
  UnitCell box{12.0, 12.0, 12.0};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
};

Structure randomStructure(std::size_t seed, std::size_t count, double smallest, double spread)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  Structure structure;
  for (std::size_t i = 0; i < count; ++i)
  {
    structure.fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    structure.radii.push_back(smallest + spread * random.uniform());
  }
  return structure;
}

// The clearance at a point, over every atom and enough of its images, by looking at them all.
double clearanceAt(const VoronoiNetwork& network, const double3& point)
{
  double clearance = std::numeric_limits<double>::max();
  for (std::size_t j = 0; j < network.atomPositionsFractional.size(); ++j)
  {
    const double3 atom = network.unitCell.cell * network.atomPositionsFractional[j];
    for (int x = -2; x <= 2; ++x)
    {
      for (int y = -2; y <= 2; ++y)
      {
        for (int z = -2; z <= 2; ++z)
        {
          const double3 image = atom + network.unitCell.cell * double3(x, y, z);
          clearance = std::min(clearance, (point - image).length() - network.atomRadii[j]);
        }
      }
    }
  }
  return clearance;
}

// How far apart two points are once the lattice is taken out.
double distanceModuloLattice(const VoronoiNetwork& network, const double3& first, const double3& second)
{
  double3 fractional = network.unitCell.inverseCell * (first - second);
  fractional.x -= std::round(fractional.x);
  fractional.y -= std::round(fractional.y);
  fractional.z -= std::round(fractional.z);
  return (network.unitCell.cell * fractional).length();
}

double3 latticeShift(const VoronoiNetwork& network, const int3& delta)
{
  return network.unitCell.cell *
         double3(static_cast<double>(delta.x), static_cast<double>(delta.y), static_cast<double>(delta.z));
}

}  // namespace

TEST(voronoi_network, every_edge_knows_where_it_is_narrowest_and_is_that_narrow_there)
{
  const Structure structure = randomStructure(5, 25, 0.6, 1.4);
  const VoronoiNetwork network =
      VoronoiNetwork::create(structure.box, structure.fractionalPositions, structure.radii);
  ASSERT_FALSE(network.edges.empty());

  std::size_t exactlyAsNarrow = 0;
  for (std::size_t e = 0; e < network.edges.size(); ++e)
  {
    const VoronoiEdge& edge = network.edges[e];
    ASSERT_TRUE(edge.hasBottleneckGeometry) << "edge " << e;

    const double3 from = network.nodes[edge.from].position;
    const double3 to = network.nodes[edge.to].position + latticeShift(network, edge.delta);
    const double3 span = to - from;

    // On the segment, and pointing along it.
    const double along = double3::dot(edge.bottleneckPosition - from, span) / double3::dot(span, span);
    EXPECT_GE(along, -1.0e-9) << "edge " << e;
    EXPECT_LE(along, 1.0 + 1.0e-9) << "edge " << e;
    EXPECT_LT((from + along * span - edge.bottleneckPosition).length(), 1.0e-9) << "edge " << e;
    EXPECT_LT((edge.bottleneckDirection - span / span.length()).length(), 1.0e-9) << "edge " << e;

    // Nowhere on the segment is narrower than there.
    const double there = clearanceAt(network, edge.bottleneckPosition);
    for (int step = 0; step <= 20; ++step)
    {
      EXPECT_GE(clearanceAt(network, from + (step / 20.0) * span), there - 1.0e-9) << "edge " << e;
    }

    // And the edge is as narrow as that. Its radius is never more, and is less only where a vertex the
    // edge ends on was merged with another close by: the radius is then also held to the segment as each
    // cell saw it before the merging, which lies within the merging distance of this one.
    EXPECT_LE(edge.radius, there + 1.0e-9) << "edge " << e;
    EXPECT_GE(edge.radius, there - 0.04) << "edge " << e;
    if (std::abs(edge.radius - there) < 1.0e-9) ++exactlyAsNarrow;
  }
  EXPECT_GT(10 * exactlyAsNarrow, 9 * network.edges.size());

  // The builder stores each edge both ways round, one after the other, and the two are the same place seen
  // from either end.
  for (std::size_t e = 0; e + 1 < network.edges.size(); e += 2)
  {
    const VoronoiEdge& forward = network.edges[e];
    const VoronoiEdge& backward = network.edges[e + 1];
    ASSERT_EQ(forward.from, backward.to);
    ASSERT_EQ(forward.to, backward.from);
    const double3 seenFromTheOtherEnd = backward.bottleneckPosition + latticeShift(network, forward.delta);
    EXPECT_LT((seenFromTheOtherEnd - forward.bottleneckPosition).length(), 1.0e-9);
    EXPECT_LT((backward.bottleneckDirection + forward.bottleneckDirection).length(), 1.0e-12);
  }
}

TEST(voronoi_network, the_window_of_the_simple_cubic_lattice_is_measured_at_the_middle_of_a_cube_edge)
{
  // The same lattice the Apollonius network is held to: the bottleneck is the middle of an edge of the cube,
  // four atoms equally near, the free sphere passing between two of them.
  const double a = 5.0;
  const double r = 1.0;
  const UnitCell box(a, a, a);
  const VoronoiNetwork network = VoronoiNetwork::create(box, {double3(0.0, 0.0, 0.0)}, {r});
  const ChannelAnalysis channels = ChannelAnalysis::compute(network, 0.5);

  const std::vector<ChannelWindow> windows = channelWindows(network, channels);
  ASSERT_EQ(windows.size(), 1uz);
  const ChannelWindow& channel = windows[0];
  EXPECT_NEAR(channel.freeSphereDiameter, std::sqrt(2.0) * a - 2.0 * r, 1.0e-9);

  ASSERT_TRUE(channel.window.measured);
  EXPECT_NEAR(2.0 * channel.window.freeRadius, channel.freeSphereDiameter, 1.0e-9);
  EXPECT_NEAR(channel.window.smallestFreeWidth, channel.freeSphereDiameter, 1.0e-6);
  EXPECT_TRUE(channel.window.clipped);
  EXPECT_GE(channel.window.boundingAtoms, 4uz);

  // Half a cell from the atom along two axes and level with it along the third, which the passage runs along.
  std::array<double, 3> offset{channel.window.position.x / a, channel.window.position.y / a,
                               channel.window.position.z / a};
  std::array<double, 3> normal{std::abs(channel.window.normal.x), std::abs(channel.window.normal.y),
                               std::abs(channel.window.normal.z)};
  std::size_t halves = 0;
  for (std::size_t axis = 0; axis < 3; ++axis)
  {
    const double fraction = offset[axis] - std::floor(offset[axis]);
    if (std::abs(fraction - 0.5) < 1.0e-9)
    {
      ++halves;
      EXPECT_NEAR(normal[axis], 0.0, 1.0e-9);
    }
    else
    {
      EXPECT_LT(std::min(fraction, 1.0 - fraction), 1.0e-9);
      EXPECT_NEAR(normal[axis], 1.0, 1.0e-9);
    }
  }
  EXPECT_EQ(halves, 2uz);

  const PoreWindow ofTheStructure = freeSphereWindow(network);
  ASSERT_TRUE(ofTheStructure.measured);
  EXPECT_NEAR(2.0 * ofTheStructure.freeRadius, std::sqrt(2.0) * a - 2.0 * r, 1.0e-9);
}

TEST(voronoi_network, on_equal_radii_its_windows_are_those_of_the_apollonius_network)
{
  // With one radius the two diagrams are the same diagram, so the widest path is limited at the same place
  // and the window across it is the same window.
  for (std::size_t seed : {2uz, 9uz, 14uz})
  {
    const Structure structure = randomStructure(seed, 22, 1.1, 0.0);
    const VoronoiNetwork radical =
        VoronoiNetwork::create(structure.box, structure.fractionalPositions, structure.radii);
    const ApolloniusPoreNetwork apollonius =
        ApolloniusPoreNetwork::create(structure.box, structure.fractionalPositions, structure.radii);

    const PoreWindow mine = freeSphereWindow(radical);
    const PoreWindow theirs = freeSphereWindow(apollonius.network);
    ASSERT_TRUE(theirs.measured) << "seed " << seed;
    ASSERT_TRUE(mine.measured) << "seed " << seed;

    EXPECT_NEAR(mine.freeRadius, theirs.freeRadius, 1.0e-4) << "seed " << seed;
    EXPECT_LT(distanceModuloLattice(radical, mine.position, theirs.position), 1.0e-2) << "seed " << seed;
    EXPECT_GT(std::abs(double3::dot(mine.normal, theirs.normal)), 1.0 - 1.0e-4) << "seed " << seed;
    EXPECT_NEAR(mine.smallestFreeWidth, theirs.smallestFreeWidth, 2.0e-2) << "seed " << seed;
    EXPECT_NEAR(mine.minorAxis, theirs.minorAxis, 5.0e-2) << "seed " << seed;
  }
}

TEST(voronoi_network, on_unequal_radii_its_windows_bracket_its_own_free_sphere)
{
  // What holds of any window, on atoms of several sizes where the channels are of no particular shape: the
  // free sphere is one of the chords through the bottleneck, so the shortest of them is at least that.
  const Structure structure = randomStructure(7, 25, 0.8, 0.9);
  const VoronoiNetwork network =
      VoronoiNetwork::create(structure.box, structure.fractionalPositions, structure.radii);
  const ChannelAnalysis channels = ChannelAnalysis::compute(network, 1.0);
  const std::vector<ChannelWindow> windows = channelWindows(network, channels);
  ASSERT_FALSE(windows.empty());

  const PoreDiameters diameters = PoreDiameters::compute(network);
  for (const ChannelWindow& channel : windows)
  {
    EXPECT_LE(channel.freeSphereDiameter, diameters.freeSphereDiameter + 1.0e-9);
    ASSERT_TRUE(channel.window.measured);
    EXPECT_NEAR(2.0 * channel.window.freeRadius, channel.freeSphereDiameter, 1.0e-9);
    EXPECT_GE(channel.window.smallestFreeWidth, channel.freeSphereDiameter - 1.0e-6);
    EXPECT_GE(channel.window.largestFreeWidth, channel.window.smallestFreeWidth - 1.0e-9);
    EXPECT_GE(channel.window.majorAxis, channel.window.minorAxis - 1.0e-9);
    EXPECT_LE(channel.window.majorAxis, channel.window.largestFreeWidth + 1.0e-6);
    EXPECT_NEAR(channel.window.normal.length(), 1.0, 1.0e-9);
  }

  const PoreWindow ofTheStructure = freeSphereWindow(network);
  ASSERT_TRUE(ofTheStructure.measured);
  EXPECT_NEAR(2.0 * ofTheStructure.freeRadius, diameters.freeSphereDiameter, 1.0e-9);
}
