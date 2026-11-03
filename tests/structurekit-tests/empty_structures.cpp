#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import voronoi_network;
import pore_accessibility;
import pore_window;

// What the diagram routes do with nothing to work on: a cell with no atoms, and a build asked for the
// unweighted diagram by being given no radii. The first is refused outright. The second has no largest radius
// to read off a list, and used to be read off the end of an empty one, and off its end again wherever an atom's
// radius was looked up.

TEST(empty_structures, a_cell_with_no_atoms_is_refused_by_the_diagram)
{
  // There is no diagram of nothing, and the refusal comes from the diagram itself, before any radius is
  // looked for: a structure with no atoms never reaches the largest radius.
  const UnitCell box(10.0, 10.0, 10.0);
  EXPECT_THROW(PoreAccessibility::create(box, {}, {}, 1.0), std::runtime_error);
}


TEST(empty_structures, an_unweighted_build_carries_a_zero_radius_per_atom)
{
  const UnitCell box(10.0, 10.0, 10.0);
  const std::vector<double3> positions{double3(0.1, 0.2, 0.3), double3(0.6, 0.7, 0.8), double3(0.4, 0.9, 0.1)};
  VoronoiNetwork network = VoronoiNetwork::create(box, positions, {});

  ASSERT_EQ(network.atomRadii.size(), positions.size());
  for (double radius : network.atomRadii) EXPECT_EQ(radius, 0.0);
  EXPECT_FALSE(network.nodes.empty());

  // The two readers of the largest radius, on that network.
  PoreAccessibility accessibility = PoreAccessibility::createFromNetwork(network, PoreAccessibility::Metric::Power, 0.0);
  EXPECT_EQ(accessibility.maximumAtomRadius, 0.0);

  PoreWindow window = PoreWindow::measure(network, double3(5.0, 5.0, 5.0), double3(0.0, 0.0, 1.0));
  EXPECT_NEAR((window.position - double3(5.0, 5.0, 5.0)).length(), 0.0, 1.0e-12);
}
