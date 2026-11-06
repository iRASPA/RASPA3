#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import voronoi_network;
import pore_accessibility;
import apollonius_accessibility;
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


namespace
{

// The answers a structure with no atoms has to give, wherever it came from: every point is free, to any
// distance, and next to nothing. Whether it is in a pore is for the nodes to say, and there are none.
void expectTheAnswersOfNoAtoms(const PoreAccessibility& accessibility)
{
  for (const double3& point : {double3(0.0, 0.0, 0.0), double3(3.1, 7.2, 9.9), double3(-14.0, 25.0, 3.0)})
  {
    EXPECT_EQ(accessibility.clearance(point), std::numeric_limits<double>::max());

    const PointClassification classification = accessibility.classify(point);
    EXPECT_FALSE(classification.inside);
    EXPECT_FALSE(classification.accessible);
    EXPECT_TRUE(classification.resample);
    EXPECT_EQ(classification.poreId, -1);

    EXPECT_FALSE(accessibility.containingNode(point).has_value());
    EXPECT_FALSE(accessibility.provablyAccessible(point));
    EXPECT_FALSE(accessibility.overlapsAtom(point, 0));
    EXPECT_TRUE(accessibility.neighbourAtoms(point, 5.0).empty());
    EXPECT_TRUE(accessibility.neighbourAtomImages(point, 5.0).empty());
  }
}

}  // namespace

TEST(empty_structures, a_network_of_no_atoms_answers_every_question_and_returns)
{
  // The radical diagram refuses a cell with no atoms, but a network need not have come from it, and the
  // classifier takes whatever network it is handed. Its searches walk outward over shells of bins until
  // they have found an atom and no nearer one can follow, which with no atom to find is never.
  VoronoiNetwork network;
  network.unitCell = UnitCell(10.0, 11.0, 12.0);
  const PoreAccessibility accessibility =
      PoreAccessibility::createFromNetwork(network, PoreAccessibility::Metric::Clearance, 1.0);

  EXPECT_TRUE(accessibility.atomPositions.empty());
  expectTheAnswersOfNoAtoms(accessibility);

  const PoreAccessibility inPowers = PoreAccessibility::createFromNetwork(network, PoreAccessibility::Metric::Power);
  expectTheAnswersOfNoAtoms(inPowers);
}

TEST(empty_structures, a_classifier_that_was_never_built_answers_the_same)
{
  // Nothing at all: no cell, no bins to walk over.
  const PoreAccessibility accessibility;
  expectTheAnswersOfNoAtoms(accessibility);
}

TEST(empty_structures, the_apollonius_route_on_no_atoms_builds_a_classifier_that_returns)
{
  // That builder does not refuse a cell with no atoms: it makes a network of no nodes and hands it straight
  // to the classifier, which is how a structure of nothing gets as far as being asked about a point.
  const UnitCell box(10.0, 10.0, 10.0);
  const ApolloniusAccessibility built = ApolloniusAccessibility::create(box, {}, {}, 1.0);

  EXPECT_TRUE(built.accessibility.network.nodes.empty());
  EXPECT_TRUE(built.accessibility.atomPositions.empty());
  expectTheAnswersOfNoAtoms(built.accessibility);
}
