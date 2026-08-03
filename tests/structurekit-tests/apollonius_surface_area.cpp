#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import pore_accessibility;
import voronoi_surface_area;
import apollonius_accessibility;
import apollonius_surface_area;

// The measured accessible surface area: the boundary of the union of the probe-inflated atoms, taken
// patch by patch instead of by throwing points at it.
//
// What can be checked exactly is the area, so that is what most of these do. A configuration of a few
// spheres has a closed form, and the routine has to return it to round-off rather than to a tolerance
// chosen to accommodate it. The split into channel and pocket is the classifier's answer and not this
// routine's, so it is checked only against the sampled estimate, which asks the same classifier.

namespace
{

// A classifier over bare geometry: the atoms and the cell, and no network at all. Every point then comes
// back undecided, which is what is wanted here -- the area is all in `total()` and none of it has been
// through a diagram -- and it keeps these cases free of the question of whether a diagram of four atoms
// in a large box has any vertices.
PoreAccessibility bareGeometry(const UnitCell& box, const std::vector<double3>& fractionalPositions,
                                 const std::vector<double>& radii)
{
  VoronoiNetwork network;
  network.unitCell = box;
  network.atomRadii = radii;
  for (const double3& fractional : fractionalPositions)
  {
    network.atomPositionsFractional.push_back(double3::fract(fractional));
  }
  network.atomNodeVectors.assign(fractionalPositions.size(), {});
  return PoreAccessibility::createFromNetwork(std::move(network), PoreAccessibility::Metric::Clearance);
}


// Area of the union of two overlapping spheres: each loses the cap its neighbour cuts out of it, and the
// cap's area is set by the polar angle of the circle of intersection.
double unionAreaOfTwoSpheres(double firstRadius, double secondRadius, double distance)
{
  double firstCosine =
      (firstRadius * firstRadius + distance * distance - secondRadius * secondRadius) / (2.0 * firstRadius * distance);
  double secondCosine = (secondRadius * secondRadius + distance * distance - firstRadius * firstRadius) /
                        (2.0 * secondRadius * distance);
  return 2.0 * std::numbers::pi * firstRadius * firstRadius * (1.0 + firstCosine) +
         2.0 * std::numbers::pi * secondRadius * secondRadius * (1.0 + secondCosine);
}

}  // namespace


// A sphere on its own, with nothing near it, has the area of a sphere. It is the one case where the
// latitude integral runs over the whole sphere with no circle to cut it, and it fixes the quadrature
// and the frame independently of everything else.
TEST(apollonius_surface_area, lone_sphere)
{
  double a = 40.0;
  UnitCell box(a, a, a);
  const double radius = 1.7;

  PoreAccessibility geometry = bareGeometry(box, {double3(0.5, 0.5, 0.5)}, {radius});
  MeasuredPatches area = exactAccessibleSurfaceAreaByPore(geometry);

  EXPECT_NEAR(area.total(), 4.0 * std::numbers::pi * radius * radius, 1.0e-11);
}


// Two spheres against the closed form, over radii and separations that put the circle of intersection
// on either side of both centres. The tolerance is the one the closed form itself deserves.
TEST(apollonius_surface_area, two_spheres_against_the_closed_form)
{
  double a = 40.0;
  UnitCell box(a, a, a);

  const std::array<std::array<double, 3>, 5> cases = {{{2.0, 1.0, 2.5},
                                                       {3.0, 2.5, 1.0},
                                                       {1.5, 1.5, 2.9},
                                                       {4.0, 0.6, 3.5},
                                                       {2.2, 2.0, 0.3}}};

  for (const std::array<double, 3>& one : cases)
  {
    double firstRadius = one[0];
    double secondRadius = one[1];
    double distance = one[2];

    double3 first = double3(0.5, 0.5, 0.5);
    double3 second = first + double3(distance / a, 0.0, 0.0);
    PoreAccessibility geometry = bareGeometry(box, {first, second}, {firstRadius, secondRadius});

    MeasuredPatches area = exactAccessibleSurfaceAreaByPore(geometry);
    double expected = unionAreaOfTwoSpheres(firstRadius, secondRadius, distance);
    EXPECT_NEAR(area.total(), expected, 1.0e-10 * expected);
  }
}


// A sphere inside another carries no surface, and the one that swallowed it keeps all of its own. Both
// halves of that matter: an atom whose sphere is buried has to be left out of the total, and it must not
// take a bite out of its neighbour on the way.
TEST(apollonius_surface_area, a_buried_sphere_carries_no_area)
{
  double a = 40.0;
  UnitCell box(a, a, a);
  const double large = 3.0;
  const double small = 0.5;

  PoreAccessibility geometry =
      bareGeometry(box, {double3(0.5, 0.5, 0.5), double3(0.5 + 1.0 / a, 0.5, 0.5)}, {large, small});
  MeasuredPatches area = exactAccessibleSurfaceAreaByPore(geometry);

  EXPECT_NEAR(area.total(), 4.0 * std::numbers::pi * large * large, 1.0e-11);
}


// One atom in a cell too small to hold it: its own periodic images cut into it, and the area left is
// known. With the cell edge between the radius and the radius times the square root of two, the six
// face images each take a cap and no two of those caps meet, so what is left is the sphere less six
// caps of height r - a/2. An atom is a neighbour of its own images like any other, and this is what says
// so; getting it wrong shows up nowhere else, because in a cell of ordinary size it never happens.
TEST(apollonius_surface_area, an_atom_is_cut_by_its_own_images)
{
  const double a = 10.0;
  const double radius = 6.0;  // between a/2 and a/sqrt(2), so exactly the six face images cut it
  UnitCell box(a, a, a);

  PoreAccessibility geometry = bareGeometry(box, {double3(0.5, 0.5, 0.5)}, {radius});
  MeasuredPatches area = exactAccessibleSurfaceAreaByPore(geometry);

  double capHeight = radius - 0.5 * a;
  double expected = 4.0 * std::numbers::pi * radius * radius -
                    6.0 * 2.0 * std::numbers::pi * radius * capHeight;
  EXPECT_NEAR(area.total(), expected, 1.0e-10 * expected);
}


// The answer must not depend on how finely the pieces of the latitude integral are cut up, which is the
// statement that the pieces are the right ones: the integrand is smooth inside each of them, so the
// quadrature has nothing left to resolve. A random packing is used because it puts the corners of the
// patches in general position, and a lattice with equal radii is added because it puts them all in the
// same place, three and four spheres meeting in one point.
TEST(apollonius_surface_area, refining_the_quadrature_does_not_move_the_answer)
{
  double a = 14.0;
  UnitCell box(a, a, a);

  RandomNumber random{std::optional<std::size_t>(31)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 40; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(1.4 + 1.4 * random.uniform());
  }

  PoreAccessibility packing = bareGeometry(box, fractionalPositions, radii);
  double coarse = exactAccessibleSurfaceAreaByPore(packing, 1).total();
  double fine = exactAccessibleSurfaceAreaByPore(packing, 5).total();
  EXPECT_GT(coarse, 0.0);
  EXPECT_NEAR(coarse, fine, 1.0e-8 * coarse);

  // Body-centred cubic of equal spheres, every contact degenerate.
  std::vector<double3> lattice;
  std::vector<double> equalRadii;
  for (double x : {0.0, 0.5})
  {
    for (double y : {0.0, 0.5})
    {
      for (double z : {0.0, 0.5})
      {
        lattice.push_back(double3(x, y, z));
        equalRadii.push_back(4.2);
      }
    }
  }
  PoreAccessibility crystal = bareGeometry(box, lattice, equalRadii);
  double coarseCrystal = exactAccessibleSurfaceAreaByPore(crystal, 1).total();
  double fineCrystal = exactAccessibleSurfaceAreaByPore(crystal, 5).total();
  EXPECT_GT(coarseCrystal, 0.0);
  EXPECT_NEAR(coarseCrystal, fineCrystal, 1.0e-8 * coarseCrystal);
}


// The same surface as the sampled estimate measures, and the same division of it. The area is the one
// thing here the two routes work out independently -- one by integrating the uncovered latitudes, the
// other by counting points that no other atom buries -- so agreeing to the noise of the sampling is a
// check of both. The split is asked of one classifier by both, so it has to agree as well, and does
// not have to be right for that to mean something.
TEST(apollonius_surface_area, agrees_with_the_sampled_estimate)
{
  double a = 12.0;
  UnitCell box(a, a, a);
  const double probe = 0.4;

  RandomNumber random{std::optional<std::size_t>(7)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 18; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(1.2 + 0.8 * random.uniform());
  }

  ApolloniusAccessibility classifier = ApolloniusAccessibility::create(box, fractionalPositions, radii, probe);
  ASSERT_TRUE(classifier.diagram.networkIsComplete());

  MeasuredPatches measured = exactAccessibleSurfaceAreaByPore(classifier.accessibility);
  SurfaceAreaSample sampled = sampleAccessibleSurfaceArea(classifier.accessibility, 400);

  double total = measured.total();
  ASSERT_GT(total, 0.0);
  EXPECT_NEAR(measured.undecided, 0.0, 1.0e-8 * total);

  // 400 points per Å² over some 500 Å² of surface: the standard error of the total is a few tenths of a
  // percent, and of each part of it no worse.
  EXPECT_NEAR(sampled.accessible + sampled.inaccessible, total, 0.01 * total);
  EXPECT_NEAR(sampled.accessible, measured.accessible, 0.01 * total);
  EXPECT_NEAR(sampled.inaccessible, measured.inaccessible, 0.01 * total);
}


// Every point of the boundary of the union belongs to exactly one atom's patch, so the area of the whole
// is the sum over the atoms of theirs, and the classifier is what decides which side of the division
// each patch falls on. That the parts add up to the whole is arithmetic; what this checks is that the
// classifier was asked about all of it and that none of it went missing between the two.
TEST(apollonius_surface_area, the_parts_add_up_to_the_whole)
{
  double a = 13.0;
  UnitCell box(a, a, a);

  RandomNumber random{std::optional<std::size_t>(19)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 22; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(1.3 + 0.9 * random.uniform());
  }

  ApolloniusAccessibility classifier = ApolloniusAccessibility::create(box, fractionalPositions, radii, 0.5);
  MeasuredPatches split = exactAccessibleSurfaceAreaByPore(classifier.accessibility);

  // The bare geometry inflates nothing of its own, so it is handed the inflated radii directly, and then
  // it is the same surface the classifier's is.
  std::vector<double> inflated(radii.size());
  for (std::size_t i = 0; i < radii.size(); ++i) inflated[i] = radii[i] + 0.5;
  MeasuredPatches whole = exactAccessibleSurfaceAreaByPore(bareGeometry(box, fractionalPositions, inflated));

  ASSERT_GT(whole.total(), 0.0);
  EXPECT_NEAR(split.accessible + split.inaccessible + split.undecided, whole.total(), 1.0e-9 * whole.total());
  EXPECT_GT(split.accessible, 0.0);
}


// The area is the union of the atoms and nothing else, so measuring it against the radical network and
// against the Apollonius network has to return the same number to the last bit. What the network decides
// is where that number is divided, and the two do divide it differently: the radical one places the
// threshold too low and calls reachable what a probe cannot get to. Sampled numbers cannot separate the
// two errors, since they differ in the geometry as well; these can, and that is the reason for running the
// measurement against either network.
TEST(apollonius_surface_area, the_total_is_the_same_on_either_network)
{
  double a = 14.0;
  UnitCell box(a, a, a);
  const double probe = 0.6;

  RandomNumber random{std::optional<std::size_t>(101)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 24; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(1.1 + 1.1 * random.uniform());
  }

  ApolloniusAccessibility apollonius = ApolloniusAccessibility::create(box, fractionalPositions, radii, probe);
  PoreAccessibility radical = PoreAccessibility::create(box, fractionalPositions, radii, probe);

  MeasuredPatches fromApollonius = exactAccessibleSurfaceAreaByPore(apollonius.accessibility);
  MeasuredPatches fromRadical = exactAccessibleSurfaceAreaByPore(radical);

  ASSERT_GT(fromApollonius.total(), 0.0);

  // The same geometry measured the same way: the totals are equal to round-off, not to a tolerance.
  EXPECT_NEAR(fromRadical.area, fromApollonius.area, 1.0e-12 * fromApollonius.area);
  EXPECT_NEAR(fromRadical.total(), fromApollonius.total(), 1.0e-12 * fromApollonius.total());
  EXPECT_EQ(fromRadical.diagnostics.numberOfArcs, fromApollonius.diagnostics.numberOfArcs);
}
