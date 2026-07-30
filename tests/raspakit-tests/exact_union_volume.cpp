#include <gtest/gtest.h>

import std;

import double3;
import simulationbox;
import randomnumbers;
import voronoi_network;
import pore_accessibility;
import exact_union_volume;

// The measured volume of the union of the probe-inflated atoms, which is the cell less its void.
//
// The volume is cut into one piece per atom, the part of an atom inside its own cell, and each piece is
// closed: the atom's exposed patch weighted by its radius, and the flat faces of the cell weighted by
// their distance from the centre. Configurations of a few spheres have closed forms and the routine has
// to return them to round-off; a cell packed solid has to come out as the whole cell; and against a
// dense random arrangement, where atoms are buried by several neighbours at once rather than by any one
// of them, it is checked against the sampled volume, to the accuracy the sampling has.
//
// Round-off here means a part in ten thousand million of the volume, not of a cubic angstrom: the flat
// faces are closed form but the patch areas come from a quadrature good to about that, and the volume
// carries it. So the closed forms below are held to a relative tolerance, and a deep overlap, where the
// faces are large and their distances nearly cancel, is expected to use most of it.

namespace
{

PoreAccessibility bareGeometry(const SimulationBox& box, const std::vector<double3>& fractionalPositions,
                                 const std::vector<double>& radii)
{
  VoronoiNetwork network;
  network.simulationBox = box;
  network.atomRadii = radii;
  for (const double3& fractional : fractionalPositions)
  {
    network.atomPositionsFractional.push_back(double3::fract(fractional));
  }
  network.atomNodeVectors.assign(fractionalPositions.size(), {});
  return PoreAccessibility::createFromNetwork(std::move(network), PoreAccessibility::Metric::Clearance);
}


double ballVolume(double radius) { return 4.0 / 3.0 * std::numbers::pi * radius * radius * radius; }


// Volume of the union of two overlapping spheres: the two balls, less the lens they share, which is two
// spherical caps meeting on the circle of intersection.
double unionVolumeOfTwoSpheres(double firstRadius, double secondRadius, double distance)
{
  if (distance >= firstRadius + secondRadius) return ballVolume(firstRadius) + ballVolume(secondRadius);
  if (distance + firstRadius <= secondRadius) return ballVolume(secondRadius);
  if (distance + secondRadius <= firstRadius) return ballVolume(firstRadius);

  double toPlane = (distance * distance + firstRadius * firstRadius - secondRadius * secondRadius) / (2.0 * distance);
  double firstHeight = firstRadius - toPlane;
  double secondHeight = secondRadius - (distance - toPlane);
  double firstCap = std::numbers::pi * firstHeight * firstHeight * (3.0 * firstRadius - firstHeight) / 3.0;
  double secondCap = std::numbers::pi * secondHeight * secondHeight * (3.0 * secondRadius - secondHeight) / 3.0;
  return ballVolume(firstRadius) + ballVolume(secondRadius) - firstCap - secondCap;
}


// The same volume by throwing points into the cell, for the cases that have no closed form.
double sampledUnionVolume(const SimulationBox& box, const std::vector<double3>& fractionalPositions,
                          const std::vector<double>& radii, std::size_t samples, std::size_t seed)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  std::size_t hits = 0;
  for (std::size_t sample = 0; sample < samples; ++sample)
  {
    double3 point = box.cell * double3(random.uniform(), random.uniform(), random.uniform());
    bool inside = false;
    for (std::size_t i = 0; i < radii.size() && !inside; ++i)
    {
      double3 centre = box.cell * fractionalPositions[i];
      double3 delta = point - centre;
      delta = box.cell * double3::fract(box.inverseCell * delta + double3(0.5, 0.5, 0.5)) -
              box.cell * double3(0.5, 0.5, 0.5);
      inside = delta.length() <= radii[i];
    }
    if (inside) ++hits;
  }
  return box.volume * static_cast<double>(hits) / static_cast<double>(samples);
}

}  // namespace


// One sphere, with room around it, is a ball. It is the case in which every atom's cell is the whole box
// and no face contributes, so it tests the sphere half of the sum on its own.
TEST(exact_union_volume, lone_sphere)
{
  double a = 40.0;
  SimulationBox box(a, a, a);
  const double radius = 1.7;

  PoreAccessibility geometry = bareGeometry(box, {double3(0.5, 0.5, 0.5)}, {radius});

  EXPECT_NEAR(unionOfBallsVolume(geometry), ballVolume(radius), 1.0e-9 * ballVolume(radius));
}


// Two spheres against the closed form, at separations that put the plane between them on either side of
// either centre, so that the face enters the sum with both signs.
TEST(exact_union_volume, two_spheres_against_the_closed_form)
{
  double a = 40.0;
  SimulationBox box(a, a, a);

  const std::array<std::array<double, 3>, 6> cases = {{{2.0, 1.0, 2.5},
                                                       {3.0, 2.5, 1.0},
                                                       {1.5, 1.5, 2.9},
                                                       {4.0, 0.6, 3.5},
                                                       {2.2, 2.0, 0.3},
                                                       {3.0, 1.0, 2.2}}};

  for (const std::array<double, 3>& one : cases)
  {
    double firstRadius = one[0];
    double secondRadius = one[1];
    double distance = one[2];

    PoreAccessibility geometry = bareGeometry(
        box, {double3(0.5, 0.5, 0.5), double3(0.5 + distance / a, 0.5, 0.5)}, {firstRadius, secondRadius});

    double expected = unionVolumeOfTwoSpheres(firstRadius, secondRadius, distance);
    EXPECT_NEAR(unionOfBallsVolume(geometry), expected, 1.0e-9 * expected)
        << "radii " << firstRadius << ", " << secondRadius << " at " << distance;
  }
}


// A sphere swallowed by a larger one adds nothing: its cell holds none of it, and the volume is the
// larger ball's. Without that the buried sphere would be counted twice over.
TEST(exact_union_volume, a_buried_sphere_adds_nothing)
{
  double a = 40.0;
  SimulationBox box(a, a, a);
  const double large = 3.0;
  const double small = 0.8;

  PoreAccessibility geometry =
      bareGeometry(box, {double3(0.5, 0.5, 0.5), double3(0.5 + 1.0 / a, 0.5, 0.5)}, {large, small});

  EXPECT_NEAR(unionOfBallsVolume(geometry), ballVolume(large), 1.0e-9 * ballVolume(large));
}


// An atom large enough to cover its own cell in every direction leaves no void, and the union has to
// come out as the cell itself. Every face of the cell is then in play at once, and each is cut back by
// all the others, so it is the one case where the flat half of the sum is checked on its own and exactly.
TEST(exact_union_volume, a_cell_packed_solid_is_the_cell)
{
  const double a = 3.0;
  SimulationBox box(a, a, a);

  // Half the body diagonal is the farthest a corner of the cell can be from the atom at its centre.
  const double radius = 0.5 * std::sqrt(3.0) * a + 0.2;
  PoreAccessibility geometry = bareGeometry(box, {double3(0.5, 0.5, 0.5)}, {radius});

  EXPECT_NEAR(unionOfBallsVolume(geometry), box.volume, 1.0e-9 * box.volume);
}


// Two atoms in a cell they overfill between them, each cut by the other and by the images of both. The
// answer is again the whole cell, which no amount of image bookkeeping may exceed or fall short of.
TEST(exact_union_volume, images_fill_the_cell_too)
{
  const double a = 4.0;
  SimulationBox box(a, a, a);
  const double radius = 0.5 * std::sqrt(3.0) * a;

  PoreAccessibility geometry =
      bareGeometry(box, {double3(0.2, 0.3, 0.7), double3(0.6, 0.8, 0.1)}, {radius, radius});

  EXPECT_NEAR(unionOfBallsVolume(geometry), box.volume, 1.0e-9 * box.volume);
}


// A ball larger than the cell it sits in, cut by its own images on all sides, against the sampled volume.
// The union then has no closed form but it does have a fraction of the cell, and nothing here may claim
// volume the cell does not have.
TEST(exact_union_volume, an_atom_is_cut_by_its_own_images)
{
  const double a = 3.0;
  SimulationBox box(a, a, a);
  const double radius = 1.9;  // more than half the cell, less than half its diagonal

  std::vector<double3> positions = {double3(0.5, 0.5, 0.5)};
  std::vector<double> radii = {radius};
  PoreAccessibility geometry = bareGeometry(box, positions, radii);

  double measured = unionOfBallsVolume(geometry);
  double sampled = sampledUnionVolume(box, positions, radii, 4000000, 17);

  EXPECT_LT(measured, box.volume);
  EXPECT_GT(measured, 0.5 * box.volume);
  EXPECT_NEAR(measured, sampled, 0.01 * box.volume);
}


// Dense random arrangements, where an atom is commonly buried by several neighbours together rather than
// by any one of them, against the sampled volume. This is the case that the closed forms above cannot
// reach and the one in which a decomposition that misses or repeats a piece shows up as a bias.
TEST(exact_union_volume, agrees_with_the_sampled_volume)
{
  const double a = 8.0;
  SimulationBox box(a, a, a);
  RandomNumber random{std::optional<std::size_t>(4321)};

  for (std::size_t trial = 0; trial < 4; ++trial)
  {
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::size_t i = 0; i < 40; ++i)
    {
      positions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
      radii.push_back(0.9 + 1.2 * random.uniform());
    }

    PoreAccessibility geometry = bareGeometry(box, positions, radii);
    double measured = unionOfBallsVolume(geometry);
    double sampled = sampledUnionVolume(box, positions, radii, 2000000, 991 + trial);

    // The sampling's own error at two million points over this cell is a few parts in ten thousand of
    // the cell, so the comparison is held to a little more than that.
    EXPECT_NEAR(measured, sampled, 0.003 * box.volume) << "trial " << trial;
  }
}


// The quadrature behind the patch areas is the only approximation in the sum, and refining it must not
// move the answer. A cage-like arrangement is included because that is where the latitudes at which the
// integrand turns are degenerate and a missed one would show.
TEST(exact_union_volume, refining_the_quadrature_does_not_move_the_answer)
{
  const double a = 8.0;
  SimulationBox box(a, a, a);
  RandomNumber random{std::optional<std::size_t>(2024)};

  std::vector<double3> positions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 30; ++i)
  {
    positions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(1.0 + random.uniform());
  }

  {
    PoreAccessibility geometry = bareGeometry(box, positions, radii);
    double coarse = unionOfBallsVolume(geometry, 1);
    double fine = unionOfBallsVolume(geometry, 4);
    EXPECT_NEAR(coarse, fine, 1.0e-8 * std::abs(fine));
  }

  // Eight atoms on the corners of a cube in the middle of the cell: every sphere meets its neighbours in
  // the same latitudes and three spheres meet in a point.
  {
    std::vector<double3> corners;
    std::vector<double> equal;
    for (std::size_t bits = 0; bits < 8; ++bits)
    {
      corners.push_back(double3(0.4 + 0.2 * static_cast<double>(bits & 1u),
                                0.4 + 0.2 * static_cast<double>((bits >> 1) & 1u),
                                0.4 + 0.2 * static_cast<double>((bits >> 2) & 1u)));
      equal.push_back(1.2);
    }

    PoreAccessibility geometry = bareGeometry(box, corners, equal);
    double coarse = unionOfBallsVolume(geometry, 1);
    double fine = unionOfBallsVolume(geometry, 6);
    EXPECT_NEAR(coarse, fine, 1.0e-8 * std::abs(fine));
  }
}
