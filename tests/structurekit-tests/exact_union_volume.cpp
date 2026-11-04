#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import pore_accessibility;
import exact_union_volume;
import exact_boundary_components;

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
double sampledUnionVolume(const UnitCell& box, const std::vector<double3>& fractionalPositions,
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
  UnitCell box(a, a, a);
  const double radius = 1.7;

  PoreAccessibility geometry = bareGeometry(box, {double3(0.5, 0.5, 0.5)}, {radius});

  EXPECT_NEAR(unionOfBallsVolume(geometry), ballVolume(radius), 1.0e-9 * ballVolume(radius));
}


// Two spheres against the closed form, at separations that put the plane between them on either side of
// either centre, so that the face enters the sum with both signs.
TEST(exact_union_volume, two_spheres_against_the_closed_form)
{
  double a = 40.0;
  UnitCell box(a, a, a);

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
  UnitCell box(a, a, a);
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
  UnitCell box(a, a, a);

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
  UnitCell box(a, a, a);
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
  UnitCell box(a, a, a);
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
  UnitCell box(a, a, a);
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
  UnitCell box(a, a, a);
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

// Two atoms on the very same spot with the same radius are one ball, not two: the tie is broken on the index, so
// that the volume counts the ball once, as the area does.
TEST(exact_union_volume, coincident_equal_spheres_count_once)
{
  double a = 40.0;
  UnitCell box(a, a, a);
  const double radius = 2.0;
  PoreAccessibility geometry = bareGeometry(box, {double3(0.5, 0.5, 0.5), double3(0.5, 0.5, 0.5)}, {radius, radius});
  EXPECT_NEAR(unionOfBallsVolume(geometry), ballVolume(radius), 1.0e-9 * ballVolume(radius));

  // and a third, on the same spot of an image of the cell, is the same ball again
  PoreAccessibility images =
      bareGeometry(box, {double3(0.5, 0.5, 0.5), double3(1.5, 0.5, 0.5), double3(0.5, 0.5, 0.5)}, {radius, radius, radius});
  EXPECT_NEAR(unionOfBallsVolume(images), ballVolume(radius), 1.0e-9 * ballVolume(radius));
}

// The same tie in the decomposition of the boundary: the one with the higher index is buried, and the surface is
// the one sphere's.
TEST(exact_union_volume, coincident_equal_spheres_leave_one_surface)
{
  double a = 40.0;
  UnitCell box(a, a, a);
  const double radius = 2.0;
  PoreAccessibility geometry =
      bareGeometry(box, {double3(0.5, 0.5, 0.5), double3(0.5, 0.5, 0.5), double3(1.5, 0.5, 0.5)}, {radius, radius, radius});

  BoundaryComponents components = boundaryComponents(geometry);
  ASSERT_EQ(components.atoms.size(), 3uz);
  EXPECT_FALSE(components.atoms[0].buried);
  EXPECT_TRUE(components.atoms[1].buried);
  EXPECT_TRUE(components.atoms[2].buried);
  EXPECT_EQ(components.numberOfComponents, 1uz);
  EXPECT_EQ(components.atoms[0].numberOfPatches, 1uz);
  EXPECT_TRUE(components.atoms[0].circles.empty());
}


// A third sphere cutting the two of a duplicated atom sees one sphere, whichever of the three comes where in the
// input. It used to be handed a cap from each of the pair and to keep the one the tie had buried, which left
// the seam between it and the survivor joined to nothing, and it cut the cell of the third atom twice along the
// one plane.
TEST(exact_union_volume, a_third_sphere_sees_a_coincident_pair_as_one)
{
  const double a = 40.0;
  UnitCell box(a, a, a);
  const double radius = 2.0;
  const double distance = 3.0;
  const double3 twice(0.5, 0.5, 0.5);
  const double3 third(0.5 + distance / a, 0.5, 0.5);

  const double lens = std::numbers::pi * (4.0 * radius + distance) * (2.0 * radius - distance) *
                      (2.0 * radius - distance) / 12.0;
  const double expected = 2.0 * ballVolume(radius) - lens;

  for (const std::vector<double3>& order : {std::vector<double3>{twice, twice, third}, std::vector<double3>{third, twice, twice},
                                            std::vector<double3>{twice, third, twice}})
  {
    PoreAccessibility geometry = bareGeometry(box, order, {radius, radius, radius});
    BoundaryComponents components = boundaryComponents(geometry);

    EXPECT_EQ(components.numberOfComponents, 1uz);
    EXPECT_EQ(components.looseEdges, 0uz);
    EXPECT_NEAR(unionOfBallsVolume(geometry), expected, 1.0e-9 * expected);
  }
}

// The same for a pair with a larger sphere on the same spot, which buries both: only the larger one is a cap.
TEST(exact_union_volume, a_third_sphere_sees_only_the_largest_of_a_coincident_group)
{
  const double a = 40.0;
  UnitCell box(a, a, a);
  const double radius = 2.0;
  const double distance = 3.0;
  const double3 twice(0.5, 0.5, 0.5);
  const double3 third(0.5 + distance / a, 0.5, 0.5);

  const double lens = std::numbers::pi * (4.0 * radius + distance) * (2.0 * radius - distance) *
                      (2.0 * radius - distance) / 12.0;
  const double expected = 2.0 * ballVolume(radius) - lens;

  PoreAccessibility geometry = bareGeometry(box, {twice, twice, twice, third}, {1.0, radius, 1.0, radius});
  BoundaryComponents components = boundaryComponents(geometry);
  EXPECT_EQ(components.numberOfComponents, 1uz);
  EXPECT_EQ(components.looseEdges, 0uz);
  EXPECT_NEAR(unionOfBallsVolume(geometry), expected, 1.0e-9 * expected);
}

// Duplicating atoms of a real arrangement changes nothing that can be measured, wherever in the input the copies
// land. An arrangement dense enough that spheres meet in threes is the case with seams to lose.
TEST(exact_union_volume, duplicating_atoms_of_an_arrangement_changes_nothing)
{
  UnitCell box(14.0, 15.0, 16.0);

  for (std::size_t seed = 0; seed < 6; ++seed)
  {
    RandomNumber random{std::optional<std::size_t>(100 + seed)};
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::size_t i = 0; i < 14; ++i)
    {
      positions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
      radii.push_back(1.6 + random.uniform());
    }

    std::vector<double3> doubledPositions;
    std::vector<double> doubledRadii;
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
      order.push_back(i);
      if (i % 2 == 0) order.push_back(i);
    }
    // A shuffle of the input, so that a copy need not sit next to its original.
    for (std::size_t i = order.size() - 1; i > 0; --i)
    {
      std::swap(order[i], order[static_cast<std::size_t>(random.uniform() * static_cast<double>(i + 1))]);
    }
    for (std::size_t i : order)
    {
      doubledPositions.push_back(positions[i]);
      doubledRadii.push_back(radii[i]);
    }

    PoreAccessibility plain = bareGeometry(box, positions, radii);
    PoreAccessibility doubled = bareGeometry(box, doubledPositions, doubledRadii);

    BoundaryComponents plainComponents = boundaryComponents(plain);
    BoundaryComponents doubledComponents = boundaryComponents(doubled);
    EXPECT_EQ(doubledComponents.looseEdges, plainComponents.looseEdges) << "seed " << seed;
    EXPECT_EQ(doubledComponents.numberOfComponents, plainComponents.numberOfComponents) << "seed " << seed;

    const double expected = unionOfBallsVolume(plain);
    EXPECT_NEAR(unionOfBallsVolume(doubled), expected, 1.0e-8 * expected) << "seed " << seed;
  }
}


namespace
{

// A sphere through the circle two others meet in, with its centre on the line between theirs. All three cut one
// another in that one circle: on either side of its plane the middle one lies inside one of the outer two, so it
// has no surface and takes up no room, although neither of them holds the whole of it.
double middleSphereRadius(const double3& first, double firstRadius, const double3& second, double secondRadius,
                          const double3& middle)
{
  const double separation = (second - first).length();
  const double plane = (separation * separation + firstRadius * firstRadius - secondRadius * secondRadius) /
                       (2.0 * separation);
  const double along = (middle - first).length();
  return std::sqrt(firstRadius * firstRadius - plane * plane + (along - plane) * (along - plane));
}

PoreAccessibility cartesianGeometry(const UnitCell& box, const std::vector<double3>& positions,
                                    const std::vector<double>& radii)
{
  std::vector<double3> fractional;
  for (const double3& position : positions) fractional.push_back(box.inverseCell * position);
  return bareGeometry(box, fractional, radii);
}

}  // namespace

// The closed form: a cap of one ball and a cap of the other, meeting at the plane. Asked along a direction with no
// exact representation, so that nothing about the three is equal to the last bit, and with the middle sphere in
// every place in the input.
TEST(exact_union_volume, a_sphere_through_the_circle_of_two_others_adds_nothing)
{
  UnitCell box(40.0, 40.0, 40.0);
  const double3 direction = double3::normalize(double3(1.0, 0.7, 0.3));
  const double3 first(20.0, 20.0, 20.0);
  const double firstRadius = 2.0;
  const double3 second = first + direction * 4.0;
  const double secondRadius = std::sqrt(8.0);
  const double3 middle = first + direction * 3.0;
  const double middleRadius = middleSphereRadius(first, firstRadius, second, secondRadius, middle);

  const auto cap = [](double radius, double height)
  { return std::numbers::pi * height * height * (3.0 * radius - height) / 3.0; };
  const double plane = 1.5;
  const double expected = ballVolume(firstRadius) - cap(firstRadius, firstRadius - plane) +
                          ballVolume(secondRadius) - cap(secondRadius, secondRadius - (4.0 - plane));

  EXPECT_NEAR(unionOfBallsVolume(cartesianGeometry(box, {first, second}, {firstRadius, secondRadius})), expected,
              1.0e-9 * expected);

  const std::vector<std::vector<std::size_t>> orders = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  const std::array<double3, 3> centres{first, middle, second};
  const std::array<double, 3> radii{firstRadius, middleRadius, secondRadius};
  for (const std::vector<std::size_t>& order : orders)
  {
    std::vector<double3> positions;
    std::vector<double> sizes;
    for (std::size_t i : order)
    {
      positions.push_back(centres[i]);
      sizes.push_back(radii[i]);
    }
    PoreAccessibility geometry = cartesianGeometry(box, positions, sizes);
    BoundaryComponents components = boundaryComponents(geometry);

    EXPECT_EQ(components.numberOfComponents, 1uz);
    EXPECT_EQ(components.looseEdges, 0uz);
    EXPECT_NEAR(unionOfBallsVolume(geometry), expected, 1.0e-9 * expected);
  }
}

// The same inside a real arrangement, where other spheres cut the three and so meet the circle at corners shared
// by all of them, which is where the walk round the edges has the most to get wrong.
TEST(exact_union_volume, a_sphere_through_the_circle_of_two_others_changes_nothing_in_an_arrangement)
{
  UnitCell box(14.0, 15.0, 16.0);
  std::size_t added = 0;

  for (std::size_t seed = 0; seed < 8; ++seed)
  {
    RandomNumber random{std::optional<std::size_t>(200 + seed)};
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::size_t i = 0; i < 12; ++i)
    {
      positions.push_back(box.cell * double3(0.2 + 0.6 * random.uniform(), 0.2 + 0.6 * random.uniform(),
                                             0.2 + 0.6 * random.uniform()));
      radii.push_back(1.6 + random.uniform());
    }

    std::vector<double3> morePositions = positions;
    std::vector<double> moreRadii = radii;
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
      for (std::size_t j = i + 1; j < positions.size(); ++j)
      {
        const double separation = (positions[j] - positions[i]).length();
        const double plane = (separation * separation + radii[i] * radii[i] - radii[j] * radii[j]) / (2.0 * separation);
        if (separation >= radii[i] + radii[j] || std::abs(plane) >= radii[i] || added % 3 != seed % 3)
        {
          if (separation < radii[i] + radii[j]) ++added;
          continue;
        }
        ++added;
        const double3 middle = positions[i] + (positions[j] - positions[i]) * 0.5;
        const std::size_t place = static_cast<std::size_t>(random.uniform() * static_cast<double>(morePositions.size() + 1));
        morePositions.insert(morePositions.begin() + static_cast<std::ptrdiff_t>(place), middle);
        moreRadii.insert(moreRadii.begin() + static_cast<std::ptrdiff_t>(place),
                         middleSphereRadius(positions[i], radii[i], positions[j], radii[j], middle));
      }
    }
    ASSERT_GT(morePositions.size(), positions.size());

    PoreAccessibility plain = cartesianGeometry(box, positions, radii);
    PoreAccessibility more = cartesianGeometry(box, morePositions, moreRadii);
    BoundaryComponents plainComponents = boundaryComponents(plain);
    BoundaryComponents moreComponents = boundaryComponents(more);

    EXPECT_EQ(moreComponents.numberOfComponents, plainComponents.numberOfComponents) << "seed " << seed;
    EXPECT_EQ(moreComponents.looseEdges, plainComponents.looseEdges) << "seed " << seed;
    const double expected = unionOfBallsVolume(plain);
    EXPECT_NEAR(unionOfBallsVolume(more), expected, 1.0e-8 * expected) << "seed " << seed;
  }
}
