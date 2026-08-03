#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import pore_accessibility;
import voronoi_accessible_volume;
import exact_surface_patches;
import exact_union_volume;
import exact_void_split;
import exact_boundary_components;

// The void divided between the channels and the pockets without sampling it.
//
// A pocket is bounded and its boundary is exactly the arcs the surface sweep labels as facing it, so the
// divergence theorem gives its volume from those arcs. Two things have to hold for that to be worth
// anything, and both are checked here rather than assumed: the volume must be the pocket's own volume,
// against an independent measurement of it; and the boundary must have been assembled in one frame, which
// is what a pocket straddling the periodic boundary tests and what the closure defect reports on.

namespace
{

// A sealed cage: `count` spheres with centres on a sphere of radius `shellRadius` about `centre`, each of
// radius `ballRadius`. Taking the ball radius between half the neighbour spacing and the shell radius
// leaves the neighbours overlapping, so the cage is closed, and leaves the middle empty, so there is a
// cavity in it. The centres are the Fibonacci spiral, which spaces them evenly without any symmetry the
// diagram could trip over.
std::vector<double3> cageCentres(const double3& centre, double shellRadius, std::size_t count)
{
  std::vector<double3> positions;
  const double golden = std::numbers::pi * (1.0 + std::sqrt(5.0));
  for (std::size_t i = 0; i < count; ++i)
  {
    double fraction = (static_cast<double>(i) + 0.5) / static_cast<double>(count);
    double cosinePolar = 1.0 - 2.0 * fraction;
    double sinePolar = std::sqrt(std::max(0.0, 1.0 - cosinePolar * cosinePolar));
    double azimuth = golden * (static_cast<double>(i) + 0.5);
    positions.push_back(centre + double3(shellRadius * sinePolar * std::cos(azimuth),
                                         shellRadius * sinePolar * std::sin(azimuth), shellRadius * cosinePolar));
  }
  return positions;
}


PoreAccessibility geometryWithRadii(const UnitCell& box, const std::vector<double3>& cartesianPositions,
                                       const std::vector<double>& radii)
{
  std::vector<double3> fractionalPositions;
  for (const double3& position : cartesianPositions)
  {
    fractionalPositions.push_back(double3::fract(box.inverseCell * position));
  }

  // The radii are already the inflated ones, so the probe adds nothing: what is wanted here is the void
  // of this arrangement and not of a smaller one seen by a probe.
  return PoreAccessibility::create(box, fractionalPositions, radii, 0.0);
}


PoreAccessibility cageGeometry(const UnitCell& box, const std::vector<double3>& cartesianPositions,
                                  double ballRadius)
{
  return geometryWithRadii(box, cartesianPositions, std::vector<double>(cartesianPositions.size(), ballRadius));
}


// The volume of the cavity inside a cage, by throwing points at the cavity itself rather than at the
// cell. The cavity is the part of the ball of radius `shellRadius` about the centre that lies outside
// every sphere of the cage: the shell of centres is covered by the spheres, so nothing inside it connects
// to the outside, and nothing outside it is in the cavity.
double sampledCavityVolume(const std::vector<double3>& cartesianPositions, const double3& centre,
                           double shellRadius, double ballRadius, std::size_t samples, std::size_t seed)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  const double side = 2.0 * shellRadius;
  std::size_t hits = 0;
  for (std::size_t sample = 0; sample < samples; ++sample)
  {
    double3 point = centre + double3(side * (random.uniform() - 0.5), side * (random.uniform() - 0.5),
                                     side * (random.uniform() - 0.5));
    if ((point - centre).length() >= shellRadius) continue;

    bool covered = false;
    for (const double3& position : cartesianPositions)
    {
      if ((point - position).length() < ballRadius)
      {
        covered = true;
        break;
      }
    }
    if (!covered) ++hits;
  }
  return side * side * side * static_cast<double>(hits) / static_cast<double>(samples);
}


// How far the cavity gets from a point, by the same sampling: the largest distance from `from` to a point that
// landed in it. An extreme rather than a mean, so it approaches its limit from below and slowly, which is what
// makes it worth comparing against a measured reach: it can only fail by exceeding it.
double sampledCavityReach(const std::vector<double3>& cartesianPositions, const double3& centre,
                          double shellRadius, double ballRadius, const double3& from, std::size_t samples,
                          std::size_t seed)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  const double side = 2.0 * shellRadius;
  double furthest = 0.0;
  for (std::size_t sample = 0; sample < samples; ++sample)
  {
    double3 point = centre + double3(side * (random.uniform() - 0.5), side * (random.uniform() - 0.5),
                                     side * (random.uniform() - 0.5));
    if ((point - centre).length() >= shellRadius) continue;

    bool covered = false;
    for (const double3& position : cartesianPositions)
    {
      if ((point - position).length() < ballRadius)
      {
        covered = true;
        break;
      }
    }
    if (!covered) furthest = std::max(furthest, (point - from).length());
  }
  return furthest;
}


// Where that cavity is, by the same sampling: the mean of the points that landed in it. Sampling the cavity
// rather than the cell is what makes this worth comparing against, the mean of a few thousand hits settling
// to a hundredth of an Ångström while the same points would give the volume to no better than a percent.
double3 sampledCavityCentre(const std::vector<double3>& cartesianPositions, const double3& centre,
                            double shellRadius, double ballRadius, std::size_t samples, std::size_t seed)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  const double side = 2.0 * shellRadius;
  double3 sum(0.0, 0.0, 0.0);
  std::size_t hits = 0;
  for (std::size_t sample = 0; sample < samples; ++sample)
  {
    double3 point = centre + double3(side * (random.uniform() - 0.5), side * (random.uniform() - 0.5),
                                     side * (random.uniform() - 0.5));
    if ((point - centre).length() >= shellRadius) continue;

    bool covered = false;
    for (const double3& position : cartesianPositions)
    {
      if ((point - position).length() < ballRadius)
      {
        covered = true;
        break;
      }
    }
    if (covered) continue;
    sum += point;
    ++hits;
  }
  return sum * (1.0 / static_cast<double>(std::max<std::size_t>(1, hits)));
}

}  // namespace


// The cavity of one sealed cage, against a direct measurement of that cavity. This is the whole claim:
// the arcs facing a pocket know its volume, with no cell and no connected-component search anywhere.
TEST(exact_void_split, a_sealed_cage_gives_its_cavity)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  PoreAccessibility geometry = cageGeometry(box, positions, ballRadius);

  ExactVoidSplit split = exactVoidSplit(geometry, box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_EQ(split.numberOfPockets, 1uz);
  EXPECT_EQ(split.undecidedArea, 0.0);

  // The boundary closes to round-off, not to a tolerance: the arcs and their moments are closed form, and
  // the only approximation between them and this number is the latitude quadrature. This cage leaves about
  // 8e-10 of the area, its spheres meeting at every angle and none of the latitudes falling together; the
  // frameworks are an order or three better than that, and the split is rejected only above 1e-6.
  EXPECT_LT(split.closureDefect, 1.0e-8);

  double sampled = sampledCavityVolume(positions, centre, shellRadius, ballRadius, 8000000, 31);
  EXPECT_NEAR(split.inaccessibleVolume, sampled, 0.02 * sampled);

  // The void is the union's complement and the pockets are a part of it, so the channels are the rest by
  // subtraction and the three have to close exactly.
  EXPECT_NEAR(split.accessibleVolume + split.inaccessibleVolume, split.voidVolume, 1.0e-9 * split.voidVolume);
}


// The same cage moved onto a corner of the cell, so that its cavity is cut into eight pieces by the
// periodic boundary and each piece is found on atoms in a different part of the home cell. The pocket
// volume may not notice. This is the test that the lattice offsets are right: the arcs of one pocket lie
// against different lifts of it, and summed in the frames they are found in they do not describe a closed
// surface at all.
TEST(exact_void_split, a_pocket_across_the_periodic_boundary_is_the_same_pocket)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  ExactVoidSplit centred =
      exactVoidSplit(cageGeometry(box, cageCentres(double3(0.5 * a, 0.5 * a, 0.5 * a), shellRadius, 32), ballRadius),
                     box.volume);
  ExactVoidSplit corner =
      exactVoidSplit(cageGeometry(box, cageCentres(double3(0.0, 0.0, 0.0), shellRadius, 32), ballRadius), box.volume);

  ASSERT_TRUE(centred.reliable) << centred.rejection;
  ASSERT_TRUE(corner.reliable) << corner.rejection;
  EXPECT_EQ(corner.numberOfPockets, 1uz);
  EXPECT_LT(corner.closureDefect, 1.0e-8);

  // Nothing about the arrangement changed but where it was placed, so the two agree to the accuracy of
  // the quadrature and not to the accuracy of the bookkeeping.
  EXPECT_NEAR(corner.inaccessibleVolume, centred.inaccessibleVolume, 1.0e-6 * centred.inaccessibleVolume);
  EXPECT_NEAR(corner.voidVolume, centred.voidVolume, 1.0e-9 * centred.voidVolume);
}


// A cage in a cell barely large enough to hold it, so that the copies of its cavity are as close to it as
// the cavity is wide. Which copy an arc faces then cannot be settled by asking which node of the pore is
// nearest in some other sense, and a frame taken from the node that happened to classify the arc rather
// than from the copy the arc bounds puts part of the boundary a whole cell away from the rest. Sodalite is
// this case, and it is the one that shows the difference.
TEST(exact_void_split, a_cage_that_nearly_fills_its_cell)
{
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  // The cage reaches shellRadius + ballRadius from its centre, so this leaves a tenth of an ångström
  // between it and its own image.
  const double a = 2.0 * (shellRadius + ballRadius) + 0.1;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  ExactVoidSplit split = exactVoidSplit(cageGeometry(box, positions, ballRadius), box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_EQ(split.numberOfPockets, 1uz);
  EXPECT_LT(split.closureDefect, 1.0e-8);

  double sampled = sampledCavityVolume(positions, centre, shellRadius, ballRadius, 8000000, 53);
  EXPECT_NEAR(split.inaccessibleVolume, sampled, 0.02 * sampled);
}


// Two cages in one cell: two pockets, each with its own frame and its own origin, summed into one total.
// A single origin for both would still be valid in exact arithmetic, but each pocket has to close on its
// own, so this checks that they are kept apart.
TEST(exact_void_split, two_cages_give_two_pockets)
{
  const double a = 20.0;
  UnitCell box(a, a, a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  std::vector<double3> positions = cageCentres(double3(5.0, 5.0, 5.0), shellRadius, 32);
  for (const double3& position : cageCentres(double3(14.0, 14.0, 14.0), shellRadius, 32))
  {
    positions.push_back(position);
  }

  PoreAccessibility geometry = cageGeometry(box, positions, ballRadius);
  ExactVoidSplit split = exactVoidSplit(geometry, box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_EQ(split.numberOfPockets, 2uz);
  EXPECT_LT(split.closureDefect, 1.0e-8);

  ExactVoidSplit one = exactVoidSplit(
      cageGeometry(box, cageCentres(double3(5.0, 5.0, 5.0), shellRadius, 32), ballRadius), box.volume);
  EXPECT_NEAR(split.inaccessibleVolume, 2.0 * one.inaccessibleVolume, 1.0e-6 * split.inaccessibleVolume);
}


// A structure with nothing sealed off has no inaccessible volume at all, and the channels take the whole
// void. Zero has to come out as zero here, not as the small residue of a sum of arcs.
TEST(exact_void_split, an_open_structure_has_no_pockets)
{
  const double a = 12.0;
  UnitCell box(a, a, a);

  PoreAccessibility geometry = cageGeometry(box, {double3(6.0, 6.0, 6.0)}, 2.0);
  ExactVoidSplit split = exactVoidSplit(geometry, box.volume);

  EXPECT_EQ(split.numberOfPockets, 0uz);
  EXPECT_EQ(split.inaccessibleVolume, 0.0);
  EXPECT_NEAR(split.accessibleVolume, split.voidVolume, 1.0e-12 * split.voidVolume);
}


// Refining the latitude quadrature must not move the pocket volume: the arcs and their moments are closed
// form, and the quadrature is the only approximation in either.
TEST(exact_void_split, refining_the_quadrature_does_not_move_the_pocket)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);

  PoreAccessibility geometry = cageGeometry(box, cageCentres(centre, 3.0, 32), 1.8);

  ExactVoidSplit coarse = exactVoidSplit(geometry, box.volume, 1);
  ExactVoidSplit fine = exactVoidSplit(geometry, box.volume, 4);

  ASSERT_TRUE(coarse.reliable) << coarse.rejection;
  ASSERT_TRUE(fine.reliable) << fine.rejection;
  EXPECT_NEAR(coarse.inaccessibleVolume, fine.inaccessibleVolume, 1.0e-6 * fine.inaccessibleVolume);
  EXPECT_LT(fine.closureDefect, coarse.closureDefect * 2.0);
}


// Against the sampled split, which asks the same classifier the same question and so has to agree with it
// to the accuracy of the sampling. This is the comparison that would catch a pocket counted twice or
// missed altogether, which the cases above could not: they know what the answer is, and this one does not.
TEST(exact_void_split, agrees_with_the_sampled_split)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);

  PoreAccessibility geometry = cageGeometry(box, cageCentres(centre, 3.0, 32), 1.8);

  ExactVoidSplit split = exactVoidSplit(geometry, box.volume);
  VolumeSample sampled = sampleAccessibleVolume(geometry, 4000000);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_NEAR(split.inaccessibleVolume / box.volume, sampled.inaccessibleFraction, 0.05 * split.inaccessibleVolume / box.volume);
  EXPECT_NEAR(split.accessibleVolume / box.volume, sampled.accessibleFraction, 0.01 * split.accessibleVolume / box.volume);
}


// The same split taken over the connected surfaces of the boundary instead of over the classifier's pores.
//
// What is being checked is a different thing from the volume, which the tests above already pin down. It is
// that the pieces the divergence theorem is applied to are the right pieces: every arc has to be found on
// the surface it lies on, so that the surface is whole, and each surface has to be asked once which side of
// it the void is on. So these look at the lookup failing nowhere, at the surfaces being asked about at a
// point where the answer can be proved, and at the volume coming out where the other route puts it.

// One sealed cage. Two surfaces: the wall of the cavity and the outside of the cage, of which only the
// first bounds any void.
TEST(exact_void_split, by_components_a_sealed_cage_gives_its_cavity)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  PoreAccessibility geometry = cageGeometry(box, positions, ballRadius);

  BoundaryComponents components = boundaryComponents(geometry);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(geometry, components);
  MeasuredPatches patches =
      exactAccessibleSurfaceAreaByComponent(geometry, components, verdicts, 1, SurfaceMoments::andCentre);

  // Every arc was placed on a patch. This is the whole of what the lookup has to do, and it is exact: an
  // arc ends where a bounding circle crosses the latitude, and the arcs of the circles carry their patch.
  EXPECT_EQ(patches.diagnostics.unplacedArcs, 0uz);
  EXPECT_EQ(patches.diagnostics.unplacedArea, 0.0);

  ExactVoidSplit split = exactVoidSplitByComponents(geometry, components, verdicts, patches, box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_EQ(split.numberOfSurfaces, 2uz);
  EXPECT_EQ(split.numberOfPockets, 1uz);
  EXPECT_EQ(split.undecidedArea, 0.0);
  EXPECT_LT(split.closureDefect, 1.0e-8);

  // Both surfaces were settled by reaching a free ball from them rather than by inferring anything: the
  // cavity's wall reaches the node in the middle of the cavity, the outside reaches a node of the open void.
  EXPECT_EQ(split.provedSurfaces, 2uz);

  // And it is the same cavity the arc-by-arc route finds, which is where the sampled comparison already is.
  // Not to the last digit: the two take their moments about different points, and what survives of that
  // choice is the closure defect times the distance between them, some parts in 1e9 here.
  ExactVoidSplit byPore = exactVoidSplit(geometry, box.volume);
  ASSERT_TRUE(byPore.reliable) << byPore.rejection;
  EXPECT_NEAR(split.inaccessibleVolume, byPore.inaccessibleVolume, 1.0e-7 * byPore.inaccessibleVolume);
  EXPECT_NEAR(split.voidVolume, byPore.voidVolume, 1.0e-12 * byPore.voidVolume);

  // The area divides the same way too, the cavity's wall being the inaccessible part of the surface.
  MeasuredPatches byArc = exactAccessibleSurfaceAreaByPore(geometry);
  EXPECT_NEAR(patches.inaccessible, byArc.inaccessible, 1.0e-9 * byArc.inaccessible);
  EXPECT_NEAR(patches.accessible, byArc.accessible, 1.0e-9 * byArc.accessible);
}


// The cage on a corner of the cell, its cavity cut into eight pieces by the periodic boundary. Here the
// translations come from the decomposition rather than from the classifier, so this is the test that the
// surface is assembled in one frame by construction.
TEST(exact_void_split, by_components_a_pocket_across_the_periodic_boundary_is_the_same_pocket)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  ExactVoidSplit centred = exactVoidSplitByComponents(
      cageGeometry(box, cageCentres(double3(0.5 * a, 0.5 * a, 0.5 * a), shellRadius, 32), ballRadius), box.volume);
  ExactVoidSplit corner = exactVoidSplitByComponents(
      cageGeometry(box, cageCentres(double3(0.0, 0.0, 0.0), shellRadius, 32), ballRadius), box.volume);

  ASSERT_TRUE(centred.reliable) << centred.rejection;
  ASSERT_TRUE(corner.reliable) << corner.rejection;
  EXPECT_EQ(corner.numberOfPockets, 1uz);
  EXPECT_LT(corner.closureDefect, 1.0e-8);
  EXPECT_NEAR(corner.inaccessibleVolume, centred.inaccessibleVolume, 1.0e-6 * centred.inaccessibleVolume);
}


// A cage in a cell barely large enough to hold it, where the copies of the cavity are as close as the cavity
// is wide and no argument from nearness can say which copy an arc bounds.
TEST(exact_void_split, by_components_a_cage_that_nearly_fills_its_cell)
{
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;
  const double a = 2.0 * (shellRadius + ballRadius) + 0.1;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  PoreAccessibility geometry = cageGeometry(box, positions, ballRadius);

  ExactVoidSplit split = exactVoidSplitByComponents(geometry, box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_EQ(split.numberOfPockets, 1uz);
  EXPECT_LT(split.closureDefect, 1.0e-8);

  double sampled = sampledCavityVolume(positions, centre, shellRadius, ballRadius, 8000000, 53);
  EXPECT_NEAR(split.inaccessibleVolume, sampled, 0.02 * sampled);
}


// Two cages, four surfaces, two of which bound a cavity.
TEST(exact_void_split, by_components_two_cages_give_two_pockets)
{
  const double a = 20.0;
  UnitCell box(a, a, a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  std::vector<double3> positions = cageCentres(double3(5.0, 5.0, 5.0), shellRadius, 32);
  for (const double3& position : cageCentres(double3(14.0, 14.0, 14.0), shellRadius, 32))
  {
    positions.push_back(position);
  }

  ExactVoidSplit split = exactVoidSplitByComponents(cageGeometry(box, positions, ballRadius), box.volume);
  ExactVoidSplit one = exactVoidSplitByComponents(
      cageGeometry(box, cageCentres(double3(5.0, 5.0, 5.0), shellRadius, 32), ballRadius), box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  EXPECT_EQ(split.numberOfSurfaces, 4uz);
  EXPECT_EQ(split.numberOfPockets, 2uz);
  EXPECT_LT(split.closureDefect, 1.0e-8);
  EXPECT_NEAR(split.inaccessibleVolume, 2.0 * one.inaccessibleVolume, 1.0e-6 * split.inaccessibleVolume);
}


// A ball floating loose inside the cavity, touching nothing. Its own surface is a third connected surface,
// bounded, and facing the same pocket -- but with the void outside it rather than inside, so what it
// encloses is solid and has to come off the pocket's volume rather than onto it. The divergence theorem
// does this by itself: the normals of that surface point out of the solid, and the sum comes out negative.
// Which is also what tells the two apart, so the ball is counted as the cluster it is and not as a pocket.
TEST(exact_void_split, by_components_a_ball_in_a_cavity_takes_up_room)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;
  const double innerRadius = 0.8;  // 0.8 + 1.8 < 3.0, so it reaches none of the cage

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  std::vector<double> radii(positions.size(), ballRadius);
  positions.push_back(centre);
  radii.push_back(innerRadius);

  ExactVoidSplit withBall = exactVoidSplitByComponents(geometryWithRadii(box, positions, radii), box.volume);
  ExactVoidSplit empty = exactVoidSplitByComponents(
      cageGeometry(box, cageCentres(centre, shellRadius, 32), ballRadius), box.volume);

  ASSERT_TRUE(withBall.reliable) << withBall.rejection;
  ASSERT_TRUE(empty.reliable) << empty.rejection;
  EXPECT_EQ(withBall.numberOfSurfaces, 3uz);
  EXPECT_EQ(withBall.numberOfPockets, 1uz);
  EXPECT_EQ(withBall.numberOfEnclosedSolids, 1uz);
  EXPECT_LT(withBall.closureDefect, 1.0e-8);

  const double ballVolume = 4.0 / 3.0 * std::numbers::pi * innerRadius * innerRadius * innerRadius;
  EXPECT_NEAR(withBall.inaccessibleVolume, empty.inaccessibleVolume - ballVolume,
              1.0e-6 * withBall.inaccessibleVolume);
}


// The area follows the same rule as the volume, surface by surface. Around the ball in the cavity there are
// three surfaces: the outside of the cage, which runs away through the crystal and is reachable; the wall of
// the cavity, which encloses void and so cannot be reached; and the ball, which encloses solid and is
// reachable only if the void around it is, which here it is not. So the whole of the cavity, ball included,
// counts as out of reach, and no area is left undecided.
TEST(exact_void_split, by_components_the_area_of_a_sealed_cluster_is_out_of_reach)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const double innerRadius = 0.8;

  std::vector<double3> positions = cageCentres(centre, 3.0, 32);
  std::vector<double> radii(positions.size(), 1.8);
  positions.push_back(centre);
  radii.push_back(innerRadius);

  PoreAccessibility geometry = geometryWithRadii(box, positions, radii);
  MeasuredPatches patches = exactAccessibleSurfaceAreaByComponent(geometry);

  EXPECT_EQ(patches.diagnostics.unplacedArcs, 0uz);
  EXPECT_EQ(patches.undecided, 0.0);
  EXPECT_NEAR(patches.accessible + patches.inaccessible, patches.area, 1.0e-9 * patches.area);

  // The ball is a sphere on its own, so the area it adds is known in closed form, and it has to land on the
  // unreachable side along with the wall that seals it in.
  MeasuredPatches empty = exactAccessibleSurfaceAreaByComponent(cageGeometry(box, cageCentres(centre, 3.0, 32), 1.8));
  const double ballArea = 4.0 * std::numbers::pi * innerRadius * innerRadius;
  EXPECT_NEAR(patches.inaccessible, empty.inaccessible + ballArea, 1.0e-9 * patches.inaccessible);
  EXPECT_NEAR(patches.accessible, empty.accessible, 1.0e-9 * empty.accessible);
}


// A chain of balls winding round one of them, close enough that consecutive links overlap. What is left of the
// middle sphere is a strip following the chain, and a strip is where the surfaces of the boundary are hardest
// to assemble: a strip cut in two by mistake has its halves on different surfaces, and a surface missing a
// piece does not close. So the closure of every bounded surface is the whole of what is asked here, which is
// also all that has to hold for the volumes to mean anything.
TEST(exact_void_split, by_components_a_winding_surface_still_closes)
{
  UnitCell box(40.0, 40.0, 40.0);
  const double3 centre(20.0, 20.0, 20.0);
  const double radius = 4.0;

  std::vector<double3> positions{centre};
  const std::size_t links = 11;
  for (std::size_t i = 0; i < links; ++i)
  {
    double polar = 0.3 + 2.5 * static_cast<double>(i) / static_cast<double>(links - 1);
    double azimuth = 2.2 * static_cast<double>(i);
    double3 axis(std::sin(polar) * std::cos(azimuth), std::sin(polar) * std::sin(azimuth), std::cos(polar));
    positions.push_back(centre + axis * (radius + 1.4));
  }

  ExactVoidSplit split = exactVoidSplitByComponents(cageGeometry(box, positions, radius), box.volume);

  EXPECT_TRUE(split.reliable) << split.rejection;
  EXPECT_LT(split.closureDefect, 1.0e-8);
}


// Nothing sealed off: one surface, percolating, and no inaccessible volume at all.
TEST(exact_void_split, by_components_an_open_structure_has_no_pockets)
{
  const double a = 12.0;
  UnitCell box(a, a, a);

  ExactVoidSplit split = exactVoidSplitByComponents(cageGeometry(box, {double3(6.0, 6.0, 6.0)}, 2.0), box.volume);

  EXPECT_EQ(split.numberOfSurfaces, 1uz);
  EXPECT_EQ(split.numberOfPockets, 0uz);
  EXPECT_EQ(split.inaccessibleVolume, 0.0);
  EXPECT_NEAR(split.accessibleVolume, split.voidVolume, 1.0e-12 * split.voidVolume);
}


// Refining the latitude quadrature must not move the pocket, and here it must not move which surface an arc
// was found on either: the panels are cut differently, so the arcs are different arcs.
TEST(exact_void_split, by_components_refining_the_quadrature_does_not_move_the_pocket)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);

  PoreAccessibility geometry = cageGeometry(box, cageCentres(centre, 3.0, 32), 1.8);

  ExactVoidSplit coarse = exactVoidSplitByComponents(geometry, box.volume, 1);
  ExactVoidSplit fine = exactVoidSplitByComponents(geometry, box.volume, 4);

  ASSERT_TRUE(coarse.reliable) << coarse.rejection;
  ASSERT_TRUE(fine.reliable) << fine.rejection;
  EXPECT_NEAR(coarse.inaccessibleVolume, fine.inaccessibleVolume, 1.0e-6 * fine.inaccessibleVolume);
}


// Where a pocket is, from the same arcs that say how large it is. The centroid is a moment of the enclosed
// region, so the divergence theorem reaches it as well, and this is it against the mean of points thrown into
// the cavity itself. A ball centred there and no wider than the nearest atom lies wholly in the pocket, which
// is what makes the pair a blocking sphere, so the free radius has to be a real length and not zero.
TEST(exact_void_split, by_components_a_pocket_knows_where_it_is)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  ExactVoidSplit split = exactVoidSplitByComponents(cageGeometry(box, positions, ballRadius), box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  ASSERT_EQ(split.pockets.size(), 1uz);

  double3 sampled = sampledCavityCentre(positions, centre, shellRadius, ballRadius, 4000000, 71);
  EXPECT_NEAR(split.pockets[0].centre.x, sampled.x, 0.02);
  EXPECT_NEAR(split.pockets[0].centre.y, sampled.y, 0.02);
  EXPECT_NEAR(split.pockets[0].centre.z, sampled.z, 0.02);

  EXPECT_NEAR(split.pockets[0].volume, split.inaccessibleVolume, 1.0e-9 * split.inaccessibleVolume);
  EXPECT_GT(split.pockets[0].freeRadius, 0.5);

  // The three radii are nested, and not by convention: the free ball lies inside the pocket and the covering
  // ball holds it, so their volumes bracket the pocket's own and the radii come in that order.
  EXPECT_LE(split.pockets[0].freeRadius, split.pockets[0].equivalentRadius + 1.0e-9);
  EXPECT_LE(split.pockets[0].equivalentRadius, split.pockets[0].coveringRadius + 1.0e-9);
}


// How far the pocket reaches from its own centre, against points thrown into the cavity. The reach is a
// maximum over the surface, taken in closed form on each patch, so what has to hold is that it covers the
// cavity -- no sampled point may lie outside it -- and that it is the least radius that does, which is the
// sampled maximum creeping up to it from below rather than stopping short.
TEST(exact_void_split, by_components_a_pockets_reach_covers_it_and_no_more)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const double shellRadius = 3.0;
  const double ballRadius = 1.8;

  std::vector<double3> positions = cageCentres(centre, shellRadius, 32);
  ExactVoidSplit split = exactVoidSplitByComponents(cageGeometry(box, positions, ballRadius), box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  ASSERT_EQ(split.pockets.size(), 1uz);
  const PocketGeometry& pocket = split.pockets[0];

  double sampled =
      sampledCavityReach(positions, centre, shellRadius, ballRadius, pocket.centre, 4000000, 89);
  EXPECT_LE(sampled, pocket.coveringRadius);
  EXPECT_GT(sampled, pocket.coveringRadius - 0.05);
}


// The same cage, moved. A centroid that came out of a boundary assembled in the wrong frame would not follow
// the pocket across the periodic boundary, so this is the check that the moment is taken in one frame and
// brought home afterwards rather than averaged over the copies of the cell.
TEST(exact_void_split, by_components_a_pockets_centre_moves_with_the_pocket)
{
  const double a = 14.0;
  UnitCell box(a, a, a);
  const double3 shift(0.37 * a, -0.62 * a, 1.24 * a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);

  ExactVoidSplit centred =
      exactVoidSplitByComponents(cageGeometry(box, cageCentres(centre, 3.0, 32), 1.8), box.volume);
  ExactVoidSplit moved =
      exactVoidSplitByComponents(cageGeometry(box, cageCentres(centre + shift, 3.0, 32), 1.8), box.volume);

  ASSERT_EQ(centred.pockets.size(), 1uz);
  ASSERT_EQ(moved.pockets.size(), 1uz);

  double3 expected = double3::fract(box.inverseCell * (centred.pockets[0].centre + shift));
  double3 difference = moved.pockets[0].centreFractional - expected;
  difference -= double3(std::round(difference.x), std::round(difference.y), std::round(difference.z));
  EXPECT_LT((box.cell * difference).length(), 1.0e-6);
}


// Two cages, each with its own centre. A single centroid taken over both pockets at once would land halfway
// between them, in the middle of the framework, which is the mistake this rules out.
TEST(exact_void_split, by_components_two_pockets_have_their_own_centres)
{
  const double a = 20.0;
  UnitCell box(a, a, a);
  const double3 first(5.0, 5.0, 5.0);
  const double3 second(14.0, 13.0, 15.0);

  std::vector<double3> positions = cageCentres(first, 3.0, 32);
  for (const double3& position : cageCentres(second, 3.0, 32)) positions.push_back(position);

  ExactVoidSplit split = exactVoidSplitByComponents(cageGeometry(box, positions, 1.8), box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  ASSERT_EQ(split.pockets.size(), 2uz);

  std::vector<double3> found{split.pockets[0].centre, split.pockets[1].centre};
  if ((found[0] - first).length() > (found[1] - first).length()) std::swap(found[0], found[1]);
  EXPECT_LT((found[0] - first).length(), 0.2);
  EXPECT_LT((found[1] - second).length(), 0.2);
}


// What a blocking sphere has to be, and the one of the two conditions that cannot be checked against the
// pocket alone: the sphere may hold no point a molecule is entitled to sit at. The channel radius is a minimum
// over the surfaces facing the accessible void, taken in closed form, and both halves of it are put to points
// thrown at the cell. No point of the accessible void may fall inside the sphere -- a single one would be a
// pore the simulation had lost -- and the nearest of them has to creep down to the radius rather than stop
// short of it, or the sphere is smaller than it could safely have been.
//
// The two halves cannot use the same points. A cage in open void has the void reaching right up to its outer
// wall, where the free space thins to nothing, and there a point can be in the accessible void without being
// provably so; the proof is what makes the first half a proof, and its absence is what keeps such points from
// coming near enough for the second. So the safe half asks the classifier and the tight half uses what this
// arrangement is known to be: outside both shells, and out of the atoms, is the void around the cages.
TEST(exact_void_split, by_components_a_blocking_sphere_holds_no_accessible_point)
{
  const double a = 20.0;
  const double shellRadius = 3.0;
  UnitCell box(a, a, a);
  const std::vector<double3> centres{double3(5.0, 5.0, 5.0), double3(14.0, 13.0, 15.0)};

  std::vector<double3> positions;
  for (const double3& centre : centres)
  {
    for (const double3& position : cageCentres(centre, shellRadius, 32)) positions.push_back(position);
  }

  PoreAccessibility geometry = cageGeometry(box, positions, 1.8);
  ExactVoidSplit split = exactVoidSplitByComponents(geometry, box.volume);

  ASSERT_TRUE(split.reliable) << split.rejection;
  ASSERT_EQ(split.pockets.size(), 2uz);

  // The cages stand in void that runs away through the crystal, so each of them has a channel to keep clear of,
  // and here it is further off than the cage is wide: one sphere covers the pocket and the cap does not bite.
  for (const PocketGeometry& pocket : split.pockets)
  {
    ASSERT_TRUE(pocket.hasChannel());
    EXPECT_TRUE(pocket.coversPocket());
    EXPECT_EQ(pocket.blockingRadius(), pocket.coveringRadius);
    EXPECT_GT(pocket.channelRadius, pocket.coveringRadius);
  }

  std::vector<double> nearest(split.pockets.size(), std::numeric_limits<double>::max());
  RandomNumber random{std::optional<std::size_t>(2027)};
  for (std::size_t sample = 0; sample < 4000000uz; ++sample)
  {
    double3 point = box.cell * double3(random.uniform(), random.uniform(), random.uniform());

    bool free = geometry.clearance(point) > 0.0;
    bool aroundTheCages =
        free && std::ranges::all_of(centres, [&](const double3& centre)
                                    { return (point - centre).length() > shellRadius; });
    if (!aroundTheCages && !geometry.provablyAccessible(point)) continue;

    for (std::size_t index = 0; index < split.pockets.size(); ++index)
    {
      double distance = (point - split.pockets[index].centre).length();
      EXPECT_GT(distance, split.pockets[index].blockingRadius());
      if (aroundTheCages) nearest[index] = std::min(nearest[index], distance);
    }
  }

  for (std::size_t index = 0; index < split.pockets.size(); ++index)
  {
    EXPECT_GE(nearest[index], split.pockets[index].channelRadius);

    // The nearest point of the void is an extremum of the sampling, so it comes down to its limit slowly and
    // from above; a few hundredths of an angstrom is as close as this many points get to a waist this narrow.
    EXPECT_LT(nearest[index], split.pockets[index].channelRadius + 0.15);
  }
}
