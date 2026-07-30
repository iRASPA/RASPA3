#include <gtest/gtest.h>

import std;

import int3;
import double3;
import simulationbox;
import pore_accessibility;
import exact_boundary_components;

// The boundary of the union of the balls, cut into its connected surfaces.
//
// Every case here has an answer known before it is run, from the arrangement rather than from the code: a
// sealed cage has an inside and an outside and nothing else, so its boundary is two surfaces and every one
// of its atoms carries one patch of each; a chain of balls round the periodic boundary is one surface that
// closes on a translate of itself; and moving any of it about the cell changes none of that.
//
// The last group of tests goes underneath those and asks the two predicates the cutting is built on directly.
// `exposedGreatCircleArc` is the closed form for whether the shortest path between two directions stays clear
// of every cap, and it is the only geometry in the whole of the patch merging; `connectedOnSphere` is that
// predicate extended to paths of several legs. A true answer from either is meant to be a proof, so what
// matters is not only that they say yes where they should but that they never say yes where they should not:
// an arc wrongly called exposed merges two patches that are not one, and no later step can undo it.

namespace
{

// Centres on a sphere, spaced by the Fibonacci spiral so that no symmetry of the arrangement can flatter
// the geometry. With the ball radius above half the neighbour spacing the cage is sealed, and below the
// shell radius it leaves a cavity.
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


PoreAccessibility geometryOf(const SimulationBox& box, const std::vector<double3>& positions, double ballRadius)
{
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const double3& position : positions)
  {
    fractionalPositions.push_back(double3::fract(box.inverseCell * position));
    radii.push_back(ballRadius);
  }
  return PoreAccessibility::create(box, fractionalPositions, radii, 0.0);
}


// The patches of one sphere counted without any of the machinery under test: a fine grid of directions,
// those outside every other ball kept, and the kept ones joined when they are neighbours on the grid. It is
// a raster, so it can miss a connection thinner than its spacing, which is why it is only ever asked about
// arrangements whose patches are broad.
std::vector<std::int32_t> lastFloodFill;
std::size_t lastFloodFillLatitudes = 0;

std::size_t patchesByFloodFill(const PoreAccessibility& accessibility, std::size_t atomIndex,
                               std::size_t latitudes = 400)
{
  const double radius = accessibility.atomRadii[atomIndex];
  const double3 centre = accessibility.atomPositions[atomIndex];
  const std::size_t longitudes = 2 * latitudes;

  std::vector<std::pair<double3, double>> neighbours;
  for (const auto& entry : accessibility.neighbourAtoms(centre, radius + accessibility.maximumAtomRadius))
  {
    if (entry.first.length() < 1.0e-12) continue;
    neighbours.push_back(entry);
  }

  auto exposedAt = [&](std::size_t i, std::size_t j)
  {
    double polar = std::numbers::pi * (static_cast<double>(i) + 0.5) / static_cast<double>(latitudes);
    double azimuth = 2.0 * std::numbers::pi * static_cast<double>(j) / static_cast<double>(longitudes);
    double3 point = centre + double3(radius * std::sin(polar) * std::cos(azimuth),
                                     radius * std::sin(polar) * std::sin(azimuth), radius * std::cos(polar));
    for (const auto& [delta, neighbourRadius] : neighbours)
    {
      if ((point - (centre + delta)).length() < neighbourRadius) return false;
    }
    return true;
  };

  std::vector<std::uint8_t> exposed(latitudes * longitudes, 0);
  for (std::size_t i = 0; i < latitudes; ++i)
  {
    for (std::size_t j = 0; j < longitudes; ++j) exposed[i * longitudes + j] = exposedAt(i, j) ? 1 : 0;
  }

  std::vector<std::int32_t> label(latitudes * longitudes, -1);
  std::size_t components = 0;
  for (std::size_t seed = 0; seed < exposed.size(); ++seed)
  {
    if (exposed[seed] == 0 || label[seed] >= 0) continue;
    std::vector<std::size_t> queue{seed};
    label[seed] = static_cast<std::int32_t>(components);
    while (!queue.empty())
    {
      std::size_t cell = queue.back();
      queue.pop_back();
      std::size_t i = cell / longitudes;
      std::size_t j = cell % longitudes;
      std::array<std::pair<std::size_t, std::size_t>, 4> around = {
          std::pair{i, (j + 1) % longitudes}, std::pair{i, (j + longitudes - 1) % longitudes},
          std::pair{(i == 0) ? 0 : i - 1, j}, std::pair{std::min(latitudes - 1, i + 1), j}};
      for (const auto& [ni, nj] : around)
      {
        std::size_t next = ni * longitudes + nj;
        if (next == cell || exposed[next] == 0 || label[next] >= 0) continue;
        label[next] = static_cast<std::int32_t>(components);
        queue.push_back(next);
      }
    }
    ++components;
  }
  lastFloodFill = label;
  lastFloodFillLatitudes = latitudes;
  return components;
}


// Which region of the last flood fill a direction fell in, or -1 where the raster has it covered.
std::int32_t floodFillRegionOf(const double3& unitDirection)
{
  const std::size_t latitudes = lastFloodFillLatitudes;
  const std::size_t longitudes = 2 * latitudes;
  double polar = std::acos(std::clamp(unitDirection.z, -1.0, 1.0));
  double azimuth = std::atan2(unitDirection.y, unitDirection.x);
  if (azimuth < 0.0) azimuth += 2.0 * std::numbers::pi;

  std::size_t i = std::min(latitudes - 1, static_cast<std::size_t>(polar / std::numbers::pi *
                                                                  static_cast<double>(latitudes)));
  std::size_t j = static_cast<std::size_t>(azimuth / (2.0 * std::numbers::pi) *
                                           static_cast<double>(longitudes)) % longitudes;
  return lastFloodFill[i * longitudes + j];
}


// A cap of this half angle about this direction, for the predicates below. Only the axis and the cosine
// are read by them, but the rest is filled so that the circles are the ones the decomposition would build.
SphereCircle capOf(const double3& axis, double halfAngle)
{
  SphereCircle circle;
  circle.axis = double3::normalize(axis);
  circle.halfAngle = halfAngle;
  circle.cosineHalfAngle = std::cos(halfAngle);
  circle.sineHalfAngle = std::sin(halfAngle);
  return circle;
}

double3 fromPolar(double polar, double azimuth)
{
  return double3(std::sin(polar) * std::cos(azimuth), std::sin(polar) * std::sin(azimuth), std::cos(polar));
}

}  // namespace


// Two balls too far apart to touch: two surfaces, each of one patch, neither closing on a translate of
// itself. The simplest thing the decomposition can be asked, and the one that fails first if the joining
// across circles ever invents an edge.
TEST(boundary_components, separate_balls_are_separate_surfaces)
{
  SimulationBox box(30.0, 30.0, 30.0);
  PoreAccessibility geometry = geometryOf(box, {double3(8.0, 8.0, 8.0), double3(20.0, 20.0, 20.0)}, 2.0);

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(components.numberOfComponents, 2uz);
  EXPECT_EQ(components.numberOfPatches, 2uz);
  EXPECT_EQ(components.atoms[0].numberOfPatches, 1uz);
  EXPECT_EQ(components.atoms[1].numberOfPatches, 1uz);
  EXPECT_EQ(components.componentPercolates[0], 0);
  EXPECT_EQ(components.componentPercolates[1], 0);
  EXPECT_EQ(components.looseEdges, 0uz);
}


// Two balls that overlap: one surface, and one patch on each of them, joined across the circle they share.
TEST(boundary_components, overlapping_balls_are_one_surface)
{
  SimulationBox box(30.0, 30.0, 30.0);
  PoreAccessibility geometry = geometryOf(box, {double3(14.0, 15.0, 15.0), double3(16.0, 15.0, 15.0)}, 2.0);

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(components.numberOfComponents, 1uz);
  EXPECT_EQ(components.numberOfPatches, 2uz);
  EXPECT_EQ(components.componentPercolates[0], 0);
  EXPECT_EQ(components.looseEdges, 0uz);
}


// A ball with two others cutting caps from opposite sides. What is left of it is a band, connected, and
// bounded by two separate loops -- the case the walk around the edges cannot settle on its own, since
// neither loop meets the other. It is one patch all the same, and the flood fill says so independently.
TEST(boundary_components, a_band_between_two_caps_is_one_patch)
{
  SimulationBox box(30.0, 30.0, 30.0);
  PoreAccessibility geometry =
      geometryOf(box, {double3(15.0, 15.0, 15.0), double3(12.5, 15.0, 15.0), double3(17.5, 15.0, 15.0)}, 2.0);

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(patchesByFloodFill(geometry, 0), 1uz);
  EXPECT_EQ(components.atoms[0].numberOfPatches, 1uz);
  EXPECT_EQ(components.numberOfComponents, 1uz);
  EXPECT_EQ(components.looseEdges, 0uz);
}


// A chain of balls winding round a sphere in a helix, overlapping their neighbours in the chain. What is left
// of the sphere is a strip that follows the chain round, so a point of it is joined to another only by a path
// that turns as often as the strip does: no single arc of a great circle lies along it, and no path of two
// arcs does either. Whether the region is followed as far as it goes decides how many patches the sphere is
// cut into, and a patch cut in two puts its halves on different surfaces of the boundary, after which neither
// of those surfaces closes.
//
// What is checked is that no two patches are one region wrongly cut apart --- each has to be shown connected
// to be joined, so a wrong join would be a wrong proof --- and the raster settles that without its resolution
// mattering, two representatives in one region being the thing that must not happen.
TEST(boundary_components, a_region_winding_round_the_sphere_is_not_cut_up)
{
  SimulationBox box(40.0, 40.0, 40.0);
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
  PoreAccessibility geometry = geometryOf(box, positions, radius);

  BoundaryComponents components = boundaryComponents(geometry);
  patchesByFloodFill(geometry, 0, 1600);

  EXPECT_EQ(components.looseEdges, 0uz);
  std::vector<std::int32_t> seen;
  for (const double3& representative : components.atoms[0].patchRepresentative)
  {
    std::int32_t region = floodFillRegionOf(representative);
    ASSERT_GE(region, 0);
    EXPECT_EQ(std::ranges::find(seen, region), seen.end()) << "two patches in region " << region;
    seen.push_back(region);
  }
}


// A sealed cage. Its boundary is two surfaces, the outside and the wall of the cavity, and every atom of it
// carries exactly one patch of each: this is the case the whole exercise is for, since an atom with a patch
// on each side of a wall is what no per-arc classification can be trusted with.
TEST(boundary_components, a_sealed_cage_has_an_inside_and_an_outside)
{
  const double a = 14.0;
  SimulationBox box(a, a, a);
  const std::size_t count = 32;
  std::vector<double3> positions = cageCentres(double3(0.5 * a, 0.5 * a, 0.5 * a), 3.0, count);
  PoreAccessibility geometry = geometryOf(box, positions, 1.8);

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(components.numberOfComponents, 2uz);
  EXPECT_EQ(components.numberOfPatches, 2 * count);
  EXPECT_EQ(components.looseEdges, 0uz);
  EXPECT_EQ(components.unjoinedPatches, 0uz);
  for (std::size_t i = 0; i < count; ++i)
  {
    EXPECT_EQ(components.atoms[i].numberOfPatches, 2uz) << "atom " << i;
    EXPECT_EQ(patchesByFloodFill(geometry, i), 2uz) << "atom " << i;
  }

  // Neither surface runs anywhere: the cage is one object in the middle of the cell.
  EXPECT_EQ(components.componentPercolates[0], 0);
  EXPECT_EQ(components.componentPercolates[1], 0);

  // The two surfaces are told apart by the room in front of them, and the one asked about the cavity is an
  // atom of the cage seen from the inside.
  EXPECT_NE(components.componentRepresentative[0].first, components.componentRepresentative[1].first * 0 - 1);
}


// Looking up which patch a direction is on, on an atom that has one patch facing a cavity and another facing
// the open. Both have to be found, and they have to come out different: an atom of a cage wall is where a
// lookup that guessed would put the two sides of the wall on the same side of it.
//
// The two directions are chosen to be the awkward ones rather than the easy ones. Radially outward from the
// cage, the sphere is exposed for a long way round and no neighbour's circle is anywhere near, which is the
// case where there is no edge to read the patch off and the answer has to be walked to.
TEST(boundary_components, a_direction_is_looked_up_on_the_patch_it_lies_on)
{
  const double a = 14.0;
  SimulationBox box(a, a, a);
  const double3 centre(0.5 * a, 0.5 * a, 0.5 * a);
  const std::size_t count = 32;
  std::vector<double3> positions = cageCentres(centre, 3.0, count);
  PoreAccessibility geometry = geometryOf(box, positions, 1.8);

  BoundaryComponents components = boundaryComponents(geometry);
  ASSERT_EQ(components.numberOfComponents, 2uz);

  for (std::size_t i = 0; i < count; ++i)
  {
    double3 outward = positions[i] - centre;
    outward = outward * (1.0 / outward.length());

    std::int32_t away = components.patchOfDirection(i, outward);
    std::int32_t toward = components.patchOfDirection(i, outward * -1.0);
    ASSERT_GE(away, 0) << "atom " << i;
    ASSERT_GE(toward, 0) << "atom " << i;
    EXPECT_NE(away, toward) << "atom " << i;

    // And the two are on the two surfaces, the cavity's wall and the outside of the cage.
    EXPECT_NE(components.componentOfPatch[i][static_cast<std::size_t>(away)],
              components.componentOfPatch[i][static_cast<std::size_t>(toward)])
        << "atom " << i;
  }

  // The whole cage seen from outside is one surface, so the outward direction of every atom must land on the
  // same one of the two.
  std::int32_t outsideLabel = -1;
  for (std::size_t i = 0; i < count; ++i)
  {
    double3 outward = positions[i] - centre;
    outward = outward * (1.0 / outward.length());
    std::int32_t patch = components.patchOfDirection(i, outward);
    std::int32_t label = components.componentOfPatch[i][static_cast<std::size_t>(patch)];
    if (outsideLabel < 0) outsideLabel = label;
    EXPECT_EQ(label, outsideLabel) << "atom " << i;
  }
}


// The same cage on a corner of the cell, so that its patches are found on atoms scattered through the home
// cell and joined only through the images. Nothing about the answer may change.
TEST(boundary_components, a_cage_on_the_boundary_is_the_same_cage)
{
  const double a = 14.0;
  SimulationBox box(a, a, a);
  const std::size_t count = 32;
  PoreAccessibility geometry = geometryOf(box, cageCentres(double3(0.0, 0.0, 0.0), 3.0, count), 1.8);

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(components.numberOfComponents, 2uz);
  EXPECT_EQ(components.numberOfPatches, 2 * count);
  EXPECT_EQ(components.looseEdges, 0uz);
  EXPECT_EQ(components.componentPercolates[0], 0);
  EXPECT_EQ(components.componentPercolates[1], 0);
}


// Two cages: four surfaces, and the two cavities are not one another.
TEST(boundary_components, two_cages_give_four_surfaces)
{
  const double a = 24.0;
  SimulationBox box(a, a, a);
  std::vector<double3> positions = cageCentres(double3(6.0, 6.0, 6.0), 3.0, 32);
  for (const double3& position : cageCentres(double3(18.0, 18.0, 18.0), 3.0, 32)) positions.push_back(position);
  PoreAccessibility geometry = geometryOf(box, positions, 1.8);

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(components.numberOfComponents, 4uz);
  EXPECT_EQ(components.looseEdges, 0uz);
  for (std::size_t component = 0; component < 4; ++component) EXPECT_EQ(components.componentPercolates[component], 0);
}


// A chain of overlapping balls round the periodic boundary. It is one surface, and following it brings one
// back to where one started a cell away, which is what percolation is and what a bounded surface cannot do.
// Nothing but the accumulated translations can tell this from a ring closing on itself.
TEST(boundary_components, a_chain_round_the_cell_closes_on_a_translate)
{
  const double a = 20.0;
  SimulationBox box(a, a, a);
  std::vector<double3> positions;
  const std::size_t count = 10;
  for (std::size_t i = 0; i < count; ++i)
  {
    positions.push_back(double3(a * static_cast<double>(i) / static_cast<double>(count), 10.0, 10.0));
  }
  PoreAccessibility geometry = geometryOf(box, positions, 1.5);  // spacing 2.0, so consecutive balls overlap

  BoundaryComponents components = boundaryComponents(geometry);

  EXPECT_EQ(components.numberOfComponents, 1uz);
  EXPECT_EQ(components.looseEdges, 0uz);
  EXPECT_EQ(components.componentPercolates[0], 1);
}


// The nesting rule on the case the walk round the edges cannot settle: a band bounded by two loops that never
// meet. Each loop lies inside the other -- from a point of the one, the nearest of the other's circles is
// reached on the exposed side of it -- so the band is one patch, decided rather than searched for.
TEST(boundary_components, the_nesting_rule_makes_the_band_one_patch)
{
  SimulationBox box(30.0, 30.0, 30.0);
  PoreAccessibility geometry =
      geometryOf(box, {double3(15.0, 15.0, 15.0), double3(12.5, 15.0, 15.0), double3(17.5, 15.0, 15.0)}, 2.0);

  BoundaryComponents components = boundaryComponents(geometry, LoopMerge::nesting);

  EXPECT_EQ(patchesByFloodFill(geometry, 0), 1uz);
  EXPECT_EQ(components.atoms[0].numberOfPatches, 1uz);
  EXPECT_EQ(components.numberOfComponents, 1uz);
  EXPECT_EQ(components.looseEdges, 0uz);
}


// The nesting rule where the paths one had to be taught to keep going: a strip winding round the sphere behind
// a helical chain. As before, no two patches may fall in one region of the raster, since each join has to be
// earned; and the count may not come out above what the paths rule reaches, a rule that decides being of no
// use if it decides less often than a search.
TEST(boundary_components, the_nesting_rule_follows_a_winding_region)
{
  SimulationBox box(40.0, 40.0, 40.0);
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
  PoreAccessibility geometry = geometryOf(box, positions, radius);

  BoundaryComponents byNesting = boundaryComponents(geometry, LoopMerge::nesting);
  BoundaryComponents byPaths = boundaryComponents(geometry, LoopMerge::paths);
  patchesByFloodFill(geometry, 0, 1600);

  EXPECT_EQ(byNesting.looseEdges, 0uz);
  EXPECT_LE(byNesting.atoms[0].numberOfPatches, byPaths.atoms[0].numberOfPatches);

  std::vector<std::int32_t> seen;
  for (const double3& representative : byNesting.atoms[0].patchRepresentative)
  {
    std::int32_t region = floodFillRegionOf(representative);
    ASSERT_GE(region, 0);
    EXPECT_EQ(std::ranges::find(seen, region), seen.end()) << "two patches in region " << region;
    seen.push_back(region);
  }
}


// The two rules on the arrangements above, side by side. They answer the same question by different means --
// one looks for a path across the region, the other asks which side of a loop a point is on -- so where both
// are right they cut the spheres the same way, and a cage is where being wrong shows: merging the two sides of
// a wall would leave one surface where there are two.
TEST(boundary_components, the_two_merge_rules_agree)
{
  const double a = 14.0;
  SimulationBox box(a, a, a);
  const std::size_t count = 32;

  std::vector<std::pair<std::string, PoreAccessibility>> cases;
  cases.emplace_back("sealed cage", geometryOf(box, cageCentres(double3(0.5 * a, 0.5 * a, 0.5 * a), 3.0, count), 1.8));
  cases.emplace_back("cage on the corner", geometryOf(box, cageCentres(double3(0.0, 0.0, 0.0), 3.0, count), 1.8));

  std::vector<double3> chain;
  for (std::size_t i = 0; i < 10; ++i) chain.push_back(double3(a * static_cast<double>(i) / 10.0, 7.0, 7.0));
  cases.emplace_back("chain round the cell", geometryOf(box, chain, 1.5));

  for (const auto& [name, geometry] : cases)
  {
    BoundaryComponents byPaths = boundaryComponents(geometry, LoopMerge::paths);
    BoundaryComponents byNesting = boundaryComponents(geometry, LoopMerge::nesting);

    EXPECT_EQ(byNesting.numberOfPatches, byPaths.numberOfPatches) << name;
    EXPECT_EQ(byNesting.numberOfComponents, byPaths.numberOfComponents) << name;
    EXPECT_EQ(byNesting.looseEdges, byPaths.looseEdges) << name;
    EXPECT_EQ(byNesting.unjoinedPatches, byPaths.unjoinedPatches) << name;
    for (std::size_t component = 0; component < byNesting.numberOfComponents; ++component)
    {
      EXPECT_EQ(byNesting.componentPercolates[component], byPaths.componentPercolates[component])
          << name << ", surface " << component;
    }
  }
}


// The images a query reports have to be the images they claim to be, and the two atoms sharing a circle have
// to see it in opposite copies of the cell, or nothing joined across a circle is joined to the right thing.
TEST(boundary_components, neighbour_images_are_consistent_and_symmetric)
{
  const double a = 9.0;  // small enough that the images of an atom reach its own sphere
  SimulationBox box(a, a, a);
  PoreAccessibility geometry =
      geometryOf(box, {double3(1.0, 2.0, 3.0), double3(5.0, 6.0, 7.0), double3(8.5, 0.5, 4.0)}, 2.2);

  for (std::size_t i = 0; i < geometry.atomPositions.size(); ++i)
  {
    double3 centre = geometry.atomPositions[i];
    for (const NeighbourImage& neighbour : geometry.neighbourAtomImages(centre, 6.0))
    {
      // The image is where it says it is.
      double3 claimed = geometry.atomPositions[neighbour.index] +
                        box.cell * double3(static_cast<double>(neighbour.image.x),
                                           static_cast<double>(neighbour.image.y),
                                           static_cast<double>(neighbour.image.z));
      EXPECT_NEAR((claimed - (centre + neighbour.delta)).length(), 0.0, 1.0e-9);

      // And the neighbour sees this atom in the opposite copy, at the same distance.
      bool mirrored = false;
      for (const NeighbourImage& back :
           geometry.neighbourAtomImages(geometry.atomPositions[neighbour.index], 6.0))
      {
        if (back.index == i && back.image.x == -neighbour.image.x && back.image.y == -neighbour.image.y &&
            back.image.z == -neighbour.image.z)
        {
          mirrored = true;
          EXPECT_NEAR(back.delta.length(), neighbour.delta.length(), 1.0e-9);
        }
      }
      if (neighbour.delta.length() > 1.0e-12) EXPECT_TRUE(mirrored);
    }
  }
}


// How many independent directions a surface runs away in, on arrangements built to run away in a known number
// of them. The rank is the one over the translations the surface closes on itself by, so it is decided in
// integers and each of these has an answer that can be read off the arrangement: a cage closes on itself and
// on nothing else, a chain of balls round the cell closes on itself one cell along and nowhere else, a sealed
// sheet closes on itself in the two directions it spreads in, and a scaffold of rods along all three axes
// closes on itself in all three.
TEST(boundary_components, a_surface_knows_how_many_directions_it_runs_away_in)
{
  const double a = 16.0;
  SimulationBox box(a, a, a);

  // A cage: two surfaces, an inside and an outside, neither going anywhere.
  {
    BoundaryComponents cage = boundaryComponents(geometryOf(box, cageCentres(double3(8.0, 8.0, 8.0), 3.0, 32), 1.8));
    ASSERT_EQ(cage.numberOfComponents, 2uz);
    EXPECT_EQ(cage.componentDimensionality[0], 0);
    EXPECT_EQ(cage.componentDimensionality[1], 0);
  }

  // A chain of overlapping balls round the cell along x. One surface, closing on itself one cell along x.
  {
    std::vector<double3> positions;
    for (std::size_t i = 0; i < 8; ++i) positions.push_back(double3(2.0 * static_cast<double>(i), 8.0, 8.0));

    BoundaryComponents chain = boundaryComponents(geometryOf(box, positions, 1.5));
    ASSERT_EQ(chain.numberOfComponents, 1uz);
    EXPECT_EQ(chain.componentDimensionality[0], 1);

    // And the direction is the one it was built along, not merely some direction.
    ASSERT_FALSE(chain.componentTranslations[0].empty());
    for (const int3& translation : chain.componentTranslations[0])
    {
      EXPECT_NE(translation.x, 0);
      EXPECT_EQ(translation.y, 0);
      EXPECT_EQ(translation.z, 0);
    }
  }

  // A sheet of balls, close enough together to leave no way through it. Two surfaces, one above and one
  // below, each spreading in x and y and neither reaching the other.
  {
    std::vector<double3> positions;
    for (std::size_t i = 0; i < 8; ++i)
    {
      for (std::size_t j = 0; j < 8; ++j)
      {
        positions.push_back(double3(2.0 * static_cast<double>(i), 2.0 * static_cast<double>(j), 8.0));
      }
    }

    BoundaryComponents sheet = boundaryComponents(geometryOf(box, positions, 1.5));
    ASSERT_EQ(sheet.numberOfComponents, 2uz);
    EXPECT_EQ(sheet.componentDimensionality[0], 2);
    EXPECT_EQ(sheet.componentDimensionality[1], 2);
  }

  // A scaffold of rods along all three axes, joined at the corners of a coarser lattice. One surface, and it
  // comes back to itself along each of the three.
  {
    std::vector<double3> positions;
    for (std::size_t i = 0; i < 2; ++i)
    {
      for (std::size_t j = 0; j < 2; ++j)
      {
        for (std::size_t k = 0; k < 2; ++k)
        {
          double3 node(8.0 * static_cast<double>(i), 8.0 * static_cast<double>(j), 8.0 * static_cast<double>(k));
          positions.push_back(node);
          for (double along : {2.0, 4.0, 6.0})
          {
            positions.push_back(node + double3(along, 0.0, 0.0));
            positions.push_back(node + double3(0.0, along, 0.0));
            positions.push_back(node + double3(0.0, 0.0, along));
          }
        }
      }
    }

    BoundaryComponents scaffold = boundaryComponents(geometryOf(box, positions, 1.2));
    ASSERT_EQ(scaffold.numberOfComponents, 1uz);
    EXPECT_EQ(scaffold.componentDimensionality[0], 3);
  }
}


// Where the rank of a surface is not the dimensionality of the pore behind it, which is the case the reading
// has to be honest about. A single chain of balls threaded through the cell is a rod, and its surface closes
// on itself along the rod and nowhere else, so the surface is one-dimensional; but the void it stands in is
// the whole of the rest of the cell and runs away in all three directions. The rank is a lower bound on the
// pore's, and here it is a strict one, because the pore is walled by every copy of this surface at once and
// no surface can say that of itself. Only a pore network joins them up.
TEST(boundary_components, a_rod_in_open_void_bounds_more_than_it_says)
{
  const double a = 16.0;
  SimulationBox box(a, a, a);

  std::vector<double3> positions;
  for (std::size_t i = 0; i < 8; ++i) positions.push_back(double3(2.0 * static_cast<double>(i), 8.0, 8.0));

  BoundaryComponents rod = boundaryComponents(geometryOf(box, positions, 1.5));

  ASSERT_EQ(rod.numberOfComponents, 1uz);
  EXPECT_EQ(rod.componentDimensionality[0], 1);
  EXPECT_EQ(rod.componentPercolates[0], 1);
}

TEST(boundary_components, a_sphere_with_no_caps_joins_everything)
{
  const std::vector<SphereCircle> none;
  EXPECT_TRUE(exposedGreatCircleArc(none, double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0)));
  EXPECT_TRUE(exposedGreatCircleArc(none, double3(1.0, 0.0, 0.0), double3(0.0, 0.0, 1.0)));
}

TEST(boundary_components, a_direction_is_joined_to_itself)
{
  const std::vector<SphereCircle> circles = {capOf(double3(0.0, 0.0, 1.0), 0.5)};
  const double3 point = fromPolar(2.0, 0.3);
  EXPECT_TRUE(exposedGreatCircleArc(circles, point, point));

  // And so is a direction inside a cap, the arc between the two being of no length at all. The predicate is
  // asked only about points of the exposed region, so this is a statement about the degenerate case and not
  // about coverage.
  EXPECT_TRUE(exposedGreatCircleArc(circles, double3(0.0, 0.0, 1.0), double3(0.0, 0.0, 1.0)));
}

TEST(boundary_components, two_antipodal_directions_are_joined_by_no_arc_at_all)
{
  // No great circle through them is determined, so there is no arc to test and the answer has to be no.
  const std::vector<SphereCircle> none;
  EXPECT_FALSE(exposedGreatCircleArc(none, double3(1.0, 0.0, 0.0), double3(-1.0, 0.0, 0.0)));
}

TEST(boundary_components, an_arc_running_into_a_cap_is_not_exposed)
{
  // A cap over the north pole, and two points of the equator a third of the way round from one another. The
  // arc between them stays on the equator, which the cap does not reach.
  const std::vector<SphereCircle> circles = {capOf(double3(0.0, 0.0, 1.0), 0.5)};
  EXPECT_TRUE(exposedGreatCircleArc(circles, fromPolar(0.5 * std::numbers::pi, 0.0),
                                    fromPolar(0.5 * std::numbers::pi, 2.0)));

  // The same two points, but with the cap grown until it swallows the equator: now nothing of the arc is
  // exposed, its own ends least of all.
  const std::vector<SphereCircle> swallowed = {capOf(double3(0.0, 0.0, 1.0), 2.0)};
  EXPECT_FALSE(exposedGreatCircleArc(swallowed, fromPolar(0.5 * std::numbers::pi, 0.0),
                                     fromPolar(0.5 * std::numbers::pi, 2.0)));
}

TEST(boundary_components, an_arc_over_a_cap_that_neither_end_touches_is_caught)
{
  // The case the interior turning point is there for. Both ends stand clear of a cap on the equator, and the
  // arc between them passes straight over it: testing the ends alone would call this exposed.
  const std::vector<SphereCircle> circles = {capOf(double3(1.0, 0.0, 0.0), 0.4)};
  const double3 from = fromPolar(0.5 * std::numbers::pi, -1.0);
  const double3 to = fromPolar(0.5 * std::numbers::pi, 1.0);

  EXPECT_GT(std::acos(double3::dot(from, double3(1.0, 0.0, 0.0))), 0.4);
  EXPECT_GT(std::acos(double3::dot(to, double3(1.0, 0.0, 0.0))), 0.4);
  EXPECT_FALSE(exposedGreatCircleArc(circles, from, to));

  // Turned the long way round instead, the same two points are joined: the arc is the short one, so the pair
  // has to be taken past the cap rather than round it.
  const double3 far = fromPolar(0.5 * std::numbers::pi, std::numbers::pi);
  EXPECT_TRUE(exposedGreatCircleArc(circles, from, far));
  EXPECT_TRUE(exposedGreatCircleArc(circles, far, to));
}

TEST(boundary_components, the_predicate_does_not_depend_on_which_end_the_arc_is_taken_from)
{
  const std::vector<SphereCircle> circles = {capOf(double3(1.0, 0.0, 0.0), 0.4), capOf(double3(0.0, 1.0, 1.0), 0.6),
                                             capOf(double3(-1.0, 0.3, -0.2), 0.9)};

  // Over a grid of pairs, in both directions. The height along the arc is written from one end and swept
  // towards the other, so the two are different arithmetic on the same geometry.
  std::size_t exposed = 0;
  for (std::size_t i = 0; i < 12; ++i)
  {
    for (std::size_t j = 0; j < 12; ++j)
    {
      const double3 from = fromPolar(0.2 + 0.25 * static_cast<double>(i), 0.5 * static_cast<double>(j));
      const double3 to = fromPolar(0.3 + 0.2 * static_cast<double>(j), 0.4 * static_cast<double>(i) + 1.1);
      const bool forwards = exposedGreatCircleArc(circles, from, to);
      EXPECT_EQ(forwards, exposedGreatCircleArc(circles, to, from));
      if (forwards) ++exposed;
    }
  }

  // And the configuration is one where both answers occur, or the agreement above would say nothing.
  EXPECT_GT(exposed, 0uz);
  EXPECT_LT(exposed, 144uz);
}

TEST(boundary_components, an_arc_called_exposed_is_clear_of_every_cap_along_its_whole_length)
{
  // The claim the merging rests on, checked by walking the arc. Where the predicate says yes, no point of the
  // arc may lie inside a cap; where it says no, some point has to.
  const std::vector<SphereCircle> circles = {capOf(double3(1.0, 0.2, 0.0), 0.5), capOf(double3(0.1, 1.0, 0.4), 0.7),
                                             capOf(double3(-0.4, 0.2, 1.0), 0.45), capOf(double3(0.0, -1.0, 0.3), 0.8)};

  auto coveredSomewhere = [&](const double3& from, const double3& to)
  {
    constexpr std::size_t steps = 4000;
    for (std::size_t s = 0; s <= steps; ++s)
    {
      const double t = static_cast<double>(s) / static_cast<double>(steps);
      const double3 along = double3::normalize(from * (1.0 - t) + to * t);
      for (const SphereCircle& circle : circles)
      {
        if (double3::dot(along, circle.axis) > circle.cosineHalfAngle) return true;
      }
    }
    return false;
  };

  std::size_t agreed = 0;
  std::size_t checked = 0;
  for (std::size_t i = 0; i < 10; ++i)
  {
    for (std::size_t j = 0; j < 10; ++j)
    {
      const double3 from = fromPolar(0.15 + 0.3 * static_cast<double>(i), 0.6 * static_cast<double>(j));
      const double3 to = fromPolar(0.25 + 0.28 * static_cast<double>(j), 0.55 * static_cast<double>(i) + 2.0);

      // The endpoints have to be exposed for the question to be the one the merging asks.
      bool endsExposed = true;
      for (const SphereCircle& circle : circles)
      {
        if (double3::dot(from, circle.axis) > circle.cosineHalfAngle) endsExposed = false;
        if (double3::dot(to, circle.axis) > circle.cosineHalfAngle) endsExposed = false;
      }
      if (!endsExposed) continue;

      ++checked;
      if (exposedGreatCircleArc(circles, from, to) == !coveredSomewhere(from, to)) ++agreed;
    }
  }

  EXPECT_GT(checked, 20uz);
  EXPECT_EQ(agreed, checked);
}

TEST(boundary_components, a_path_through_a_way_point_joins_what_one_arc_cannot)
{
  // Two points either side of a cap, so that no arc between them is exposed, and a way point standing clear
  // of the cap on the other side. The path through it is two legs, and both of them are exposed.
  SphereBoundary boundary;
  boundary.circles = {capOf(double3(1.0, 0.0, 0.0), 0.4)};

  const double3 from = fromPolar(0.5 * std::numbers::pi, -1.0);
  const double3 to = fromPolar(0.5 * std::numbers::pi, 1.0);
  ASSERT_FALSE(exposedGreatCircleArc(boundary.circles, from, to));

  // With no way points there is nothing to go through, and the pair stays unjoined. That is the fallback the
  // nesting rule leaves behind, and it is allowed to prove less.
  EXPECT_FALSE(connectedOnSphere(boundary, from, to));

  boundary.wayPoints = {fromPolar(0.5 * std::numbers::pi, std::numbers::pi)};
  boundary.wayPointGroup = {0};
  EXPECT_TRUE(connectedOnSphere(boundary, from, to));
}

TEST(boundary_components, two_ends_reaching_different_groups_of_way_points_are_not_joined)
{
  // Caps over both poles leave a band round the equator, and two more cut that band in two at azimuth zero
  // and at azimuth pi. So there are two exposed regions with nothing joining them, one either side.
  const double equator = 0.5 * std::numbers::pi;
  SphereBoundary boundary;
  boundary.circles = {capOf(double3(0.0, 0.0, 1.0), 1.4), capOf(double3(0.0, 0.0, -1.0), 1.4),
                      capOf(double3(1.0, 0.0, 0.0), 0.3), capOf(double3(-1.0, 0.0, 0.0), 0.3)};

  const double3 from = fromPolar(equator, equator);
  const double3 to = fromPolar(equator, 3.0 * equator + 0.05);
  boundary.wayPoints = {fromPolar(equator, equator + 0.15), fromPolar(equator, 3.0 * equator - 0.15)};
  boundary.wayPointGroup = {0, 1};

  // Each end reaches the way point of its own region and neither reaches the other's, so the two groups are
  // reached one apiece and there is no path through them.
  ASSERT_FALSE(exposedGreatCircleArc(boundary.circles, from, to));
  ASSERT_TRUE(exposedGreatCircleArc(boundary.circles, from, boundary.wayPoints[0]));
  ASSERT_TRUE(exposedGreatCircleArc(boundary.circles, to, boundary.wayPoints[1]));
  ASSERT_FALSE(exposedGreatCircleArc(boundary.circles, from, boundary.wayPoints[1]));
  ASSERT_FALSE(exposedGreatCircleArc(boundary.circles, to, boundary.wayPoints[0]));

  EXPECT_FALSE(connectedOnSphere(boundary, from, to));
}

TEST(boundary_components, the_single_arc_is_tried_before_the_way_points)
{
  // Way points that are of no use, on a pair that one arc already joins. The answer is yes either way; what
  // this pins down is that a boundary carrying no way points is not thereby weaker on the easy case.
  SphereBoundary boundary;
  boundary.circles = {capOf(double3(0.0, 0.0, 1.0), 0.3)};

  const double3 from = fromPolar(0.5 * std::numbers::pi, 0.0);
  const double3 to = fromPolar(0.5 * std::numbers::pi, 0.5);
  EXPECT_TRUE(connectedOnSphere(boundary, from, to));

  boundary.wayPoints = {double3(0.0, 0.0, 1.0)};
  boundary.wayPointGroup = {0};
  EXPECT_TRUE(connectedOnSphere(boundary, from, to));
}
