#include <gtest/gtest.h>

import std;

import int3;
import double3;
import simulationbox;
import voronoi_accessibility;
import exact_boundary_components;

// The boundary of the union of the balls, cut into its connected surfaces.
//
// Every case here has an answer known before it is run, from the arrangement rather than from the code: a
// sealed cage has an inside and an outside and nothing else, so its boundary is two surfaces and every one
// of its atoms carries one patch of each; a chain of balls round the periodic boundary is one surface that
// closes on a translate of itself; and moving any of it about the cell changes none of that.

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


VoronoiAccessibility geometryOf(const SimulationBox& box, const std::vector<double3>& positions, double ballRadius)
{
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const double3& position : positions)
  {
    fractionalPositions.push_back(double3::fract(box.inverseCell * position));
    radii.push_back(ballRadius);
  }
  return VoronoiAccessibility::create(box, fractionalPositions, radii, 0.0);
}


// The patches of one sphere counted without any of the machinery under test: a fine grid of directions,
// those outside every other ball kept, and the kept ones joined when they are neighbours on the grid. It is
// a raster, so it can miss a connection thinner than its spacing, which is why it is only ever asked about
// arrangements whose patches are broad.
std::vector<std::int32_t> lastFloodFill;
std::size_t lastFloodFillLatitudes = 0;

std::size_t patchesByFloodFill(const VoronoiAccessibility& accessibility, std::size_t atomIndex,
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

}  // namespace


// Two balls too far apart to touch: two surfaces, each of one patch, neither closing on a translate of
// itself. The simplest thing the decomposition can be asked, and the one that fails first if the joining
// across circles ever invents an edge.
TEST(boundary_components, separate_balls_are_separate_surfaces)
{
  SimulationBox box(30.0, 30.0, 30.0);
  VoronoiAccessibility geometry = geometryOf(box, {double3(8.0, 8.0, 8.0), double3(20.0, 20.0, 20.0)}, 2.0);

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
  VoronoiAccessibility geometry = geometryOf(box, {double3(14.0, 15.0, 15.0), double3(16.0, 15.0, 15.0)}, 2.0);

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
  VoronoiAccessibility geometry =
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
  VoronoiAccessibility geometry = geometryOf(box, positions, radius);

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
  VoronoiAccessibility geometry = geometryOf(box, positions, 1.8);

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
  VoronoiAccessibility geometry = geometryOf(box, positions, 1.8);

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
  VoronoiAccessibility geometry = geometryOf(box, cageCentres(double3(0.0, 0.0, 0.0), 3.0, count), 1.8);

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
  VoronoiAccessibility geometry = geometryOf(box, positions, 1.8);

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
  VoronoiAccessibility geometry = geometryOf(box, positions, 1.5);  // spacing 2.0, so consecutive balls overlap

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
  VoronoiAccessibility geometry =
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
  VoronoiAccessibility geometry = geometryOf(box, positions, radius);

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

  std::vector<std::pair<std::string, VoronoiAccessibility>> cases;
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
  VoronoiAccessibility geometry =
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
