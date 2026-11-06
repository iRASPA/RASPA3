#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;
import brute_force_pore_spikes;

// The pore-size spikes read off the bare void: the maximal spheres the walk uphill finds, and which families
// of them are reported.

namespace
{

// A tube of rings of six atoms about an axis through the middle of the cell. `radiusAt` gives the ring's
// radius a distance `t` along the axis, for t from `from` to `to` in steps of `spacing`; the rings' atoms all
// sit at the same six azimuths, so the axis is the middle of the tube everywhere along it.
struct Tube
{
  UnitCell box;
  double3 middle;
  double3 axis;
  std::vector<double3> positions;
  std::vector<double> radii;

  Tube(double side, double3 along) : box(side, side, side), middle(0.5 * side, 0.5 * side, 0.5 * side)
  {
    this->axis = (1.0 / along.length()) * along;
  }

  void rings(double from, double to, double spacing, double atomRadius, auto &&radiusAt)
  {
    double3 u = double3::cross(this->axis, double3(0.0, 0.0, 1.0));
    if (u.length() < 0.5) u = double3::cross(this->axis, double3(1.0, 0.0, 0.0));
    u = (1.0 / u.length()) * u;
    double3 v = double3::cross(this->axis, u);
    for (double t = from; t <= to + 1.0e-9; t += spacing)
    {
      double radius = radiusAt(t);
      for (int k = 0; k < 6; ++k)
      {
        double azimuth = k * std::numbers::pi / 3.0;
        this->positions.push_back(this->middle + t * this->axis +
                                  (radius * std::cos(azimuth)) * u + (radius * std::sin(azimuth)) * v);
        this->radii.push_back(atomRadius);
      }
    }
  }

  BruteForceStructure structure() const
  {
    return BruteForceStructure::create("tube", this->box, this->positions, this->radii);
  }
};

// The most room anywhere on the axis between `from` and `to`, by a scan refined three times over.
double mostRoomOnTheAxis(const BruteForceStructure &structure, const Tube &tube, double from, double to)
{
  double bestT = from;
  double best = -1.0;
  for (double step : {1.0e-3, 1.0e-6, 1.0e-9})
  {
    double low = (best < 0.0) ? from : bestT - 1000.0 * step;
    double high = (best < 0.0) ? to : bestT + 1000.0 * step;
    for (double t = low; t <= high; t += step)
    {
      double room = structure.clearance(tube.middle + t * tube.axis);
      if (room > best)
      {
        best = room;
        bestT = t;
      }
    }
  }
  return best;
}

// A sealed cage of `count` spheres on a shell about `centre`, placed on the Fibonacci spiral.
std::vector<double3> fibonacciShell(const double3 &centre, double shellRadius, std::size_t count)
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

}  // namespace


TEST(brute_force_pore_spikes, the_walk_climbs_a_crease_that_runs_between_the_lattice_directions)
{
  // A tube that widens slowly along an axis pointing nowhere in particular, closed off by a narrower ring at
  // its wide end. The rings are close enough that the room on the axis rises all the way to the closing ring
  // and nowhere in between, so the only maximal sphere in the tube sits where the widening meets the closing.
  // Along the axis the room rises by five parts in a hundred, and any step off the axis costs its full length
  // in room, so the only way up is along the axis to within a few degrees. No lattice direction lies that
  // close to it, and a walk confined to the lattice directions stops on the axis wherever it first reaches it.
  Tube tube(20.0, double3(1.0, 2.0, 3.0));
  tube.rings(-6.0, 6.0, 0.5, 1.5, [](double t) { return 2.5 + 0.05 * (t + 6.0); });
  tube.rings(6.5, 7.0, 0.5, 1.5, [](double) { return 2.0; });
  BruteForceStructure structure = tube.structure();
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.5);

  const double room = mostRoomOnTheAxis(structure, tube, -6.0, 7.0);
  ASSERT_GT(room, 1.0);

  BruteForcePoreSpikes spikes = BruteForcePoreSpikes::compute(structure, voxels, 100000, 8);

  bool found = false;
  for (const BruteForceSpikeFamily &family : spikes.families)
  {
    if (std::abs(family.diameter - 2.0 * room) <= 1.0e-6) found = true;
  }
  EXPECT_TRUE(found) << "no family has the diameter " << 2.0 * room << " of the sphere that fits at the top of the tube";
}


TEST(brute_force_pore_spikes, families_beyond_the_dozen_widest_are_counted_and_the_widest_above_the_floor_reported)
{
  // A tube of rings two steps apart, each a little wider than the last: the room on the axis bulges between
  // every pair of rings, and each bulge is a maximal sphere of its own, so the tube alone gives a dozen
  // families of distinct diameters. Beside it, eight sealed cages, all alike and narrower than any bulge, which
  // are one family. There are more families than are weighed in full without a first look, so the first look
  // is taken; every family is still counted, and the ones reported are the widest with any void to their name.
  Tube tube(30.0, double3(1.0, 0.0, 0.0));
  tube.rings(-11.0, 11.0, 2.0, 1.5, [](double t) { return 3.1 + 0.05 * (t + 11.0); });
  for (double3 corner : {double3(4.0, 4.0, 4.0), double3(26.0, 4.0, 4.0), double3(4.0, 26.0, 4.0),
                         double3(26.0, 26.0, 4.0), double3(4.0, 4.0, 26.0), double3(26.0, 4.0, 26.0),
                         double3(4.0, 26.0, 26.0), double3(26.0, 26.0, 26.0)})
  {
    for (const double3 &position : fibonacciShell(corner, 3.0, 20))
    {
      tube.positions.push_back(position);
      tube.radii.push_back(1.5);
    }
  }
  BruteForceStructure structure = tube.structure();
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.75);

  BruteForcePoreSpikes spikes = BruteForcePoreSpikes::compute(structure, voxels, 100000, 3);
  ASSERT_GT(spikes.numberOfFamilies, 12uz) << "the bulges alone were to give a dozen families";
  ASSERT_EQ(spikes.families.size(), 3uz);
  for (std::size_t f = 0; f + 1 < spikes.families.size(); ++f)
  {
    EXPECT_GT(spikes.families[f].diameter, spikes.families[f + 1].diameter);
  }
  for (const BruteForceSpikeFamily &family : spikes.families) EXPECT_GT(family.weight, 1.0e-5);

  // Asked for every family, the cages come too: eight balls of radius one and a half, one family of diameter
  // three, the narrowest of all. Every family with void to its name is reported, which need not be every place
  // a walk ended: one that ended short of a maximum names a diameter no void has for its pore size.
  BruteForcePoreSpikes all = BruteForcePoreSpikes::compute(structure, voxels, 100000, 32);
  ASSERT_GT(all.families.size(), 12uz);
  ASSERT_LE(all.families.size(), all.numberOfFamilies);
  const BruteForceSpikeFamily &cages = all.families.back();
  EXPECT_NEAR(cages.diameter, 3.0, 1.0e-6);
  EXPECT_EQ(cages.centres.size(), 8uz);
  EXPECT_GT(cages.weight, 1.0e-3);
}


namespace
{

// A solid block of atoms on a cubic lattice, overlapping enough to leave no room anywhere, with single atoms taken out.
// What is left where one was is a hole whose size is set by the radii of the six atoms around it, and the holes are
// nothing else's business, so the void is exactly the holes and which family holds what is known beforehand.
struct Block
{
  static constexpr std::int32_t sites = 12;
  static constexpr double spacing = 3.0;
  static constexpr double solidRadius = 2.7;

  std::map<std::array<std::int32_t, 3>, double> radiusAt;
  std::set<std::array<std::int32_t, 3>> holes;

  static std::int32_t wrap(std::int32_t i) { return (i + sites) % sites; }

  void hole(std::array<std::int32_t, 3> site, double room)
  {
    holes.insert(site);
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      for (std::int32_t step : {-1, 1})
      {
        std::array<std::int32_t, 3> neighbour = site;
        neighbour[axis] = wrap(neighbour[axis] + step);
        radiusAt[neighbour] = spacing - room;
      }
    }
  }

  BruteForceStructure structure() const
  {
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::int32_t k = 0; k < sites; ++k)
    {
      for (std::int32_t j = 0; j < sites; ++j)
      {
        for (std::int32_t i = 0; i < sites; ++i)
        {
          const std::array<std::int32_t, 3> site{i, j, k};
          if (holes.contains(site)) continue;
          // A quarter of a voxel off the lattice, which puts a voxel centre on every site.
          positions.push_back(double3(spacing * i + 0.25, spacing * j + 0.25, spacing * k + 0.25));
          const auto found = radiusAt.find(site);
          radii.push_back(found == radiusAt.end() ? solidRadius : found->second);
        }
      }
    }
    const double side = spacing * static_cast<double>(sites);
    return BruteForceStructure::create("block", UnitCell(side, side, side), positions, radii);
  }
};

}  // namespace


// One wide hole, and thirteen families of narrower holes that hold more of the void each than the wide one does. The
// widest above the floor are what is reported, and how much void a family holds is no reason to pass over a wider
// one. Nor can its void go elsewhere when it is passed over: what each family holds is what it holds with all the
// wider ones in place, the same as when every family is asked for.
TEST(brute_force_pore_spikes, a_wide_family_that_holds_little_is_reported_ahead_of_narrower_ones_that_hold_more)
{
  constexpr double wideRoom = 1.3;
  Block block;
  block.hole({0, 0, 0}, wideRoom);

  constexpr std::size_t holesPerFamily = 4;
  std::int32_t slot = 1;
  for (std::size_t family = 0; family < 13; ++family)
  {
    const double room = 1.2 - 0.015 * static_cast<double>(family);
    for (std::size_t h = 0; h < holesPerFamily; ++h, ++slot)
    {
      block.hole({3 * (slot % 4), 3 * ((slot / 4) % 4), 3 * (slot / 16)}, room);
    }
  }

  BruteForceStructure structure = block.structure();
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.5);

  constexpr std::size_t points = 300000;
  BruteForcePoreSpikes all = BruteForcePoreSpikes::compute(structure, voxels, points, 64);
  ASSERT_GE(all.numberOfFamilies, 14uz);
  ASSERT_GE(all.families.size(), 14uz);
  EXPECT_NEAR(all.families[0].diameter, 2.0 * wideRoom, 1.0e-6);
  EXPECT_EQ(all.families[0].centres.size(), 1uz);

  // Enough of the points land in the void for the weights to mean something, and the wide family really is lighter
  // than the dozen it is wider than, which is the whole of the case.
  ASSERT_GT(all.voidFraction * static_cast<double>(points), 1000.0);
  std::size_t heavier = 0;
  double total = 0.0;
  for (const BruteForceSpikeFamily &family : all.families)
  {
    total += family.weight;
    if (family.weight > all.families[0].weight + 2.0 * family.weightError) ++heavier;
  }
  ASSERT_GE(heavier, 12uz);
  EXPECT_LE(total, 1.0 + 1.0e-12) << "a point of the void belongs to one family at most";

  BruteForcePoreSpikes widest = BruteForcePoreSpikes::compute(structure, voxels, points, 3);
  ASSERT_EQ(widest.families.size(), 3uz);
  EXPECT_NEAR(widest.voidVolume, all.voidVolume, 1.0e-9 * all.voidVolume);
  for (std::size_t f = 0; f < 3; ++f)
  {
    EXPECT_NEAR(widest.families[f].diameter, all.families[f].diameter, 1.0e-6) << "family " << f;

    // The same stream of points and the same families down to this one, so the same count to the point.
    EXPECT_DOUBLE_EQ(widest.families[f].weight, all.families[f].weight) << "family " << f;
  }
}


namespace
{

// Whether a sphere wider than `radius` that fits in the void holds `point`, by a search that owes nothing to the
// grid or to the way the spikes are weighed: centres are thrown into the ball of that radius about the point, and
// the roomiest few are climbed from without leaving the ball.
bool heldByAWiderSphere(const BruteForceStructure &structure, const double3 &point, double radius,
                        RandomNumber &random)
{
  struct Centre
  {
    double3 offset;
    double room{-1.0e300};
  };
  std::array<Centre, 4> roomiest{};

  for (std::size_t n = 0; n < 300; ++n)
  {
    double3 offset;
    do
    {
      offset = double3(2.0 * random.uniform() - 1.0, 2.0 * random.uniform() - 1.0, 2.0 * random.uniform() - 1.0);
    } while (double3::dot(offset, offset) > 1.0);
    offset = radius * offset;

    const double room = structure.clearance(point + offset);
    if (room > radius) return true;
    auto worst = std::ranges::min_element(roomiest, {}, &Centre::room);
    if (room > worst->room) *worst = Centre{offset, room};
  }

  auto intoTheBall = [&](const double3 &offset)
  {
    const double length = offset.length();
    return length > radius ? (radius / length) * offset : offset;
  };
  for (Centre centre : roomiest)
  {
    for (double step = 0.2; step > 1.0e-4;)
    {
      bool moved = false;
      for (const double3 &direction : {double3(1.0, 0.0, 0.0), double3(-1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0),
                                       double3(0.0, -1.0, 0.0), double3(0.0, 0.0, 1.0), double3(0.0, 0.0, -1.0)})
      {
        const double3 offset = intoTheBall(centre.offset + step * direction);
        const double room = structure.clearance(point + offset);
        if (room > radius) return true;
        if (room > centre.room)
        {
          centre = Centre{offset, room};
          moved = true;
        }
      }
      if (!moved) step *= 0.5;
    }
  }
  return false;
}

}  // namespace


// Two sealed cavities that open into one another through a neck narrower than either: the larger three and a half
// Ångström in radius, the smaller two and a half, the neck a little under that. The largest sphere of each is a
// family. A sphere the size of the smaller cavity's also fits on the larger one's side of the neck, where it
// reaches through and holds part of the smaller cavity, and so does every sphere between that size and the larger
// cavity's own. The void those hold has a pore size above the smaller family's, and is not that family's, though
// the only maximal sphere that holds it is. Charged to it, the smaller family is a quarter too heavy.
TEST(brute_force_pore_spikes, void_a_wider_sphere_reaches_through_a_neck_is_not_the_narrower_cavitys)
{
  const double3 large(5.2, 5.5, 5.5);
  const double3 small(9.0, 5.5, 5.5);
  const double shellOfLarge = 4.1;
  const double shellOfSmall = 3.1;
  const double atomRadius = 0.6;

  // Each cavity is walled by small spheres on a shell about its middle, left out where they would stand inside
  // the other cavity's shell, which is what opens the neck.
  std::vector<double3> positions;
  for (const double3 &position : fibonacciShell(large, shellOfLarge, 340))
  {
    if ((position - small).length() > shellOfSmall) positions.push_back(position);
  }
  for (const double3 &position : fibonacciShell(small, shellOfSmall, 200))
  {
    if ((position - large).length() > shellOfLarge) positions.push_back(position);
  }
  std::vector<double> radii(positions.size(), atomRadius);

  BruteForceStructure structure =
      BruteForceStructure::create("two cavities", UnitCell(13.4, 11.0, 11.0), positions, radii);
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.25);

  BruteForcePoreSpikes spikes = BruteForcePoreSpikes::compute(structure, voxels, 120000, 32);

  // The two families, by where their one sphere stands.
  const BruteForceSpikeFamily *ofLarge = nullptr;
  const BruteForceSpikeFamily *ofSmall = nullptr;
  for (const BruteForceSpikeFamily &family : spikes.families)
  {
    for (const double3 &centre : family.centres)
    {
      if ((centre - large).length() < 0.2) ofLarge = &family;
      if ((centre - small).length() < 0.2) ofSmall = &family;
    }
  }
  ASSERT_NE(ofLarge, nullptr) << "no family has its sphere in the middle of the larger cavity";
  ASSERT_NE(ofSmall, nullptr) << "no family has its sphere in the middle of the smaller cavity";
  ASSERT_EQ(ofLarge->centres.size(), 1uz);
  ASSERT_EQ(ofSmall->centres.size(), 1uz);
  const double radiusOfLarge = 0.5 * ofLarge->diameter;
  const double radiusOfSmall = 0.5 * ofSmall->diameter;
  EXPECT_NEAR(radiusOfLarge, shellOfLarge - atomRadius, 0.05);
  EXPECT_NEAR(radiusOfSmall, shellOfSmall - atomRadius, 0.05);

  // The larger family holds its whole sphere: nothing wider fits in there.
  const double ball = 4.0 * std::numbers::pi / 3.0;
  const double heldByLarge = ofLarge->weight * spikes.voidVolume;
  EXPECT_NEAR(heldByLarge, ball * std::pow(radiusOfLarge, 3), 0.03 * ball * std::pow(radiusOfLarge, 3));

  // The smaller family holds what is left of its sphere once the points a wider sphere holds are taken out,
  // which is asked of each point of the sphere in turn.
  RandomNumber random{7};
  constexpr std::size_t thrown = 4000;
  std::size_t leftToTheSmaller = 0;
  const double3 middle = ofSmall->centres.front();
  for (std::size_t n = 0; n < thrown; ++n)
  {
    double3 offset;
    do
    {
      offset = double3(2.0 * random.uniform() - 1.0, 2.0 * random.uniform() - 1.0, 2.0 * random.uniform() - 1.0);
    } while (double3::dot(offset, offset) > 1.0);
    if (!heldByAWiderSphere(structure, middle + radiusOfSmall * offset, radiusOfSmall + 1.0e-4, random))
    {
      ++leftToTheSmaller;
    }
  }
  const double sphereOfSmall = ball * std::pow(radiusOfSmall, 3);
  const double expected = sphereOfSmall * static_cast<double>(leftToTheSmaller) / static_cast<double>(thrown);
  const double heldBySmall = ofSmall->weight * spikes.voidVolume;
  EXPECT_NEAR(heldBySmall, expected, 0.08 * expected)
      << "the smaller family holds " << heldBySmall << " cubic Angstrom and a search point by point leaves it "
      << expected;

  // And that is well short of what the two maximal spheres alone would say: the smaller sphere less the lens it
  // shares with the larger.
  const double apart = (ofLarge->centres.front() - middle).length();
  const double lens = std::numbers::pi * std::pow(radiusOfLarge + radiusOfSmall - apart, 2) *
                      (apart * apart + 2.0 * apart * (radiusOfLarge + radiusOfSmall) -
                       3.0 * std::pow(radiusOfLarge - radiusOfSmall, 2)) /
                      (12.0 * apart);
  EXPECT_LT(heldBySmall, 0.9 * (sphereOfSmall - lens))
      << "the smaller family holds " << heldBySmall << " cubic Angstrom, and its sphere less the lens is "
      << sphereOfSmall - lens;
}


// The atoms near a point, which is what the search for a wider sphere asks the room of in place of the whole
// structure: every atom of every image whose surface is within reach, and none other, in a cell far from square.
TEST(brute_force_structure, the_atoms_within_reach_are_those_a_search_over_every_image_finds)
{
  const double degrees = std::numbers::pi / 180.0;
  UnitCell cell(9.0, 10.0, 11.0, 70.0 * degrees, 80.0 * degrees, 65.0 * degrees);

  RandomNumber random{3};
  std::vector<double3> positions;
  std::vector<double> radii;
  for (std::size_t n = 0; n < 12; ++n)
  {
    positions.push_back(cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(0.8 + 0.6 * random.uniform());
  }
  BruteForceStructure structure = BruteForceStructure::create("skewed", cell, positions, radii);

  std::vector<BruteForceStructure::Neighbour> found;
  for (double reach : {0.5, 3.0, 7.5, 12.0})
  {
    for (std::size_t n = 0; n < 5; ++n)
    {
      const double3 point = cell.cell * double3(random.uniform(), random.uniform(), random.uniform());
      ASSERT_TRUE(structure.atomsWithin(point, reach, found));

      std::vector<double> surfaces;
      for (const BruteForceStructure::Neighbour &atom : found)
      {
        surfaces.push_back(atom.fromCentre.length() - atom.radius);
      }
      std::ranges::sort(surfaces);

      std::vector<double> expected;
      for (std::size_t j = 0; j < positions.size(); ++j)
      {
        for (std::int32_t a = -8; a <= 8; ++a)
        {
          for (std::int32_t b = -8; b <= 8; ++b)
          {
            for (std::int32_t c = -8; c <= 8; ++c)
            {
              const double3 image = positions[j] + cell.cell * double3(static_cast<double>(a), static_cast<double>(b),
                                                                       static_cast<double>(c));
              const double surface = (point - image).length() - radii[j];
              if (surface <= reach) expected.push_back(surface);
            }
          }
        }
      }
      std::ranges::sort(expected);

      ASSERT_EQ(surfaces.size(), expected.size()) << "reach " << reach;
      for (std::size_t i = 0; i < surfaces.size(); ++i) EXPECT_NEAR(surfaces[i], expected[i], 1.0e-9);

      // and the room among them is the room, wherever it is no more than the reach
      if (!surfaces.empty()) EXPECT_NEAR(surfaces.front(), structure.clearance(point), 1.0e-9);
    }
  }

  // Past what the images held are known to cover, the list would be short, and none is given.
  EXPECT_FALSE(structure.atomsWithin(positions.front(), 1.0e3, found));
  EXPECT_TRUE(found.empty());
}


// A ray leaves the atoms it starts inside in the direction it was sent. The start is put off the centre of an
// atom, so that the way out is shorter one way than the other and the two cannot be taken for each other: sent
// along +x from a point one Ångström that side of the centre of a sphere of radius two, the ray is out after
// one, and sent along -x after three. It used to be the other way round. The displacement worked out for each
// atom runs from the atom to the ray's start, and it was turned round before the ray was laid along it, which
// is the ray laid the other way.
TEST(brute_force_structure, a_ray_leaves_the_atoms_the_way_it_was_sent)
{
  const UnitCell cell(20.0, 20.0, 20.0);
  const double3 centre(5.0, 5.0, 5.0);

  // One sphere, the start a step off its centre along each axis in turn.
  {
    BruteForceStructure structure = BruteForceStructure::create("one sphere", cell, {centre}, {2.0});
    const std::array<double3, 3> axes{double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0)};
    for (const double3 &axis : axes)
    {
      const double3 start = centre + 1.0 * axis;
      const std::optional<double> forward = structure.firstExit(start, axis, 10.0);
      const std::optional<double> backward = structure.firstExit(start, -1.0 * axis, 10.0);
      ASSERT_TRUE(forward.has_value());
      ASSERT_TRUE(backward.has_value());
      EXPECT_NEAR(forward.value(), 1.0, 1.0e-12);
      EXPECT_NEAR(backward.value(), 3.0, 1.0e-12);
    }

    // In any direction, where it comes out is on the sphere and it was inside all the way there.
    RandomNumber random{11};
    for (std::size_t n = 0; n < 200; ++n)
    {
      const double3 start = centre + double3(0.9 * (random.uniform() - 0.5), 1.1 * (random.uniform() - 0.5),
                                             1.3 * (random.uniform() - 0.5));
      const double3 direction = random.randomVectorOnUnitSphere();
      const std::optional<double> out = structure.firstExit(start, direction, 10.0);
      ASSERT_TRUE(out.has_value());
      EXPECT_NEAR((start + out.value() * direction - centre).length(), 2.0, 1.0e-9) << "ray " << n;
      EXPECT_LT((start + 0.5 * out.value() * direction - centre).length(), 2.0) << "ray " << n;
    }

    // A start outside every atom is out already, whichever way it looks.
    EXPECT_EQ(structure.firstExit(centre + double3(3.0, 0.0, 0.0), double3(-1.0, 0.0, 0.0), 10.0).value(), 0.0);

    // And one that is not out within the length allowed is not out.
    EXPECT_FALSE(structure.firstExit(centre + double3(1.0, 0.0, 0.0), double3(-1.0, 0.0, 0.0), 2.5).has_value());
    EXPECT_TRUE(structure.firstExit(centre + double3(1.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), 2.5).has_value());
  }

  // Two spheres in a row that overlap, and a third further on that does not: sent towards them the ray is
  // inside the first two without a break and comes out between the second and the third; sent away it comes
  // out of the first alone.
  {
    BruteForceStructure structure = BruteForceStructure::create(
        "a row", cell, {centre, centre + double3(3.0, 0.0, 0.0), centre + double3(9.0, 0.0, 0.0)}, {2.0, 2.0, 2.0});
    const std::optional<double> along = structure.firstExit(centre, double3(1.0, 0.0, 0.0), 15.0);
    const std::optional<double> away = structure.firstExit(centre, double3(-1.0, 0.0, 0.0), 15.0);
    ASSERT_TRUE(along.has_value());
    ASSERT_TRUE(away.has_value());
    EXPECT_NEAR(along.value(), 5.0, 1.0e-12);
    EXPECT_NEAR(away.value(), 2.0, 1.0e-12);
  }

  // The copies of an atom count: a sphere across the face of the cell is met by a ray that leaves the cell.
  {
    BruteForceStructure structure =
        BruteForceStructure::create("across the face", cell, {double3(1.0, 5.0, 5.0), double3(18.0, 5.0, 5.0)},
                                    {2.0, 2.0});
    // From inside the first, towards -x: out of it at x = -1, which is inside the copy of the second at x = -2
    // reaching from -4 to 0, so the ray is not out until x = -4.
    const std::optional<double> out = structure.firstExit(double3(1.5, 5.0, 5.0), double3(-1.0, 0.0, 0.0), 15.0);
    ASSERT_TRUE(out.has_value());
    EXPECT_NEAR(out.value(), 5.5, 1.0e-12);
    // The other way there is nothing but the first sphere.
    const std::optional<double> other = structure.firstExit(double3(1.5, 5.0, 5.0), double3(1.0, 0.0, 0.0), 15.0);
    ASSERT_TRUE(other.has_value());
    EXPECT_NEAR(other.value(), 1.5, 1.0e-12);
  }
}

// The displacement an atom is reported with runs from the atom to the point, so the point less it is where
// the atom is.
TEST(brute_force_structure, an_atom_within_reach_is_reported_from_its_centre_to_the_point)
{
  const UnitCell cell(12.0, 13.0, 14.0);
  const std::vector<double3> positions{double3(1.0, 2.0, 3.0), double3(11.5, 12.0, 0.5)};
  BruteForceStructure structure = BruteForceStructure::create("two", cell, positions, {1.0, 1.0});

  const double3 point(0.5, 0.5, 13.5);
  std::vector<BruteForceStructure::Neighbour> found;
  ASSERT_TRUE(structure.atomsWithin(point, 3.0, found));
  ASSERT_FALSE(found.empty());
  for (const BruteForceStructure::Neighbour &atom : found)
  {
    // Some copy of one of the two atoms stands exactly there.
    const double3 where = point - atom.fromCentre;
    double nearest = std::numeric_limits<double>::max();
    for (const double3 &position : positions)
    {
      nearest = std::min(nearest, structure.nearestImage(position, where).length());
    }
    EXPECT_LT(nearest, 1.0e-9);
  }
}
