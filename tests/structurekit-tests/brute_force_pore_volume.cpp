#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;
import brute_force_pore_volume;

// The volume by points thrown at the cell, and the list of lattice images the checks all walk.

// The points are dealt out in equal lanes, so asking for fewer points than there are lanes still has to throw a point
// per lane. Throwing none has nothing to divide by.
TEST(brute_force_pore_volume, fewer_points_than_lanes_still_throws_a_point_per_lane)
{
  BruteForceStructure structure =
      BruteForceStructure::create("one", UnitCell(20.0, 20.0, 20.0), {double3(10.0, 10.0, 10.0)}, {2.0});
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 1.0);

  for (std::size_t asked : {0uz, 1uz, 63uz})
  {
    BruteForcePoreVolume volume = BruteForcePoreVolume::compute(structure, voxels, asked);
    ASSERT_EQ(volume.numberOfPoints, 64uz) << "asked for " << asked;
    EXPECT_GE(volume.voidFraction, 0.0);
    EXPECT_LE(volume.voidFraction, 1.0);
  }
}

// Every lattice translation that carries an atom into contact with another has to be on the list: from the nearest
// image of the other atom, on to anywhere within the two radii. The pairs of a solvent-excluded surface are found by
// walking that list, and a cell thin along one axis has a good many such images along it, with plenty of void left
// around them all.
TEST(brute_force_pore_volume, the_image_list_reaches_every_image_in_contact)
{
  struct Case
  {
    const char *name;
    UnitCell cell;
    double radius;
  };
  const std::vector<Case> cases = {
      {"radius past the spread", UnitCell(2.0, 2.0, 2.0), 4.0},
      {"far past the spread", UnitCell(2.0, 2.0, 2.0), 12.0},
      {"thin cell, wide atom", UnitCell(2.0, 40.0, 40.0), 20.0},
      {"ordinary", UnitCell(10.0, 11.0, 12.0), 2.0}};

  RandomNumber random{std::optional<std::size_t>(31)};
  for (const Case &entry : cases)
  {
    BruteForceStructure structure =
        BruteForceStructure::create("cell", entry.cell, {0.5 * (entry.cell.cell * double3(1.0, 1.0, 1.0))}, {entry.radius});
    const double3 a = entry.cell.cell * double3(1.0, 0.0, 0.0);
    const double3 b = entry.cell.cell * double3(0.0, 1.0, 0.0);
    const double3 c = entry.cell.cell * double3(0.0, 0.0, 1.0);
    const double contact = 2.0 * entry.radius;

    std::set<std::array<long long, 3>> listed;
    for (const BruteForceStructure::Image &image : structure.images())
    {
      listed.insert({std::llround(1.0e6 * image.translation.x), std::llround(1.0e6 * image.translation.y),
                     std::llround(1.0e6 * image.translation.z)});
    }

    std::size_t wanted = 0;
    std::size_t present = 0;
    for (std::size_t pair = 0; pair < 8; ++pair)
    {
      const double3 elsewhere = entry.cell.cell * double3(random.uniform(), random.uniform(), random.uniform());
      const double3 base = structure.nearestImage(double3(0.0, 0.0, 0.0), elsewhere);

      const std::int32_t far = 40;
      for (std::int32_t k = -far; k <= far; ++k)
      {
        for (std::int32_t j = -far; j <= far; ++j)
        {
          for (std::int32_t i = -far; i <= far; ++i)
          {
            const double3 translation = static_cast<double>(i) * a + static_cast<double>(j) * b + static_cast<double>(k) * c;
            if ((base + translation).length() >= contact) continue;
            ++wanted;
            if (listed.contains({std::llround(1.0e6 * translation.x), std::llround(1.0e6 * translation.y),
                                 std::llround(1.0e6 * translation.z)}))
              ++present;
          }
        }
      }
    }
    EXPECT_GT(wanted, 0uz) << entry.name;
    EXPECT_EQ(present, wanted) << entry.name;
  }
}
