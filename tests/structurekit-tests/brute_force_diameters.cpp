#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;
import brute_force_diameters;

// The widest way through a crystal, and what is read off the piece of the void it belongs to.
//
// The crystal is a stack of layers of spheres on a square lattice, wide enough to touch along the sides of
// the squares. Each square of four is a window, the stack of windows over one another is a tube along z, and
// between the layers the tubes see each other only through narrower windows across. So at the width of its
// windows a tube is a channel that runs one way, and only at a lower width does the void run every way.

namespace
{

constexpr double sphereRadius = 4.5;
constexpr double capRadius = 2.2;

// `across` tubes along x and along y. The tubes listed in `capped` have a sphere in one of their windows,
// which seals them: a tube with a window shut is a pocket, from that window to its next image.
BruteForceStructure tubes(std::int32_t across, const std::vector<int3> &capped)
{
  double side = 8.0 * static_cast<double>(across);
  std::vector<double3> positions;
  std::vector<double> radii;
  for (std::int32_t a = 0; a < across; ++a)
  {
    for (std::int32_t b = 0; b < across; ++b)
    {
      for (double z : {3.0, 9.0})
      {
        positions.push_back(double3(4.0 + 8.0 * a, 4.0 + 8.0 * b, z));
        radii.push_back(sphereRadius);
      }
    }
  }
  for (const int3 &tube : capped)
  {
    positions.push_back(double3(8.0 * tube.x, 8.0 * tube.y, 3.0));
    radii.push_back(capRadius);
  }
  return BruteForceStructure::create("tubes", UnitCell(side, side, 12.0), positions, radii);
}

BruteForceDiameters diametersOf(const BruteForceStructure &structure)
{
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.5);
  return BruteForceDiameters::compute(structure, voxels);
}

// The room at the centre of a tube, between two layers: the centre of a cube of eight spheres.
const double tubeRoom = std::sqrt(41.0) - sphereRadius;
// A window as the grid sees it: the voxel centres sit a quarter spacing off the tube's axis.
const double gridWindow = std::sqrt(3.75 * 3.75 + 3.75 * 3.75 + 0.25 * 0.25) - sphereRadius;

}  // namespace


TEST(brute_force_diameters, a_channel_that_runs_one_way_at_its_bottleneck_is_reported_as_running_one_way)
{
  BruteForceDiameters diameters = diametersOf(tubes(2, {}));

  ASSERT_TRUE(diameters.percolates);
  EXPECT_EQ(diameters.dimensionality, 1uz);
  EXPECT_NEAR(diameters.freeSphereDiameter, 2.0 * gridWindow, 1.0e-5);
  EXPECT_NEAR(diameters.includedSphereDiameter, 2.0 * tubeRoom, 1.0e-6);
  EXPECT_NEAR(diameters.includedAlongFreePathDiameter, 2.0 * tubeRoom, 1.0e-6);
}


TEST(brute_force_diameters, a_pocket_as_wide_at_one_point_as_the_bottleneck_is_not_where_the_path_goes)
{
  // Fifteen tubes sealed, one open. Each sealed tube keeps a window exactly as wide as the open tube's, and
  // which of all those tied hops sorts last is nobody's choice; Dif is to come from the open tube whichever
  // it is. Which tube is the open one changes the order the hops are generated in, and with it the order
  // the ties come out in.
  for (std::int32_t open = 0; open < 4; ++open)
  {
    std::vector<int3> capped;
    for (std::int32_t a = 0; a < 4; ++a)
    {
      for (std::int32_t b = 0; b < 4; ++b)
      {
        if (!(a == open && b == open)) capped.push_back(int3(a, b, 0));
      }
    }

    BruteForceDiameters diameters = diametersOf(tubes(4, capped));

    ASSERT_TRUE(diameters.percolates) << "open tube " << open;
    EXPECT_EQ(diameters.dimensionality, 1uz) << "open tube " << open;
    EXPECT_NEAR(diameters.freeSphereDiameter, 2.0 * gridWindow, 1.0e-5) << "open tube " << open;
    EXPECT_NEAR(diameters.includedAlongFreePathDiameter, 2.0 * tubeRoom, 1.0e-6) << "open tube " << open;
  }
}


TEST(brute_force_diameters, the_guaranteed_width_allows_for_the_longest_hop_between_slanted_voxels)
{
  // One atom in a hexagonal cell, so the void runs away every way. The guaranteed width is the width less half
  // the longest hop, and the hops are the thirteen steps to neighbouring voxels, the longest of which runs
  // along a voxel's long diagonal: on a cell with an obtuse angle that is longer than the voxel's edges in
  // quadrature, and longer again than the widths across it.
  const double right = 0.5 * std::numbers::pi;
  UnitCell box(12.0, 12.0, 12.0, right, right, 2.0 * std::numbers::pi / 3.0);
  BruteForceStructure structure = BruteForceStructure::create("one atom", box, {double3(0.0, 0.0, 0.0)}, {2.0});
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.5);
  BruteForceDiameters diameters = BruteForceDiameters::compute(structure, voxels);

  ASSERT_TRUE(diameters.percolates);

  const std::array<int3, 13> steps{int3(1, 0, 0),  int3(0, 1, 0),  int3(0, 0, 1),   int3(1, 1, 0),   int3(1, -1, 0),
                                   int3(1, 0, 1),  int3(1, 0, -1), int3(0, 1, 1),   int3(0, 1, -1),  int3(1, 1, 1),
                                   int3(1, 1, -1), int3(1, -1, 1), int3(1, -1, -1)};
  double longestHop = 0.0;
  for (const int3 &step : steps)
  {
    double3 hop = box.cell * double3(static_cast<double>(step.x) / static_cast<double>(voxels.counts.x),
                                     static_cast<double>(step.y) / static_cast<double>(voxels.counts.y),
                                     static_cast<double>(step.z) / static_cast<double>(voxels.counts.z));
    longestHop = std::max(longestHop, hop.length());
  }
  EXPECT_NEAR(diameters.guaranteedFreeSphereDiameter, diameters.freeSphereDiameter - longestHop, 1.0e-9);
}


// Two tubes that never meet at the width of their windows, the same windows in both, and in one of them a small
// sphere at the middle of the room. Whichever tube that is, the other still has the whole room, and it is reachable
// at the bottleneck's width, so that is Dif. Which tube's loop the sort put first must not decide it.
TEST(brute_force_diameters, the_room_of_either_tube_tied_at_the_bottleneck_is_on_the_path)
{
  for (const int3 &smaller : {int3(0, 0, 0), int3(1, 1, 0)})
  {
    double side = 16.0;
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::int32_t a = 0; a < 2; ++a)
    {
      for (std::int32_t b = 0; b < 2; ++b)
      {
        for (double z : {3.0, 9.0})
        {
          positions.push_back(double3(4.0 + 8.0 * a, 4.0 + 8.0 * b, z));
          radii.push_back(sphereRadius);
        }
      }
    }
    // Only one tube of each pair round the lattice is left whole.
    for (const int3 &filled : {smaller, int3(smaller.x, 1 - smaller.y, 0), int3(1 - smaller.x, smaller.y, 0)})
    {
      // A tube has two rooms to a period, between the layers on either side.
      for (double z : {0.0, 6.0})
      {
        positions.push_back(double3(8.0 * filled.x, 8.0 * filled.y, z));
        radii.push_back(0.05);
      }
    }
    BruteForceStructure structure = BruteForceStructure::create("tubes", UnitCell(side, side, 12.0), positions, radii);
    BruteForceDiameters diameters = diametersOf(structure);

    ASSERT_TRUE(diameters.percolates);
    EXPECT_EQ(diameters.dimensionality, 1uz);
    EXPECT_NEAR(diameters.freeSphereDiameter, 2.0 * gridWindow, 1.0e-5);
    EXPECT_NEAR(diameters.includedAlongFreePathDiameter, 2.0 * tubeRoom, 1.0e-6) << "smaller from " << smaller.x;
  }
}

// One piece that runs one way and one that runs two, tied at the same width. The widest path can run two ways, and
// that is what is reported, wherever the sort puts them. Nor do two pieces add up: see the two straws below.
TEST(brute_force_diameters, the_dimensionality_at_a_tie_is_that_of_the_piece_that_runs_most_ways)
{
  BruteForceStructure structure = tubes(2, {});
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.5);
  const float width = 1.0f;
  auto carve = [&](auto &&isOpen)
  {
    std::size_t open = 0;
    for (std::int32_t k = 0; k < voxels.counts.z; ++k)
      for (std::int32_t j = 0; j < voxels.counts.y; ++j)
        for (std::int32_t i = 0; i < voxels.counts.x; ++i)
        {
          std::size_t voxel = voxels.indexOf(i, j, k);
          const bool keep = isOpen(i, j, k);
          voxels.regionOf[voxel] = keep ? 0 : -1;
          voxels.clearance[voxel] = keep ? width : 0.0f;
          if (keep) ++open;
        }
    voxels.numberOfVoidVoxels = open;
  };

  // Hops are met from the bottom of the grid up, so a straw along x low down closes before a sheet higher up, and a
  // sheet low down before a straw above it. Both ways round the answer is the sheet's.
  carve([&](std::int32_t i, std::int32_t j, std::int32_t k) { return (j == 2 && k == 2) || (k == 20 && (i == 8 || j == 8)); });
  EXPECT_EQ(BruteForceDiameters::compute(structure, voxels).dimensionality, 2uz);

  carve([&](std::int32_t i, std::int32_t j, std::int32_t k) { return (j == 26 && k == 26) || (k == 4 && (i == 8 || j == 8)); });
  EXPECT_EQ(BruteForceDiameters::compute(structure, voxels).dimensionality, 2uz);

  // Two straws along different directions that never touch run one way each.
  carve([&](std::int32_t i, std::int32_t j, std::int32_t k) { return (i == 2 && j == 2) || (j == 16 && k == 12); });
  EXPECT_EQ(BruteForceDiameters::compute(structure, voxels).dimensionality, 1uz);
}


// The roomiest place along the widest way through is found by walking uphill from the roomiest few voxels of
// that part of the void, as the roomiest place anywhere is from the roomiest few of all of it. The room at the
// middle of a voxel is short of the room at the peak beside it by as much as the peak is off the middle, so
// which cavity has the roomiest voxel is settled by where the grid fell and not by which is the larger, wherever
// two are within a voxel of each other in size. Walked from the one roomiest voxel alone, the roomiest place on
// the path came out smaller than the largest sphere of the crystal in crystals where that sphere is on the path:
// by a quarter of an Ångström in the last of these.
TEST(brute_force_diameters, the_roomiest_place_on_the_path_is_walked_to_from_more_than_one_voxel)
{
  // Atoms of several sizes scattered in a cell, on a grid coarse enough for neighbouring cavities to be told
  // apart wrongly by their voxels. In each the largest cavity is on the widest way through.
  struct Case
  {
    std::size_t seed;
    double fromOneVoxel;  // what the walk from the roomiest voxel alone came to, as a diameter
  };
  for (const Case &entry : {Case{23, 6.311101}, Case{38, 6.331454}, Case{44, 6.417676}})
  {
    RandomNumber random{std::optional<std::size_t>(entry.seed)};
    const UnitCell cell(11.0, 12.0, 13.0);
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::size_t i = 0; i < 26; ++i)
    {
      positions.push_back(cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
      radii.push_back(0.9 + 0.9 * random.uniform());
    }
    BruteForceStructure structure = BruteForceStructure::create("scattered", cell, positions, radii);
    BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.4);
    const BruteForceDiameters diameters = BruteForceDiameters::compute(structure, voxels);

    ASSERT_TRUE(diameters.percolates) << "seed " << entry.seed;
    EXPECT_NEAR(diameters.includedAlongFreePathDiameter, diameters.includedSphereDiameter, 1.0e-9)
        << "seed " << entry.seed;
    EXPECT_GT(diameters.includedAlongFreePathDiameter, entry.fromOneVoxel + 0.05) << "seed " << entry.seed;

    // And it is where it is said to be: the room there is half the diameter.
    EXPECT_NEAR(2.0 * structure.clearance(diameters.includedAlongFreePathCentre),
                diameters.includedAlongFreePathDiameter, 1.0e-9)
        << "seed " << entry.seed;
  }
}
