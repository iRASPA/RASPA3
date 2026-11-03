#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
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
