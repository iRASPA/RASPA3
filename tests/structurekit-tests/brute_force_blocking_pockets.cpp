#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import brute_force_structure;
import brute_force_voxels;
import brute_force_blocking_pockets;

// The blocking spheres the flooded voxels ask for. A voxel stands for a little box, so a pocket's covering
// radius has to reach the far corner of every box of the pocket, and its channel radius may reach into no box
// of a channel; what is checked here is that both hold against the boxes themselves.

namespace
{

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


TEST(brute_force_blocking_pockets, the_radii_are_measured_to_the_corners_of_slanted_voxels)
{
  // A cage in a hexagonal cell, so tight that one voxel centre is all the room it has. The pocket is then that
  // one voxel, its centre is the voxel's centre, and its covering radius is exactly half the voxel's longest
  // diagonal. On a cell with an obtuse angle that is further than the three widths across the voxel in
  // quadrature, which left the voxel's far corners outside the sphere.
  const double right = 0.5 * std::numbers::pi;
  const double spacing = 0.5;
  UnitCell box(16.0, 16.0, 16.0, right, right, 2.0 * std::numbers::pi / 3.0);

  // The grid the voxels will be built on, so that the cage can be centred on one of its voxels.
  double3 widths = box.perpendicularWidths();
  int3 counts(static_cast<std::int32_t>(std::round(widths.x / spacing)),
              static_cast<std::int32_t>(std::round(widths.y / spacing)),
              static_cast<std::int32_t>(std::round(widths.z / spacing)));
  double3 middleVoxel = box.cell * double3((0.5 * counts.x + 0.5) / counts.x, (0.5 * counts.y + 0.5) / counts.y,
                                           (0.5 * counts.z + 0.5) / counts.z);

  std::vector<double3> positions = fibonacciShell(middleVoxel, 1.6, 64);
  std::vector<double> radii(positions.size(), 1.5);
  BruteForceStructure structure = BruteForceStructure::create("tight cage", box, positions, radii);
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, spacing);
  ASSERT_EQ(voxels.counts.x, counts.x);
  ASSERT_EQ(voxels.counts.y, counts.y);
  ASSERT_EQ(voxels.counts.z, counts.z);

  BruteForceBlockingPockets blocking = BruteForceBlockingPockets::compute(structure, voxels);
  ASSERT_EQ(blocking.pockets.size(), 1uz);
  const BruteForcePocket &pocket = blocking.pockets.front();
  ASSERT_EQ(pocket.numberOfVoxels, 1uz);
  ASSERT_TRUE(pocket.hasChannel());
  EXPECT_NEAR((pocket.centre - middleVoxel).length(), 0.0, 1.0e-9);

  // The corners of the voxel boxes: the far corner of the pocket's own voxel bounds the covering radius from
  // below, and the near corner of any channel voxel bounds the channel radius from above.
  auto corners = [&](std::size_t voxel)
  {
    std::vector<double3> points;
    double3 centre = voxels.centre(structure, voxel);
    for (double sx : {-0.5, 0.5})
      for (double sy : {-0.5, 0.5})
        for (double sz : {-0.5, 0.5})
          points.push_back(centre + box.cell * double3(sx / static_cast<double>(counts.x),
                                                       sy / static_cast<double>(counts.y),
                                                       sz / static_cast<double>(counts.z)));
    return points;
  };

  double farthestOwnCorner = 0.0;
  double nearestChannelCorner = std::numeric_limits<double>::max();
  for (std::size_t voxel = 0; voxel < voxels.regionOf.size(); ++voxel)
  {
    const std::int32_t region = voxels.regionOf[voxel];
    if (region < 0) continue;
    const bool channel = voxels.regions[static_cast<std::size_t>(region)].percolates;
    for (const double3 &corner : corners(voxel))
    {
      double distance = structure.nearestImage(pocket.centre, corner).length();
      if (channel)
        nearestChannelCorner = std::min(nearestChannelCorner, distance);
      else
        farthestOwnCorner = std::max(farthestOwnCorner, distance);
    }
  }

  // The long diagonal of a voxel runs along a less b, which on this cell is root three times the step.
  double3 longDiagonal = box.cell * double3(1.0 / counts.x, -1.0 / counts.y, 1.0 / counts.z);
  EXPECT_NEAR(farthestOwnCorner, 0.5 * longDiagonal.length(), 1.0e-9);
  EXPECT_NEAR(pocket.coveringRadius, farthestOwnCorner, 1.0e-9);
  EXPECT_LE(pocket.channelRadius, nearestChannelCorner + 1.0e-9);
}
