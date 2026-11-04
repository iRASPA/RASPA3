#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import opencl;
import opencl_clearance_grid;
import opencl_void_fraction;
import opencl_pore_analysis;
import opencl_pore_size_distribution;
import opencl_blocking_spheres;
import opencl_surface_area;

// What the clearance field is given, and what it makes of having nothing to measure from.

namespace
{

constexpr double edge = 12.0;
const uint3 gridSize(8, 8, 8);

PairInteractions probeOfRadiusOne()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 1;
  interactions.names = {"X"};
  interactions.charges = {0.0};
  interactions.parameters = {PairParameters{2.0, 100.0, 0.0}};
  return interactions;
}

Crystal box()
{
  Crystal framework;
  framework.name = "grid_clearance_test";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;
  return framework;
}

Crystal oneAtom()
{
  Crystal framework = box();
  const double3 fractional(0.25, 0.5, 0.75);
  framework.atoms.push_back(CrystalAtom{fractional * edge, 0.0, 0});
  framework.fractionalPositions.push_back(fractional);
  return framework;
}

struct ReportCleanup
{
  ~ReportCleanup()
  {
    for (const std::filesystem::directory_entry &entry : std::filesystem::directory_iterator("."))
    {
      if (entry.path().filename().string().starts_with("grid_clearance_test.")) std::filesystem::remove(entry.path());
    }
  }
};

}  // namespace


// With no atom in the box nothing limits the room at any point, so a probe of any size fits everywhere: the box is
// all void, not all solid. The field says so the way the kernel would, with the value it starts every point from.
TEST(grid_clearance, an_empty_box_has_room_everywhere)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  ReportCleanup cleanup;
  PairInteractions interactions = probeOfRadiusOne();
  ClearanceGrid grid = ClearanceGrid::compute(interactions, box(), gridSize);

  ASSERT_EQ(grid.numberOfVoxels(), 512uz);
  for (std::size_t voxel = 0; voxel < grid.numberOfVoxels(); ++voxel)
  {
    EXPECT_EQ(grid.clearance[voxel], std::numeric_limits<float>::max()) << "at " << voxel;
    EXPECT_EQ(grid.closestAtom[voxel], -1) << "at " << voxel;
  }

  GridVoidFraction fraction;
  fraction.run(interactions, box(), "X", grid);
  EXPECT_DOUBLE_EQ(fraction.voidFraction, 1.0);
  EXPECT_DOUBLE_EQ(fraction.accessibleVolumeFraction, 1.0);
  EXPECT_DOUBLE_EQ(fraction.inaccessibleVolumeFraction, 0.0);
  EXPECT_EQ(fraction.numberOfChannels, 1uz);
}

// Fractional positions are carried for convenience and a crystal may come without them; the other routes work
// them out from the atoms when it does. Such a crystal is not an empty box.
TEST(grid_clearance, atoms_without_fractional_positions_are_still_atoms)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeOfRadiusOne();
  Crystal carried = oneAtom();
  Crystal derived = oneAtom();
  derived.fractionalPositions.clear();

  ClearanceGrid expected = ClearanceGrid::compute(interactions, carried, gridSize);
  ClearanceGrid grid = ClearanceGrid::compute(interactions, derived, gridSize);

  ASSERT_EQ(grid.numberOfVoxels(), expected.numberOfVoxels());
  for (std::size_t voxel = 0; voxel < grid.numberOfVoxels(); ++voxel)
  {
    EXPECT_FLOAT_EQ(grid.clearance[voxel], expected.clearance[voxel]) << "at " << voxel;
    EXPECT_EQ(grid.closestAtom[voxel], 0) << "at " << voxel;
  }
}

// Two lists that are meant to be the same atoms and are not the same length have no right answer: one way the radius
// is read past the end of its list, and the other way atoms are silently left out of the field.
TEST(grid_clearance, positions_and_atoms_of_different_lengths_are_refused)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeOfRadiusOne();

  Crystal morePositions = oneAtom();
  morePositions.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  EXPECT_THROW(ClearanceGrid::compute(interactions, morePositions, gridSize), std::runtime_error);

  Crystal moreAtoms = oneAtom();
  moreAtoms.atoms.push_back(CrystalAtom{double3(6.0, 6.0, 6.0), 0.0, 0});
  EXPECT_THROW(ClearanceGrid::compute(interactions, moreAtoms, gridSize), std::runtime_error);
}

// What the rest of the grid route makes of that field: one channel running every way with no surface and no pocket,
// and all of the void in pores larger than any size asked about.
TEST(grid_clearance, the_analyses_of_an_empty_box_agree_that_it_is_all_room)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  ReportCleanup cleanup;
  PairInteractions interactions = probeOfRadiusOne();
  Crystal framework = box();
  framework.mass = 0.0;
  ClearanceGrid grid = ClearanceGrid::compute(interactions, framework, gridSize);
  const double unlimited = 2.0 * static_cast<double>(std::numeric_limits<float>::max());

  GridPoreDiameters diameters = GridPoreDiameters::compute(grid);
  EXPECT_TRUE(diameters.percolates);
  EXPECT_EQ(diameters.dimensionalityAtThreshold, 3);
  EXPECT_DOUBLE_EQ(diameters.includedSphereDiameter, unlimited);
  EXPECT_DOUBLE_EQ(diameters.freeSphereDiameter, unlimited);
  EXPECT_DOUBLE_EQ(diameters.includedAlongFreePathDiameter, unlimited);

  GridPoreAnalysis analysis;
  analysis.run(interactions, framework, "X", grid);
  EXPECT_EQ(analysis.channels.numberOfChannels, 1uz);
  EXPECT_EQ(analysis.channels.numberOfPockets, 0uz);

  GridBlockingSpheres blocking;
  blocking.run(interactions, framework, "X", grid);
  EXPECT_EQ(blocking.numberOfPockets, 0uz);

  GridSurfaceArea surface;
  surface.run(interactions, framework, "X", grid);
  EXPECT_DOUBLE_EQ(surface.totalSurfaceArea, 0.0);

  GridPoreSizeDistribution distribution;
  distribution.run(interactions, framework, "X", grid, 20.0, 10uz);
  EXPECT_DOUBLE_EQ(distribution.voidFraction, 1.0);
  EXPECT_DOUBLE_EQ(distribution.accessibleVoidFraction, 1.0);
  EXPECT_DOUBLE_EQ(distribution.largestDiameter, unlimited);
  ASSERT_EQ(distribution.points.size(), 10uz);
  for (const GridPoreSizeDistributionPoint &point : distribution.points)
  {
    EXPECT_DOUBLE_EQ(point.cumulativeVoidFraction, 1.0) << "at " << point.diameter;
    EXPECT_DOUBLE_EQ(point.cumulativeAccessibleFraction, 1.0) << "at " << point.diameter;
    EXPECT_DOUBLE_EQ(point.distribution, 0.0) << "at " << point.diameter;
  }
}
