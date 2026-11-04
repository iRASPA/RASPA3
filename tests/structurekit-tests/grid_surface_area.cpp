#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import opencl_clearance_grid;
import opencl_surface_area;
import grid_connected_components;

// How the clearance-grid route shares out its surface: among the pores it faces and the atoms it belongs to.
//
// The field is written out by hand rather than computed, which is what lets these run without a GPU and what
// lets a voxel name whichever atom a case calls for. The surface is the level set at the probe's radius of one.

namespace
{

constexpr double edge = 10.0;
constexpr std::size_t points = 8;
constexpr float solid = 0.5f;

PairInteractions probeOfRadiusOne()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 1;
  interactions.names = {"X"};
  interactions.charges = {0.0};
  interactions.parameters = {PairParameters{2.0, 100.0, 0.0}};
  return interactions;
}

Crystal twoAtoms()
{
  Crystal framework;
  framework.name = "grid_surface_area_test";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 32.0;
  for (double3 fractional : {double3(0.0, 0.0, 0.1), double3(0.0, 0.0, 0.7)})
  {
    framework.atoms.push_back(CrystalAtom{fractional * edge, 0.0, 0});
    framework.fractionalPositions.push_back(fractional);
  }
  return framework;
}

// All solid, and all of it nearest the first atom, until a case paints something else in.
ClearanceGrid solidGrid()
{
  ClearanceGrid grid;
  grid.gridSize = uint3(points, points, points);
  grid.unitCell = UnitCell(edge, edge, edge);
  grid.clearance.assign(points * points * points, solid);
  grid.closestAtom.assign(points * points * points, 0);
  grid.atomRadii = {1.0, 1.0};
  return grid;
}

// Four channels along x around the point (j, k) = (4, 4), each diagonal to it, and that point itself open but
// touching none of them through a face. It only just clears the probe, the channels clear it by a long way,
// and every cube of the grid that holds the speck also holds a channel.
ClearanceGrid speckAmongChannels()
{
  ClearanceGrid grid = solidGrid();
  for (std::size_t i = 0; i < points; ++i)
  {
    for (std::size_t j : {std::size_t{3}, std::size_t{5}})
    {
      for (std::size_t k : {std::size_t{3}, std::size_t{5}})
      {
        grid.clearance[grid.voxelIndex(i, j, k)] = 3.0f;
      }
    }
  }
  grid.clearance[grid.voxelIndex(4, 4, 4)] = 1.5f;
  return grid;
}

// Removes the report a run leaves behind, whichever way the test ends.
struct ReportCleanup
{
  ~ReportCleanup()
  {
    std::error_code ignored;
    std::filesystem::remove("grid_surface_area_test.grid.sa.gpu.txt", ignored);
  }
};

}  // namespace


// The speck comes first in the order the corners are visited, for the cube whose lowest corner it is, and a
// search that settles for the first open corner hands that cube's surface to a pocket of one voxel.
TEST(grid_surface_area, a_cube_holding_a_speck_and_a_channel_faces_the_channel)
{
  ReportCleanup cleanup;
  ClearanceGrid grid = speckAmongChannels();

  GridComponents components = GridComponents::compute(grid.gridSize, grid.clearance, 1.0);
  ASSERT_EQ(components.numberOfChannels, 4u);
  ASSERT_EQ(components.numberOfPockets, 1u);

  GridSurfaceArea surface;
  surface.run(probeOfRadiusOne(), twoAtoms(), "X", grid);

  EXPECT_GT(surface.totalSurfaceArea, 0.0);
  EXPECT_DOUBLE_EQ(surface.inaccessibleSurfaceArea, 0.0);
  EXPECT_DOUBLE_EQ(surface.undecidedSurfaceArea, 0.0);
  EXPECT_DOUBLE_EQ(surface.accessibleSurfaceArea, surface.totalSurfaceArea);
}


// A slab of void one voxel thick, whose only open plane lies just past the sheet midway between the walls and
// so is nearest the second atom. Each wall is a flat sheet across the cell and belongs to the atom behind it,
// where an answer read off the open plane gives both of them to the second atom.
TEST(grid_surface_area, a_wall_of_a_narrow_neck_belongs_to_the_atom_behind_it)
{
  ReportCleanup cleanup;
  ClearanceGrid grid = solidGrid();
  for (std::size_t k = 3; k < points; ++k)
  {
    for (std::size_t j = 0; j < points; ++j)
    {
      for (std::size_t i = 0; i < points; ++i)
      {
        grid.closestAtom[grid.voxelIndex(i, j, k)] = 1;
        if (k == 3) grid.clearance[grid.voxelIndex(i, j, k)] = 1.5f;
      }
    }
  }

  GridSurfaceArea surface;
  surface.run(probeOfRadiusOne(), twoAtoms(), "X", grid);

  ASSERT_EQ(surface.atomSurfaceArea.size(), 2u);
  EXPECT_NEAR(surface.totalSurfaceArea, 2.0 * edge * edge, 1.0e-9);
  EXPECT_NEAR(surface.atomSurfaceArea[0], edge * edge, 1.0e-9);
  EXPECT_NEAR(surface.atomSurfaceArea[1], edge * edge, 1.0e-9);
}


TEST(grid_surface_area, running_twice_gives_the_same_answer_as_running_once)
{
  ReportCleanup cleanup;
  ClearanceGrid grid = speckAmongChannels();
  PairInteractions interactions = probeOfRadiusOne();
  Crystal framework = twoAtoms();

  GridSurfaceArea once;
  once.run(interactions, framework, "X", grid);

  GridSurfaceArea twice;
  twice.run(interactions, framework, "X", grid);
  twice.run(interactions, framework, "X", grid);

  ASSERT_GT(once.totalSurfaceArea, 0.0);
  EXPECT_DOUBLE_EQ(twice.totalSurfaceArea, once.totalSurfaceArea);
  EXPECT_DOUBLE_EQ(twice.accessibleSurfaceArea, once.accessibleSurfaceArea);
  EXPECT_DOUBLE_EQ(twice.inaccessibleSurfaceArea, once.inaccessibleSurfaceArea);
  EXPECT_DOUBLE_EQ(twice.undecidedSurfaceArea, once.undecidedSurfaceArea);
  EXPECT_DOUBLE_EQ(twice.curvature.total(), once.curvature.total());
  EXPECT_DOUBLE_EQ(twice.accessibleCurvature.total(), once.accessibleCurvature.total());
  EXPECT_DOUBLE_EQ(twice.inaccessibleCurvature.total(), once.inaccessibleCurvature.total());
}
