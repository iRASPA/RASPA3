module;

export module opencl_clearance_grid;

import std;

import int3;
import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;

// The clearance field of a framework, sampled on a regular grid over the unit cell.
//
// `clearance` holds, at each grid point, the room there is at that point for the centre of a hard sphere:
// the distance from the point to the nearest atom surface, min_i(|x - p_i| - r_i), with r_i half the atom's
// Lennard-Jones size parameter, which is the radius the Voronoi and Apollonius routes measure with. The
// value is negative at a point inside an atom.
//
// A probe of radius r fits at a point exactly where the clearance is at least r, so one field answers every
// probe: everything downstream only compares the field against a number. That is what makes the grid route
// cheap for a sweep over probe sizes, or for a search over probe sizes such as the free-sphere diameter.
//
// `closestAtom` names the atom whose surface is the nearest one. It is the discrete nearest-surface
// (Apollonius) tessellation of the cell, and it is what divides an area or a volume among individual atoms.
export struct ClearanceGrid
{
  uint3 gridSize{0, 0, 0};
  UnitCell unitCell;

  // Both fields are indexed by `voxelIndex`, x varying fastest.
  std::vector<float> clearance;
  std::vector<std::int32_t> closestAtom;

  // The radii the field was built from, in the order of `framework.atoms`.
  std::vector<double> atomRadii;

  double seconds{0.0};

  ClearanceGrid();
  ~ClearanceGrid();

  std::size_t numberOfVoxels() const { return this->clearance.size(); }

  std::size_t voxelIndex(std::size_t i, std::size_t j, std::size_t k) const
  {
    return (k * this->gridSize.y + j) * this->gridSize.x + i;
  }

  double3 fractionalPosition(std::size_t i, std::size_t j, std::size_t k) const;
  double3 cartesianPosition(std::size_t i, std::size_t j, std::size_t k) const;

  // The cell volume one grid point stands for. The sampling is endpoint-exclusive, so every point carries
  // the same share and the shares add up to the cell.
  double voxelVolume() const;

  // The largest sphere that fits anywhere in the void has this radius, to within the grid spacing.
  double maximumClearance() const;

  // The grid spacing along each cell axis, as a distance rather than a fraction.
  double3 spacing() const;

  // Writes out the tessellation: which atom's surface is the nearest one at every grid point, and how far away
  // it is. The cells are the nearest-surface (Apollonius) cells of the atoms, since that is what the field is a
  // distance to, and they are the same cells the volume and area of this route are divided along.
  void writeTessellation(const Crystal &framework) const;

  // Builds the field on the GPU. Throws if no OpenCL device was found.
  static ClearanceGrid compute(const PairInteractions &interactions, const Crystal &framework, uint3 gridSize);

  static const char *clearanceKernelSource;
};
