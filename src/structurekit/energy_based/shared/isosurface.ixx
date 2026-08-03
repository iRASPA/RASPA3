module;

export module energy_shared_isosurface;

import std;

import uint3;
import double3;
import double3x3;
import crystal;

// The area of the surface where a field crosses a value, and how many triangles it took.
export struct IsosurfaceArea
{
  double area{0.0};  // [Å²] per unit cell
  std::size_t numberOfTriangles{0};
  std::size_t numberOfRejectedTriangles{0};
};

// A marching-cubes triangle is confined to a single voxel, so it cannot be much larger than the biggest voxel
// face. Anything beyond that is numerical debris rather than surface, and the bound has to track the grid
// spacing: a fixed cut-off would start discarding genuine triangles on coarse grids, which on a large cell
// can throw away most of the area.
export double largestPlausibleTriangleArea(const double3x3 &unitCell, uint3 gridSize);

// Sums the areas of triangles given by their corners in fractional coordinates, three corners to a triangle,
// discarding those too large to have come from one voxel. Both extractors end here, so a difference between
// them is a difference in the triangles they found and not in how the triangles were added up.
export IsosurfaceArea accumulateTriangleAreas(const double3x3 &unitCell, uint3 gridSize,
                                              std::span<const double3> corners);
