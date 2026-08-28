module;

export module energy_shared_isosurface;

import std;

import uint3;
import double3;
import double3x3;
import crystal;
import surface_curvature;

// The area of the surface where a field crosses a value, and how many triangles it took.
export struct IsosurfaceArea
{
  double area{0.0};  // [Å²] per unit cell
  std::size_t numberOfTriangles{0};
  std::size_t numberOfRejectedTriangles{0};

  // How that area divides by the shape of the surface, when the extractor handed back the gradients it stored
  // at each vertex as well as the vertices themselves. Left at zero when it did not.
  //
  // On an energy field this is a real property of the surface rather than of the grid. The field is a sum over
  // the atoms and so is smooth everywhere: the walls of neighbouring atoms blend into one another and the
  // saddle between two of them is a genuine saddle whose curvature does not move as the grid is refined. That
  // is not true of the clearance field, which is a minimum over the atoms and has a crease wherever two of them
  // are equally near.
  CurvatureAreas curvature;
  CurvatureBands curvatureBands;
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

// The same, and the curvature split besides, from the gradients the extractor stored at the same vertices, one
// to a corner and held per grid step. `sense` says which way the field grows relative to the void, which is
// what turns a gradient into an outward normal; for an energy field, larger meaning deeper into the wall, that
// is `GrowsIntoSolid`.
//
// The split is taken over exactly the triangles the area was taken over, the implausibly large ones being
// dropped from both, so the two always describe the same surface.
export IsosurfaceArea accumulateTriangleAreas(const double3x3 &unitCell, uint3 gridSize,
                                              std::span<const double3> corners,
                                              std::span<const double3> gradients, FieldSense sense);

// Writes the curvature split of an energy iso-surface, with the commentary that says how to read it. Does
// nothing when the extractor handed back no gradients. `name` labels the row, there being routes that report
// more than one surface.
export void writeIsosurfaceCurvature(std::ostream &stream, const IsosurfaceArea &surface, const char *name);
