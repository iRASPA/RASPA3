module;

export module opencl_surface_area;

import std;

import uint3;
import double3;
import crystal;
import pair_interactions;
import surface_curvature;
import opencl_clearance_grid;
import grid_connected_components;

// The surface a probe of a given radius rolls over, taken off the clearance field.
//
// The set of points the probe's centre can occupy is bounded by the surface where the clearance equals the
// probe's radius, and that surface is the accessible surface: the atoms grown by the probe's radius, which is
// the same surface the Voronoi and Apollonius routes measure. So there is nothing to sample and nothing to
// approximate about which surface is meant, only where a level set of a field on a grid lies, and marching
// cubes places that to second order in the spacing.
//
// Each piece of the surface belongs to whichever pore it faces, so the same division into channels and pockets
// that splits the void splits the area, and it belongs to whichever atom's surface is the nearest, which the
// field already records, so the area also divides among the atoms without any further work.
export struct GridSurfaceArea
{
  double probeRadius{0.0};

  double accessibleSurfaceArea{0.0};
  double inaccessibleSurfaceArea{0.0};
  double undecidedSurfaceArea{0.0};
  double totalSurfaceArea{0.0};

  // Area attributed to each atom of the unit cell, in the order of `framework.atoms`.
  std::vector<double> atomSurfaceArea;

  // The same area divided by the shape of the sheet rather than by what lies behind it. Marching cubes stores
  // the field's gradient at every vertex it places, and three of those to a triangle are enough to estimate
  // the two principal curvatures there, so the split costs a three-by-three solve per triangle and no extra
  // pass over the field. It says how much of the wall bulges into the void and how much of it is the inside of
  // a pocket, which is a different question from the one the channel/pocket split above answers.
  CurvatureAreas curvature;
  CurvatureAreas accessibleCurvature;
  CurvatureAreas inaccessibleCurvature;

  // Where the classification stops believing itself. The upper bound is set from the grid spacing, a curvature
  // whose radius is finer than a voxel being the crease of the field rather than a shape the grid resolved.
  CurvatureBands curvatureBands;

  std::size_t numberOfTriangles{0};
  double seconds{0.0};

  GridSurfaceArea();
  ~GridSurfaceArea();

  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom, uint3 gridSize);
  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid);
};
