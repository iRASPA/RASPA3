module;

export module opencl_surface_area;

import std;

import uint3;
import double3;
import framework;
import forcefield;
import opencl_clearance_grid;
import opencl_connected_components;

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

  // Area attributed to each atom of the unit cell, in the order of `framework.unitCellAtoms`.
  std::vector<double> atomSurfaceArea;

  std::size_t numberOfTriangles{0};
  double seconds{0.0};

  GridSurfaceArea();
  ~GridSurfaceArea();

  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom, uint3 gridSize);
  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid);
};
