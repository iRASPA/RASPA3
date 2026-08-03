module;

export module energy_isosurface;

import std;

import uint3;
import double3;
import crystal;

import energy_shared_isosurface;

// Extracts the surface where a field crosses a value on the processor, by marching cubes.
//
// Which field is handed in decides what the area means: the energy field of a single probe atom gives the
// surface that atom sees, a molecular free-energy field gives the surface a whole molecule sees once its
// orientations are averaged over.
//
// The cell is periodic and both extractors treat it so, counting every cube of it exactly once including the
// ones that span its far face. The GPU reaches round with a modulo on the sample index; this one repeats the
// first plane at the far end instead, which comes to the same thing.
export struct EnergyIsosurface
{
  // The triangles themselves, three corners to a triangle, in fractional coordinates. Handing these back is
  // what lets a surface be divided among the atoms rather than only measured.
  static std::vector<double3> trianglesOfIsosurface(std::span<const float> field, uint3 gridSize, double isoValue);

  static IsosurfaceArea areaOfIsosurface(const Crystal &framework, std::span<const float> field, uint3 gridSize,
                                         double isoValue);
};
