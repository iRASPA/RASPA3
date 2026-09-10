module;

export module energy_opencl_lewiner_isosurface;

import std;

import uint3;
import double3;

// GPU marching cubes using Lewiner's tables and the same face / interior tests as the
// processor extractor. The older histo-pyramid path used Bourke's 256-case table (at most
// five triangles, no ambiguity tests) and a power-of-two image pyramid capped at 512³.
//
// This one walks the real grid, including the cubes that wrap the periodic faces, and
// writes the same edge triples the CPU `MarchingCubes` would. Positions are still
// interpolated in float, so they will not be bit-identical, but the mesh is the Lewiner
// mesh: the same cubes, the same tilings, the same interior vertices.
export std::vector<double3> trianglesOfLewinerIsosurface(std::span<const float> field, uint3 gridSize,
                                                         double isoValue,
                                                         std::vector<double3> *gradients = nullptr);
