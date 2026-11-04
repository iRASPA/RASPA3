module;

export module integration_surface_area;

import std;

import sampled_structure;

// The same accessible surface area as the Monte Carlo route, with the random points replaced by a regular
// grid of latitudes and longitudes over each atom's contact sphere.
//
// Each point of that grid stands for an equal share of the sphere: the stacks are equal steps of height along the
// axis, which by Archimedes are equal steps of area, each taken at its middle, and the points of a stack are equally
// spaced round it. The answer is still only as fine as the grid, since a point is either buried or not. What it is
// besides is deterministic: the same structure gives the same number twice, without the scatter of the random
// points, which is what makes the two routes worth running against each other.
export struct Integration_SurfaceArea
{
  double surfaceArea{0.0};  // Å²
  double seconds{0.0};

  void run(const SampledStructure &structure, const SampledProbe &probe, std::optional<std::size_t> numberOfSlices);
};
