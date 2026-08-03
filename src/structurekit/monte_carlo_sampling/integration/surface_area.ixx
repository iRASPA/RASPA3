module;

export module integration_surface_area;

import std;

import sampled_structure;

// The same accessible surface area as the Monte Carlo route, with the random points replaced by a regular
// grid of latitudes and longitudes over each atom's contact sphere.
//
// The points of that grid are not spread evenly over the sphere -- they crowd towards the poles -- so this
// is not a quadrature and the answer is not exact. What it is instead is deterministic: the same structure
// gives the same number twice, and refining `numberOfSlices` moves it in one direction rather than jittering
// it, which is what makes the two routes worth running against each other.
export struct Integration_SurfaceArea
{
  double surfaceArea{0.0};  // Å²
  double seconds{0.0};

  void run(const SampledStructure &structure, const SampledProbe &probe, std::optional<std::size_t> numberOfSlices);
};
