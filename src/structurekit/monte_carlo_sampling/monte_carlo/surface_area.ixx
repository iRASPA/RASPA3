module;

export module mc_surface_area;

import std;

import sampled_structure;

// The accessible surface area by the method of Düren et al.: points thrown at random over each atom's
// contact sphere, and the ones that fall inside another atom's thrown away.
//
// The share of a sphere that survives is the share of its area that is exposed, so the total is a sum of
// 4 pi R_i² over the atoms weighted by those shares. Each pass is one independent reading of that total;
// the running sum and sum of squares of those readings give the mean and a 95% confidence interval.
export struct MC_SurfaceArea
{
  double surfaceArea{0.0};       // Å², averaged over the passes
  double surfaceAreaError{0.0};  // Å², half-width of the 95% confidence interval of the mean
  double seconds{0.0};

  void run(const SampledStructure &structure, const SampledProbe &probe,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps);
};
