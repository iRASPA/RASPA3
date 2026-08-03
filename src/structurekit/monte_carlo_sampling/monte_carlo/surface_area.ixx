module;

export module mc_surface_area;

import std;

import sampled_structure;

// The accessible surface area by the method of Düren et al.: points thrown at random over each atom's
// contact sphere, and the ones that fall inside another atom's thrown away.
//
// The share of a sphere that survives is the share of its area that is exposed, so the total is a sum of
// 4 pi R_i² over the atoms weighted by those shares. The error is statistical and falls off as one over the
// square root of the number of points, which is what the two counts below buy: `numberOfInnerSteps` points
// per atom in a pass, and `numberOfIterations` independent passes whose spread is what the estimate's own
// uncertainty could be read off.
export struct MC_SurfaceArea
{
  double surfaceArea{0.0};  // Å², averaged over the passes
  double seconds{0.0};

  void run(const SampledStructure &structure, const SampledProbe &probe,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps);
};
