module;

module mc_opencl_window_shape;

import std;

import sampled_structure;
import sampled_roadmap;
import sampling_backend;
import mc_opencl_backend;
import mc_window_shape;

void MC_OpenCL_WindowShape::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                           std::optional<std::size_t> numberOfInnerSteps)
{
  MC_WindowShape::run(
      structure, SampledRoadmap::build(structure, samplingBackendOpenCL(), numberOfIterations, numberOfInnerSteps));
}
