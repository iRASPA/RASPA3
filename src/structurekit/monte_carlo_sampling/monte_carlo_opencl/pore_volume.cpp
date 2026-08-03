module;

module mc_opencl_pore_volume;

import std;

import sampled_structure;
import sampled_roadmap;
import sampling_backend;
import mc_opencl_backend;
import mc_pore_volume;

void MC_OpenCL_PoreVolume::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                          std::optional<std::size_t> numberOfInnerSteps)
{
  MC_PoreVolume::run(
      structure, SampledRoadmap::build(structure, samplingBackendOpenCL(), numberOfIterations, numberOfInnerSteps));
}
