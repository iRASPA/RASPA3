module;

module mc_opencl_channels;

import std;

import sampled_structure;
import sampled_roadmap;
import sampling_backend;
import mc_opencl_backend;
import mc_channels;

void MC_OpenCL_Channels::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                        std::optional<std::size_t> numberOfInnerSteps)
{
  MC_Channels::run(
      structure, SampledRoadmap::build(structure, samplingBackendOpenCL(), numberOfIterations, numberOfInnerSteps));
}
