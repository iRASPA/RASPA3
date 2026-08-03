module;

module mc_opencl_blocking_pockets;

import std;

import sampled_structure;
import sampled_roadmap;
import sampling_backend;
import mc_opencl_backend;
import mc_blocking_pockets;

void MC_OpenCL_BlockingPockets::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                               std::optional<std::size_t> numberOfInnerSteps)
{
  MC_BlockingPockets::run(
      structure, SampledRoadmap::build(structure, samplingBackendOpenCL(), numberOfIterations, numberOfInnerSteps));
}
