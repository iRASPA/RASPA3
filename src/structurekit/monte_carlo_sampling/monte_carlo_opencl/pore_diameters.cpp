module;

module mc_opencl_pore_diameters;

import std;

import sampled_structure;
import sampled_roadmap;
import sampling_backend;
import mc_opencl_backend;
import mc_pore_diameters;

void MC_OpenCL_PoreDiameters::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                             std::optional<std::size_t> numberOfInnerSteps)
{
  MC_PoreDiameters::run(
      structure, SampledRoadmap::build(structure, samplingBackendOpenCL(), numberOfIterations, numberOfInnerSteps));
}
