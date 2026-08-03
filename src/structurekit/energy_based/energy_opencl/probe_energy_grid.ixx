module;

export module energy_opencl_probe_energy_grid;

import std;

import uint3;
import crystal;
import pair_interactions;

import energy_shared_probe_energy_grid;

// Builds a `ProbeEnergyGrid` on the GPU, one work item per voxel, each summing over every framework atom and
// every periodic image of it that reaches within the cutoff.
//
// The field this returns is the same object the CPU builder returns, and the two are meant to agree to within
// the difference between single and double precision. Everything downstream takes the field and cannot tell
// which of the two produced it.
export struct ProbeEnergyGridOpenCL
{
  // Throws if no OpenCL device was found, or if the probe is not in the force field.
  static ProbeEnergyGrid compute(const PairInteractions &interactions, const Crystal &framework,
                                 std::string probePseudoAtom, uint3 gridSize);

  static const char *energyKernelSource;
};
