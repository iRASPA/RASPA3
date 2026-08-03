module;

export module energy_probe_energy_grid;

import std;

import uint3;
import crystal;
import pair_interactions;

import energy_shared_probe_energy_grid;

// Builds a `ProbeEnergyGrid` on the processor, one thread per plane of voxels, each summing over every
// framework atom and every periodic image of it that reaches within the cutoff.
//
// It is the same sum the GPU builder makes, in the same order, with the same image range and the same
// truncation. What differs is the precision: this one works in double throughout, so where the two disagree
// it is the GPU that is rounding. That makes this the reference the faster route is held to, and the reason
// it is worth having at all despite being the slower of the two by a wide margin.
export struct ProbeEnergyGridCPU
{
  // Throws if the probe is not in the force field.
  static ProbeEnergyGrid compute(const PairInteractions &interactions, const Crystal &framework,
                                 std::string probePseudoAtom, uint3 gridSize);
};
