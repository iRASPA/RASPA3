module;

export module energy_electrostatic_potential_grid;

import std;

import uint3;
import crystal;
import pair_interactions;

import energy_shared_electrostatic_potential_grid;

// Builds the far part of the Ewald potential on the processor, one thread per plane of voxels, each summing
// over the wave vectors that survived the cutoff in reciprocal space.
//
// The near part is not tabulated and is not built here; whatever reads the field sums it directly. See the
// field's own description for why the split falls where it does.
export struct ElectrostaticPotentialGridCPU
{
  static ElectrostaticPotentialGrid compute(const PairInteractions &interactions, const Crystal &framework, uint3 gridSize,
                                            double relativePrecision, std::optional<double> alphaOverride = {});
};
