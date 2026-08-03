module;

export module energy_molecular_energy_grid;

import std;

import uint3;
import crystal;
import pair_interactions;

import energy_shared_linear_probe;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// Builds a `MolecularEnergyGrid` on the processor, one thread per plane of voxels, each turning the molecule
// through every orientation in the quadrature and reducing over them twice: once to the least energy, once to
// the free energy.
//
// Both reductions come out of a single pass over the orientations, because the expensive part is the sum over
// framework atoms and it is the same sum for both. This is the same work the GPU builder does, in double
// precision, and it costs the number of orientations times what the single-site field costs.
export struct MolecularEnergyGridCPU
{
  static MolecularEnergyGrid compute(const PairInteractions &interactions, const Crystal &framework,
                                     const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                     double temperature, const ElectrostaticPotentialGrid *potential = nullptr);

  // Which framework atom holds each point of the cell, when a molecule rather than a single site is asked.
  // See `molecularStrongestAtom` in the backend for what the label means and how it is arrived at.
  static std::vector<std::int32_t> strongestAtoms(const PairInteractions &interactions, const Crystal &framework,
                                                  const LinearProbe &probe, uint3 gridSize,
                                                  std::size_t numberOfOrientations, double temperature,
                                                  const ElectrostaticPotentialGrid *potential = nullptr);
};
