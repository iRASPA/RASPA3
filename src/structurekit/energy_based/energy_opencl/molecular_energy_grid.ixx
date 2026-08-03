module;

export module energy_opencl_molecular_energy_grid;

import std;

import uint3;
import crystal;
import pair_interactions;

import energy_shared_linear_probe;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// Builds a `MolecularEnergyGrid` on the GPU, one work item per voxel, each turning the molecule through every
// orientation in the quadrature and reducing over them twice: once to the least energy, once to the free
// energy.
//
// Both reductions come out of a single pass over the orientations, because the expensive part is the sum over
// framework atoms and it is the same sum for both.
export struct MolecularEnergyGridOpenCL
{
  static MolecularEnergyGrid compute(const PairInteractions &interactions, const Crystal &framework,
                                     const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                     double temperature, const ElectrostaticPotentialGrid *potential = nullptr);

  // Which framework atom holds each point of the cell, when a molecule rather than a single site is asked.
  // See `molecularStrongestAtom` in the backend for what the label means and how it is arrived at.
  //
  // This is the one quantity here that needs the energy of every orientation at once rather than a running
  // reduction of them, since an orientation's weight is only known once the deepest is. Holding them all for
  // every voxel would be a gigabyte at 128³ over a hundred and twenty-eight orientations, so the cell is
  // taken a slab of z-planes at a time and the two kernels are run over each slab in turn.
  static std::vector<std::int32_t> strongestAtoms(const PairInteractions &interactions, const Crystal &framework,
                                                  const LinearProbe &probe, uint3 gridSize,
                                                  std::size_t numberOfOrientations, double temperature,
                                                  const ElectrostaticPotentialGrid *potential = nullptr);

  static const char *molecularKernelSource;
  static const char *tessellationKernelSource;
};
