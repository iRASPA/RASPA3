module;

export module energy_well_field;

import std;

import uint3;
import double3;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_well_field;

// The processor's well field and the vertex refinement that slides the extracted sheet onto the analytic
// well floor. Both are the same arithmetic the GPU builder does, in double precision.
export struct WellFieldCPU
{
  static WellField compute(const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                           uint3 gridSize, std::size_t numberOfOrientations, double thermalEnergy,
                           std::span<const BlockingSphere> blockingSpheres, double blockedEnergyPerAngstrom,
                           double ceiling, const ElectrostaticPotentialGrid *potential = nullptr,
                           double coulombFactor = 0.0);

  static void refineVertices(std::vector<double3> &corners, std::vector<double> &energies,
                             const PairInteractions &interactions, const Crystal &framework,
                             const NeighbourhoodParameters &parameters,
                             std::span<const BlockingSphere> blockingSpheres, double iso);
};
