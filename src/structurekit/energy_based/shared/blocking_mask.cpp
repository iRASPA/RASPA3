module;

module energy_shared_blocking_mask;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import voronoi_blocking_spheres;
import energy_shared_energy_backend;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import structure_parallel;


double blockingSphereDistance(double3 fractional, const UnitCell &unitCell, std::span<const BlockingSphere> spheres)
{
  double nearest = 1.0e10;
  for (const BlockingSphere &sphere : spheres)
  {
    double3 dr = fractional - sphere.centerFractional;
    dr.x -= std::rint(dr.x);
    dr.y -= std::rint(dr.y);
    dr.z -= std::rint(dr.z);
    nearest = std::min(nearest, (unitCell.cell * dr).length() - sphere.radius);
  }
  return nearest;
}


void applyBlockingSpheresToEnergy(std::span<float> energy, uint3 gridSize, const UnitCell &unitCell,
                                  std::span<const BlockingSphere> spheres, double energyPerAngstrom, double ceiling)
{
  if (spheres.empty() || energy.empty()) return;
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0) return;

  const std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  if (energy.size() < numberOfVoxels) return;

  std::vector<BlockingSphere> held(spheres.begin(), spheres.end());
  const double invX = 1.0 / static_cast<double>(gridSize.x);
  const double invY = 1.0 / static_cast<double>(gridSize.y);
  const double invZ = 1.0 / static_cast<double>(gridSize.z);

  forEachBlock(gridSize.z, workersAvailable(),
               [&](std::size_t, std::size_t begin, std::size_t end)
               {
                 for (std::size_t iz = begin; iz < end; ++iz)
                 {
                   for (std::size_t iy = 0; iy < gridSize.y; ++iy)
                   {
                     for (std::size_t ix = 0; ix < gridSize.x; ++ix)
                     {
                       const double3 fractional(static_cast<double>(ix) * invX, static_cast<double>(iy) * invY,
                                                static_cast<double>(iz) * invZ);
                       const double pocket = blockingSphereDistance(fractional, unitCell, held);
                       if (pocket < 0.0)
                       {
                         const std::size_t voxel = (iz * gridSize.y + iy) * gridSize.x + ix;
                         energy[voxel] = static_cast<float>(std::min(-pocket * energyPerAngstrom, ceiling));
                       }
                     }
                   }
                 }
               });
}


std::vector<BlockingSphere> analysisBlockingSpheres(const PairInteractions &interactions, const Crystal &framework,
                                                    const std::string &probePseudoAtom)
{
  if (!interactions.findType(probePseudoAtom).has_value()) return {};

  VoronoiBlockingSpheres blocks;
  blocks.compute(interactions, framework, probePseudoAtom);
  return blocks.spheres;
}


EnergyBackend withBlockedPockets(EnergyBackend backend, std::span<const BlockingSphere> spheres,
                                 double energyPerAngstrom, double ceiling)
{
  if (spheres.empty()) return backend;

  std::vector<BlockingSphere> held(spheres.begin(), spheres.end());

  std::function<ProbeEnergyGrid(const PairInteractions &, const Crystal &, std::string, uint3)> originalProbe =
      backend.probeEnergyGrid;
  backend.probeEnergyGrid = [originalProbe, held, energyPerAngstrom, ceiling](
                                const PairInteractions &interactions, const Crystal &framework, std::string probe,
                                uint3 gridSize)
  {
    ProbeEnergyGrid grid = originalProbe(interactions, framework, probe, gridSize);
    applyBlockingSpheresToEnergy(grid.energy, grid.gridSize, grid.unitCell, held, energyPerAngstrom, ceiling);
    return grid;
  };

  std::function<MolecularEnergyGrid(const PairInteractions &, const Crystal &, const LinearProbe &, uint3, std::size_t,
                                    double, const ElectrostaticPotentialGrid *)>
      originalMolecular = backend.molecularEnergyGrid;
  backend.molecularEnergyGrid =
      [originalMolecular, held, energyPerAngstrom, ceiling](
          const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe, uint3 gridSize,
          std::size_t numberOfOrientations, double temperature, const ElectrostaticPotentialGrid *potential)
  {
    MolecularEnergyGrid grid =
        originalMolecular(interactions, framework, probe, gridSize, numberOfOrientations, temperature, potential);
    applyBlockingSpheresToEnergy(grid.freeEnergy, grid.gridSize, grid.unitCell, held, energyPerAngstrom, ceiling);
    applyBlockingSpheresToEnergy(grid.minimumEnergy, grid.gridSize, grid.unitCell, held, energyPerAngstrom, ceiling);
    return grid;
  };

  return backend;
}
