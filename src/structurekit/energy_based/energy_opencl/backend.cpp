module;

module energy_opencl_backend;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;

import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;
import energy_shared_isosurface;
import energy_shared_energy_backend;
import surface_curvature;

import energy_opencl_probe_energy_grid;
import energy_opencl_molecular_energy_grid;
import energy_opencl_electrostatic_potential_grid;
import energy_opencl_surface_area;
import energy_opencl_lewiner_isosurface;
import energy_opencl_well_field;
import energy_shared_well_field;
import blocking_spheres;
import grid_pore_size_opencl;

EnergyBackend openCLEnergyBackend()
{
  EnergyBackend backend;
  backend.name = "gpu";

  backend.probeEnergyGrid = [](const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                               uint3 gridSize)
  { return ProbeEnergyGridOpenCL::compute(interactions, framework, probePseudoAtom, gridSize); };

  backend.electrostaticPotentialGrid = [](const PairInteractions &interactions, const Crystal &framework, uint3 gridSize,
                                          double relativePrecision)
  { return ElectrostaticPotentialGridOpenCL::compute(interactions, framework, gridSize, relativePrecision); };

  backend.molecularEnergyGrid = [](const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                                   uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                                   const ElectrostaticPotentialGrid *potential)
  {
    return MolecularEnergyGridOpenCL::compute(interactions, framework, probe, gridSize, numberOfOrientations,
                                              temperature, potential);
  };

  backend.molecularStrongestAtom = [](const PairInteractions &interactions, const Crystal &framework,
                                      const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                      double temperature, const ElectrostaticPotentialGrid *potential)
  {
    return MolecularEnergyGridOpenCL::strongestAtoms(interactions, framework, probe, gridSize, numberOfOrientations,
                                                     temperature, potential);
  };

  // Lewiner's tables, compiled once and kept. The field sweep still dominates the cost.
  backend.isosurfaceArea = [](const Crystal &framework, std::span<const float> field, uint3 gridSize,
                              double isoValue)
  {
    std::vector<double3> gradients;
    std::vector<double3> corners = trianglesOfLewinerIsosurface(field, gridSize, isoValue, &gradients);
    return accumulateTriangleAreas(framework.unitCell.cell, gridSize, corners, gradients,
                                   FieldSense::GrowsIntoSolid);
  };

  backend.isosurfaceTriangles = [](std::span<const float> field, uint3 gridSize, double isoValue)
  { return trianglesOfLewinerIsosurface(field, gridSize, isoValue); };

  backend.poreRadiusField = [](uint3 gridSize, const UnitCell &unitCell, std::span<const float> distance,
                               double slack)
  { return poreRadiusFieldOpenCL(gridSize, unitCell, distance, slack); };

  backend.wellField = [](const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                         uint3 gridSize, std::size_t numberOfOrientations, double thermalEnergy,
                         std::span<const BlockingSphere> blockingSpheres, double blockedEnergyPerAngstrom,
                         double ceiling, const ElectrostaticPotentialGrid *potential, double coulombFactor)
  {
    return WellFieldOpenCL::compute(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                    blockingSpheres, blockedEnergyPerAngstrom, ceiling, potential, coulombFactor);
  };

  backend.refineWellVertices = [](std::vector<double3> &corners, std::vector<double> &energies,
                                  const PairInteractions &interactions, const Crystal &framework,
                                  const NeighbourhoodParameters &parameters,
                                  std::span<const BlockingSphere> blockingSpheres, double iso)
  {
    WellFieldOpenCL::refineVertices(corners, energies, interactions, framework, parameters, blockingSpheres, iso);
  };

  return backend;
}
