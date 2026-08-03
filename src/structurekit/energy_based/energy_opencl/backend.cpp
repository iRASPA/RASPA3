module;

module energy_opencl_backend;

import std;

import uint3;
import unit_cell;
import crystal;
import pair_interactions;

import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;
import energy_shared_isosurface;
import energy_shared_energy_backend;

import energy_opencl_probe_energy_grid;
import energy_opencl_molecular_energy_grid;
import energy_opencl_electrostatic_potential_grid;
import energy_opencl_surface_area;
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

  // The extractor holds compiled kernels, so one is built per call rather than kept alive between them. The
  // iso-surfaces are extracted a handful of times per run and the landscape they are extracted from costs far
  // more than the compilation does.
  backend.isosurfaceArea = [](const Crystal &framework, std::span<const float> field, uint3 gridSize,
                              double isoValue)
  {
    EnergyOpenCLSurfaceArea extractor;
    return extractor.areaOfIsosurface(framework, field, gridSize, isoValue);
  };

  backend.isosurfaceTriangles = [](std::span<const float> field, uint3 gridSize, double isoValue)
  {
    EnergyOpenCLSurfaceArea extractor;
    return extractor.trianglesOfIsosurface(field, gridSize, isoValue);
  };

  backend.poreRadiusField = [](uint3 gridSize, const UnitCell &unitCell, std::span<const float> distance,
                               double slack)
  { return poreRadiusFieldOpenCL(gridSize, unitCell, distance, slack); };

  return backend;
}
