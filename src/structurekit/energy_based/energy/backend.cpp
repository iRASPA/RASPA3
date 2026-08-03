module;

module energy_backend;

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

import energy_probe_energy_grid;
import energy_molecular_energy_grid;
import energy_electrostatic_potential_grid;
import energy_isosurface;
import grid_pore_size;

EnergyBackend cpuEnergyBackend()
{
  EnergyBackend backend;
  backend.name = "cpu";

  backend.probeEnergyGrid = [](const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                               uint3 gridSize)
  { return ProbeEnergyGridCPU::compute(interactions, framework, probePseudoAtom, gridSize); };

  backend.electrostaticPotentialGrid = [](const PairInteractions &interactions, const Crystal &framework, uint3 gridSize,
                                          double relativePrecision)
  { return ElectrostaticPotentialGridCPU::compute(interactions, framework, gridSize, relativePrecision); };

  backend.molecularEnergyGrid = [](const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                                   uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                                   const ElectrostaticPotentialGrid *potential)
  {
    return MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, numberOfOrientations, temperature,
                                           potential);
  };

  backend.molecularStrongestAtom = [](const PairInteractions &interactions, const Crystal &framework,
                                      const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                      double temperature, const ElectrostaticPotentialGrid *potential)
  {
    return MolecularEnergyGridCPU::strongestAtoms(interactions, framework, probe, gridSize, numberOfOrientations,
                                                  temperature, potential);
  };

  backend.isosurfaceArea = [](const Crystal &framework, std::span<const float> field, uint3 gridSize,
                              double isoValue)
  { return EnergyIsosurface::areaOfIsosurface(framework, field, gridSize, isoValue); };

  backend.isosurfaceTriangles = [](std::span<const float> field, uint3 gridSize, double isoValue)
  { return EnergyIsosurface::trianglesOfIsosurface(field, gridSize, isoValue); };

  backend.poreRadiusField = [](uint3 gridSize, const UnitCell &unitCell, std::span<const float> distance,
                               double slack) { return poreRadiusField(gridSize, unitCell, distance, slack); };

  return backend;
}
