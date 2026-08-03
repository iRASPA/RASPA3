module;

export module energy_shared_molecular_energy_barrier;

import std;

import uint3;
import crystal;
import pair_interactions;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_energy_barrier;
import energy_shared_electrostatic_potential_grid;

// The percolation barrier for a molecule that has a shape, computed twice over.
//
// The same sweep is run on both landscapes the molecular field produces. From the minimum over orientations
// comes the barrier a molecule would face if it were always turned the best way, which is the zero
// temperature limit and a lower bound. From the orientational free energy comes the barrier it actually
// faces at temperature, the one that counts the freedom it must give up to thread a narrow window.
//
// Reporting both is the point. Either alone is a number; the two together separate what the barrier costs in
// energy from what it costs in orientation, and the second of those is invisible to every geometric route
// and to the single-site energy route alike.
export struct MolecularEnergyBarrier
{
  EnergyBarrier fromFreeEnergy;
  EnergyBarrier fromMinimumEnergy;

  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  // How much higher the free-energy barrier stands than the best-oriented one. This is the orientational
  // part of the barrier, and it is positive by construction.
  double orientationalPenalty{0.0};

  MolecularEnergyBarrier();
  ~MolecularEnergyBarrier();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe, uint3 gridSize,
           std::size_t numberOfOrientations, double temperature, bool useElectrostatics = true,
           double relativePrecision = 1.0e-6);
};
