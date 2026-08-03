module;

export module energy_shared_molecular_void_fraction;

import std;

import uint3;
import crystal;
import pair_interactions;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// What a rigid linear molecule would give for the energy-based void fraction, which is the average of
// exp(-U/kT) over every place it could be put.
//
// For a molecule that average runs over orientations as well as positions, and there is nothing to compute
// for it: averaging exp(-U/kT) over orientations is the definition of the free energy already on the grid,
// so the whole property is the spatial average of exp(-A/kT). The same field the barrier climbs, reduced a
// different way. Nothing is recomputed and the two numbers cannot disagree with one another.
//
// The name deserves care. The average is a fraction only for a probe so weakly held that exp(-U/kT) never
// much exceeds one, which is why helium is the probe a void fraction is measured with. Anything that binds
// gives an average above one, sometimes far above, and then the quantity is no longer a fraction of anything
// but the excess Boltzmann factor, which is what the Henry coefficient is built from. Both readings are
// reported and which one is meant is said plainly.
export struct MolecularVoidFraction
{
  // <exp(-A/kT)> over the cell. Below one for a weakly held probe, where it is the void fraction; above one
  // for anything that binds, where it is not.
  double boltzmannAverage{0.0};

  // -kT ln of the above, the work of putting the molecule into the framework from an ideal gas.
  double excessChemicalPotential{0.0};

  // The Henry coefficient the average amounts to, in mol/kg/Pa, which is the form it can be compared against
  // measurement in.
  double henryCoefficient{0.0};

  // Whether the average came out at or below one, and so whether calling it a void fraction means anything.
  bool readsAsFraction{false};

  double temperature{0.0};
  std::size_t numberOfVoxels{0};
  double seconds{0.0};

  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  MolecularVoidFraction();
  ~MolecularVoidFraction();

  // Reduces a landscape that has already been built, for when the barrier has been asked for as well and the
  // two should come from one field.
  static MolecularVoidFraction fromGrid(const MolecularEnergyGrid &grid, const Crystal &framework,
                                        double temperature);

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe, uint3 gridSize,
           std::size_t numberOfOrientations, double temperature, bool useElectrostatics = true,
           double relativePrecision = 1.0e-6);
};
