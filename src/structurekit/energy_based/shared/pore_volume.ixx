module;

export module energy_shared_pore_volume;

import std;

import uint3;
import crystal;
import pair_interactions;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// How much room a framework offers a molecule, in Å³ and in cm³/g, answered the two ways an energy landscape
// can answer it.
//
// The geometric route counts the points where the probe's centre fits and calls the count a volume. That is
// a definite thing to do because a hard sphere either fits or it does not. On a landscape there is no such
// line, and the two ways of drawing one give genuinely different numbers that are both worth having:
//
// The *geometric reading* draws the line at a level and counts. Below the level the molecule is at least as
// well off as in the gas, so the region is the framework's offer of somewhere to be, and the count of it
// divides into the part that runs through the crystal and the part that closes on itself exactly as the
// clearance route's does. This is the number to set beside a geometric table, and it moves when the level
// moves.
//
// The *thermodynamic reading* does not draw a line at all. It integrates exp(-A/kT) over the cell, which is
// the volume the molecule would occupy if it were spread over the framework in proportion to how long it
// spends in each place. Every point contributes what it is worth and nothing is in or out, so there is no
// level to choose and the answer depends on the temperature instead. This is the number a Henry coefficient
// is built from, and the one an isotherm at low loading is a statement about.
//
// They disagree, and the disagreement is the interesting part. The thermodynamic volume exceeds the
// geometric one wherever the molecule is bound, since a point it sits at constantly counts for more than the
// space it takes up, and falls short of it wherever the room is there but the molecule has no reason to be:
// their ratio is a measure of how much the framework's room is worth to this molecule rather than how much
// of it there is. For helium at room temperature, which is barely held anywhere, the two come together, and
// that is why helium is the probe a void volume is quoted for.
export struct EnergyPoreVolume
{
  double isoValue{0.0};
  double temperature{0.0};

  // The geometric reading, from counting the points below the level.
  double voidFraction{0.0};
  double accessibleFraction{0.0};
  double inaccessibleFraction{0.0};

  double voidVolume{0.0};          // Å³
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  std::size_t numberOfChannels{0};
  std::size_t numberOfPockets{0};
  int dimensionality{0};

  // The same region divided again, by what it costs to leave rather than by whether a path exists at the
  // level. This is the division worth having on a landscape, and on a framework whose windows are tight it
  // is not the one above: nothing in DDR percolates at the contour for methane, so the split by
  // connectivity calls the whole of it inaccessible, while the barrier out of the cages is twelve kT and the
  // methane is in and out of them constantly. The reachable volume below is what a measurement would fill.
  double thresholdInKT{30.0};
  double reachableFraction{0.0};
  double sealedFraction{0.0};
  double reachableVolume{0.0};  // Å³
  double sealedVolume{0.0};     // Å³
  double gravimetricReachableVolume{0.0};

  // The largest barrier anything reachable had to pay, which says how far the answer is from the threshold.
  double largestReachableBarrierInKT{0.0};

  // The same three on the landscape the molecule would see if it were always turned the best way, which
  // contains the one above and so is never smaller. The gap between them is the room a molecule with a shape
  // cannot use because it cannot turn there.
  double minimumEnergyVoidFraction{0.0};
  double minimumEnergyVoidVolume{0.0};

  // The thermodynamic reading: the integral of exp(-A/kT) over the cell, as a fraction of it and as a
  // volume. Above the cell's own volume for anything strongly bound, which is not a contradiction, the
  // quantity not being a fraction of anything at that point.
  double boltzmannFraction{0.0};
  double boltzmannVolume{0.0};  // Å³
  bool readsAsFraction{false};

  // The gravimetric forms, in cm³/g, which is how a pore volume is usually quoted.
  double gravimetricVoidVolume{0.0};
  double gravimetricAccessibleVolume{0.0};
  double gravimetricBoltzmannVolume{0.0};

  double frameworkDensity{0.0};  // g/cm³

  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  double seconds{0.0};

  EnergyPoreVolume();
  ~EnergyPoreVolume();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double isoValue, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, double thresholdInKT = 30.0, bool useElectrostatics = true,
           double relativePrecision = 1.0e-6);
};
