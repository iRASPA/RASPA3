module;

export module energy_void_fraction;

import std;

import int3;
import double2;
import double3;
import double3x3;
import pair_interactions;
import crystal;
import blocking_spheres;

// The average of exp(-U/kT) over the cell for one probe atom, by Widom insertion: the probe is dropped at random
// points and the weights are averaged. It is the sampled estimate of what the grid routes sum on a lattice, and
// shares no arithmetic with them.
export struct EnergyVoidFraction
{
  EnergyVoidFraction() {};

  // A Boltzmann average is of a temperature, and this one is taken at the one given. A point inside a blocking
  // sphere is in a pocket the probe cannot reach, and it is given what the field routes give it there, the
  // pocket's own energy, `blockedEnergyPerAngstrom` for every Ångström it is inside the sphere and no more than
  // `ceiling`, in place of the energy of the atoms round it. Nought for either is the usual figure: the ramp
  // of `blockedEnergyPerAngstromInKelvin`, and the ceiling every energy field is held at.
  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnersteps,
           double temperature = 298.0, std::span<const BlockingSphere> blockingSpheres = {},
           double blockedEnergyPerAngstrom = 0.0, double ceiling = 0.0);

  // The average of exp(-U/kT), as last computed; also written to the output file.
  double boltzmannAverage{0.0};
  std::size_t numberOfInsertions{0};
};
