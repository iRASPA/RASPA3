module;

export module energy_shared_blocking_spheres;

import std;

import int3;
import uint3;
import double3;
import crystal;
import pair_interactions;
import grid_blocking_cover;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// One region the molecule can sit in, and what it would cost to leave.
export struct EnergyCavity
{
  std::size_t pore{0};

  // The volume it holds and the deepest point in it, the latter being where the molecule will be found.
  double volume{0.0};  // Å³
  std::size_t numberOfVoxels{0};
  double deepestEnergy{0.0};

  // The lowest energy at which a path runs from this region to somewhere that leaves the crystal, and the
  // climb to it from the bottom. Highest double for a region with no way out at any energy, which is the
  // hard pocket of the geometric route.
  double escapeEnergy{0.0};
  double escapeBarrier{0.0};

  // The same barrier in units of kT, which is the number that decides the matter.
  double barrierInKT{0.0};

  bool blocked{false};
  double3 centreFractional{0.0, 0.0, 0.0};
};

// Blocking spheres from the energy landscape, where what is sealed off is decided by a barrier rather than
// by connectivity.
//
// The geometric route calls a pocket a piece of the void that does not run anywhere, and a probe put into
// one can never leave, because a hard sphere either fits through the neck or it does not. That is the right
// answer to the question it asks and the wrong answer to the question a simulation is asking. A molecule
// gets through a window it does not fit through, by paying for it, and how often depends on the price: a
// neck a hair too narrow costs a few kT and is crossed constantly, and one much too narrow costs hundreds
// and is never crossed at all. Between those two the geometric route sees no difference whatever.
//
// So the void is divided here by what it costs to leave it. The region is the part of the landscape the
// molecule can occupy, its pieces are found as before, and for each piece the sweep says at what energy a
// path first runs from it to somewhere that leaves the crystal. The climb to that from the bottom of the
// piece is the escape barrier, and a piece is blocked when the barrier is large against kT.
//
// What makes that threshold defensible rather than arbitrary is how weakly it depends on what one means by
// large. A molecule in a well attempts the barrier at something like a lattice frequency, so it leaves after
// about exp(barrier/kT) attempts, and the barrier that holds it for a time t is kT ln(nu t). Between a
// second and a year that logarithm runs from 28 to 38. The whole range of times anyone could mean by
// "cannot get in" is thus a range of ten in the threshold, and a barrier is almost never within ten kT of
// it: the necks that seal cages are hundreds of kT and the ones that merely slow a molecule down are tens.
// The cases that fall in between are worth knowing about, and they are reported rather than decided.
//
// The comparison with the geometric answer is the point of having this. A region sealed by geometry whose
// barrier is small is a region the geometric route blocks and this one does not, and blocking it would have
// thrown away uptake the framework really has.
export struct EnergyBlockingSpheres
{
  // The energy below which the molecule is taken to be in the void. Zero, the level at which it is no better
  // off than in the gas, is the reading of contact that costs nothing to defend.
  double isoValue{0.0};
  double temperature{0.0};

  // How many kT a region must cost to leave before it is called blocked.
  double thresholdInKT{30.0};

  std::vector<GridBlockingSphere> spheres;
  std::vector<EnergyCavity> cavities;

  // Of the pieces that do not run anywhere: the ones deep enough to block, and the ones a molecule gets out
  // of anyway. The second is what this route has to say that the geometric one does not.
  std::size_t numberOfBlockedCavities{0};
  std::size_t numberOfLeakyCavities{0};
  std::size_t numberOfChannels{0};

  // The shares of the cell held by the void, by the blocked part of it, and by the part that is sealed by
  // geometry but not by the barrier. The last is the volume the geometric route would have blocked in error.
  double voidFraction{0.0};
  double blockedFraction{0.0};
  double leakyFraction{0.0};

  std::size_t numberOfClippedSpheres{0};
  std::size_t numberOfRefusedPoints{0};

  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  double sweepSeconds{0.0};
  double coverSeconds{0.0};

  EnergyBlockingSpheres();
  ~EnergyBlockingSpheres();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double isoValue, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, double thresholdInKT = 30.0, bool useElectrostatics = true,
           double relativePrecision = 1.0e-6);
};
