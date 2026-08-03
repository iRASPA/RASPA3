module;

export module energy_shared_pore_analysis;

import std;

import uint3;
import crystal;
import pair_interactions;
import grid_connected_components;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_energy_barrier;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// The pore diameters, twice over: as levels of the energy landscape, and as lengths on a contour of it.
//
// Nearly everything the geometric route computes is a functional of one scalar field, the clearance, and
// the energy route replaces that field with the landscape. The useful division of the properties is then
// between the ones that are *levels* of a field and the ones that are *sizes of a region*, because the two
// carry over quite differently.
//
// Di, Df and Dif are levels. Di is the largest value the clearance takes, Df the largest value that still
// has a path through the crystal at or above it, and Dif the largest value on that path. None of the three
// needs a boundary drawn anywhere: they are answers about the field itself. Run the identical sweep on the
// landscape and all three come back with no parameter at all, in Kelvin:
//
//   Di   the deepest well anywhere, the strongest site the grid resolves
//   Df   the percolation barrier, the least the molecule must be raised to get from a cell to the next
//   Dif  the deepest well on the percolating network, rather than one in a pocket beside it
//
// These are the primary answers, and they are better questions than their geometric namesakes. A window a
// hair too narrow is as shut as a wall to a hard sphere and is crossed constantly by a molecule, and the
// difference is nowhere in the geometry but is the whole of the barrier.
//
// The lengths are the secondary answers, and they are here for comparability rather than because they are
// better. A diameter in Ångström needs a boundary, an energy landscape has none, so a level has to be named
// and the answer is then a statement about that level as much as about the framework. What is measured is
// the largest sphere that fits inside the contour, which is the same construction the clearance route runs,
// on a region the energy drew instead of one the atoms drew. Reading them beside a geometric table is the
// thing they are good for; quoting one on its own invites it to be taken for a measurement it is not.
export struct EnergyPoreDiameters
{
  // In Ångström, on the region below the iso-value.
  double includedSphereDiameter{0.0};
  double freeSphereDiameter{0.0};
  double includedAlongFreePathDiameter{0.0};
  std::array<double, 3> freeSphereDiameterByDimension{0.0, 0.0, 0.0};

  bool percolates{false};
  int dimensionalityAtThreshold{0};
  std::size_t numberOfVoidVoxels{0};
  double seconds{0.0};
};

export struct EnergyPoreAnalysis
{
  double isoValue{0.0};
  double temperature{0.0};

  // The three levels, on the landscape the molecule actually sees and on the one it would see if it were
  // always turned the best way. The difference between them is what the molecule pays in orientation.
  EnergyBarrier levels;
  EnergyBarrier minimumEnergyLevels;
  double orientationalPenalty{0.0};

  // The three lengths, on the contour.
  EnergyPoreDiameters diameters;

  // The pieces of the region at that level, and what each of them is.
  GridComponents channels;

  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  double seconds{0.0};

  EnergyPoreAnalysis();
  ~EnergyPoreAnalysis();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double isoValue, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, bool useElectrostatics = true, double relativePrecision = 1.0e-6);
};
