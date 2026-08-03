module;

export module energy_shared_tessellation;

import std;

import uint3;
import crystal;
import pair_interactions;
import energy_shared_energy_backend;
import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;

// What each framework atom is answerable for, once the cell is divided by which atom pulls hardest on the
// probe.
//
// The geometric route divides the cell by whose surface is nearest, which is the Apollonius tessellation,
// and hands each atom the volume and the area that fall in its cell. That is a statement about where the
// atoms are and nothing else: the probe enters only through its radius, and two probes of the same size get
// the same answer however differently they interact.
//
// This divides the same cell by whose attraction is strongest instead. It is a statement about the atoms and
// the probe together, and it moves: a probe that binds strongly to oxygen and weakly to silicon will find
// the oxygens claiming ground that a hard sphere of the same size would have given to the silicons. Which of
// the two divisions is wanted depends on the question. For packing, the geometric one; for anything about
// where a molecule sits and what holds it there, this one.
//
// The two agree wherever the nearest atom is also the strongest, which is most of a pocket. They part in the
// windows, and in a framework of more than one element.
//
// A division of this kind is decided by an argument of a minimum, so the boundaries between cells are where
// two atoms tie and the winner is settled by rounding. In a symmetric framework, where a great many atoms
// are exactly equivalent, those ties are everywhere along the boundaries, and the processor and the GPU
// break them differently: the fields agree to a thousandth of a Kelvin, and still about one grid point in a
// hundred is handed to a different atom at 48³. It is confined to a surface, so its share of the cell falls
// off as the reciprocal of the grid size, from about 1.9% of the points at 32³ to 0.6% at 96³. It is not a
// disagreement about the physics and the geometric tessellation has exactly the same property; only per-atom
// shares are affected, and totals are not.
export struct EnergyAtomShare
{
  std::size_t numberOfVoxels{0};

  // The share of the cell this atom holds, and the share of the accessible void within it.
  double volume{0.0};
  double volumeFraction{0.0};
  double voidVolume{0.0};

  // The part of the iso-surface that falls in this atom's cell, in Å² per unit cell and per unit mass.
  double area{0.0};
  double gravimetricArea{0.0};

  // What the probe feels over this atom's cell, summed where it is held at all. It is the atom's share of
  // the binding the framework offers, which is the one column here that has no geometric counterpart: a
  // nearest-surface cell divides space, and space does not add up to an energy.
  double bindingEnergy{0.0};

  // The deepest single point in this atom's cell, which is the site the atom offers.
  double deepestEnergy{0.0};
};

export struct EnergyTessellation
{
  double isoValue{0.0};

  std::vector<EnergyAtomShare> atoms;

  // Totals, so that the shares can be checked to add up to what the undivided route reports.
  double totalArea{0.0};
  double totalGravimetricArea{0.0};
  double totalVoidVolume{0.0};
  double totalBindingEnergy{0.0};

  // Surface that fell in no atom's cell, which should be none of it and is reported so that it can be seen
  // to be none of it.
  double undecidedArea{0.0};

  std::size_t numberOfTriangles{0};
  std::size_t numberOfRejectedTriangles{0};

  // Only one of these holds a field, according to which of the two runs below was called. The molecular one
  // is empty for a single site and the single-site one is empty for a molecule.
  ProbeEnergyGrid grid;
  MolecularEnergyGrid molecularGrid;
  bool isMolecular{false};

  double seconds{0.0};

  EnergyTessellation();
  ~EnergyTessellation();

  // The cell divided for one probe atom, which is the cheaper of the two and the one whose label the field
  // already carries.
  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           std::string probePseudoAtom, double isoValue, uint3 gridSize);

  // The same division for a molecule, which has to be turned every way at every point before it can be said
  // which atom pulls hardest on it. See `molecularStrongestAtom` in the backend for what the label means.
  //
  // It costs about twice what the landscape costs, since the energy of every orientation is wanted once to
  // weigh the orientations and once more to divide the weighted energy among the atoms. The landscape itself
  // is shared with the other molecular properties and is not built again here.
  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double isoValue, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, bool useElectrostatics = true, double relativePrecision = 1.0e-6);
};
