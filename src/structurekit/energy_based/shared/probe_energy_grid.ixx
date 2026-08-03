module;

export module energy_shared_probe_energy_grid;

import std;

import int3;
import uint3;
import double3;
import unit_cell;
import crystal;

// The interaction energy of a single probe atom with the framework, sampled on a regular grid over the unit
// cell. It is the energetic counterpart of the clearance field, and the two are laid out alike so that what
// is built on one can be run on the other.
//
// `energy` holds the sum over framework atoms and their periodic images within the cutoff of the
// Lennard-Jones pair energy, with the size and strength parameters mixed for the probe by the force field's
// own rule. The value is large and positive where the probe overlaps an atom, negative in the pockets and
// channels where it is held, and near zero in the middle of a wide pore.
//
// Where the clearance field asks whether a hard sphere fits, this one asks what it would cost to be there,
// and the two answer differently in the place that matters: a window a probe does not geometrically pass is
// often one a real molecule crosses many times a second, because the barrier is a few kT rather than
// infinite. That is the whole reason for computing it.
//
// The two fields differ in cost as well as in meaning. The clearance is a minimum, so an atom that cannot
// win is dismissed unseen and a framework of any size is nearly all dismissed. An energy is a sum, and there
// is no dismissing a term of a sum: every atom within the cutoff has to be evaluated, and the cutoff reaches
// further than one cell. This field is accordingly the more expensive of the two by a wide margin.
//
// The field is a value, not a computation. Which arithmetic filled it in is the business of the builder that
// returned it, and nothing downstream of here can tell or needs to.
export struct ProbeEnergyGrid
{
  uint3 gridSize{0, 0, 0};
  UnitCell unitCell;

  // Indexed by `voxelIndex`, x varying fastest, in the force field's internal units of energy.
  std::vector<float> energy;

  // Which framework atom holds the point, in the order of `framework.atoms`: the one whose own
  // contribution to the sum, over its periodic images within the cutoff, is the most negative. It is the
  // energetic counterpart of the clearance field's nearest-surface label, and it divides the cell into cells
  // in the same way, so that a volume or an area can be handed out among the atoms.
  //
  // Nearest and strongest are not the same question, and where they differ they differ for a reason. Deep in
  // a pocket the atom that pulls hardest is almost always the nearest one and the two tessellations agree.
  // In a window between two cavities they need not: a nearby atom sitting near the bottom of its well pulls
  // weakly, while one a little further out on the steep flank of the potential can dominate. It is those
  // places, which are the ones that decide what a molecule can do, that this label describes and the
  // geometric one does not.
  //
  // The rule needs a tie-break where no atom is attractive at all, which is the inside of a wall. There the
  // nearest atom is taken instead, which is what the geometric route would say about that region anyway.
  // Energy cannot decide it: every term is held at the ceiling inside a wall, so the atoms come out exactly
  // equal and the strongest of them is whichever was looked at first, which is not the same on two machines.
  std::vector<std::int32_t> strongestAtom;

  // The probe the field was built for, and the cutoff and image range it was built with.
  std::string probeName;
  double probeEpsilon{0.0};
  double probeSigma{0.0};
  double cutOff{0.0};
  int3 numberOfImageShells{0, 0, 0};

  // What the energy is held down to where the probe overlaps an atom, so that the field stays finite and can
  // be sorted. It is far above any barrier a molecule crosses, so nothing that is asked of the field is
  // decided by it.
  double ceiling{0.0};

  // Which builder filled the field in, as it appears in the names of the files written from it. Nothing
  // computed from the field depends on this; it is here so that running both routes leaves two sets of
  // results side by side rather than one on top of the other.
  std::string backend{"unknown"};

  double seconds{0.0};

  ProbeEnergyGrid();
  ~ProbeEnergyGrid();

  std::size_t numberOfVoxels() const { return this->energy.size(); }

  std::size_t voxelIndex(std::size_t i, std::size_t j, std::size_t k) const
  {
    return (k * this->gridSize.y + j) * this->gridSize.x + i;
  }

  double3 fractionalPosition(std::size_t i, std::size_t j, std::size_t k) const;
  double3 cartesianPosition(std::size_t i, std::size_t j, std::size_t k) const;

  // The deepest point of the landscape, which is the strongest adsorption site the grid resolves.
  double minimumEnergy() const;

  // The grid spacing along each cell axis, as a distance rather than a fraction.
  double3 spacing() const;

  // Writes out the tessellation: which atom holds every grid point and what it is worth there. The cells are
  // the strongest-attraction cells of the atoms, which is what the volume and the area of this route are
  // divided along, and they are the energetic counterpart of the nearest-surface cells the clearance route
  // writes out in the same form.
  void writeTessellation(const Crystal &framework) const;
};

// Well above any barrier a molecule crosses and well below what a float cannot hold, so that the overlap
// region stays finite and orderable without any real answer depending on where it was put. Every builder
// uses this same figure, or the fields they return would not be comparable.
export inline constexpr double probeEnergyCeilingInKelvin = 1.0e10;
