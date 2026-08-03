module;

export module energy_shared_electrostatic_potential_grid;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;

// The framework's electrostatic potential, in energy per unit charge, so that a site of charge q sitting at a
// point feels q times what is stored here.
//
// Ewald splits the lattice sum in two. The near part falls off fast enough to be summed over neighbours
// directly, and the far part is smooth enough to be summed over a handful of wave vectors instead. Only the
// far part is tabulated here. That is deliberate: it is smooth, being built from waves no shorter than a few
// Ångström, so reading it off between grid points costs almost nothing in accuracy. The near part has a 1/r
// in it and is anything but smooth, so it is left to be summed exactly wherever it is wanted, which is cheap
// because whatever needs it is already walking over the same neighbours for the dispersion.
//
// Where the split is put is arbitrary, and that is worth something: the total cannot depend on it. Moving
// alpha and finding the answer unchanged tests the near part against the far part, and neither can be wrong
// on its own without the agreement failing.
export struct ElectrostaticPotentialGrid
{
  uint3 gridSize{0, 0, 0};
  UnitCell unitCell;

  // The far part of the potential, plus the constant that goes with a charged cell, over the grid.
  std::vector<float> smoothPotential;

  // Where the split between near and far was put, in inverse Ångström. The near part must be summed with
  // this same value, or the two halves will not add up to anything.
  double alpha{0.0};

  double cutOff{0.0};
  double relativePrecision{0.0};
  double largestWaveVector{0.0};
  std::size_t numberOfWaveVectors{0};

  // The cell's total charge. A framework that does not add up to zero is neutralised by a uniform background,
  // which is the usual thing to do but leaves the potential shifted by a constant that no longer cancels.
  double netCharge{0.0};

  // The largest charge on any framework atom. A structure read from a file with no charges in it leaves this
  // at zero, and then the potential is zero everywhere and the electrostatics are doing nothing at all.
  double largestFrameworkCharge{0.0};

  // Which builder filled the field in, as it appears in the names of the files written from it. Nothing
  // computed from the field depends on this; it is here so that running both routes leaves two sets of
  // results side by side rather than one on top of the other.
  std::string backend{"unknown"};

  double seconds{0.0};

  ElectrostaticPotentialGrid();
  ~ElectrostaticPotentialGrid();

  std::size_t numberOfVoxels() const { return this->smoothPotential.size(); }

  // The whole potential at one point, near part and far, summed on the processor in double precision. Slow,
  // and meant for checking rather than for use. It reads the tabulated far part off this grid and sums the
  // near part itself, so it is the same arithmetic whichever builder filled the grid in.
  double exactPotentialAt(const Crystal &framework, double3 fractionalPosition) const;

  // Recomputes the whole potential at a few points for several splits. The answers must agree; how well they
  // do is the measure of whether the two halves are right. It is built on `exactPotentialAt` alone, so it
  // tests the Ewald arithmetic rather than either backend, and both of them offer the same check.
  static std::string splitIndependenceCheck(const PairInteractions &interactions, const Crystal &framework);
};
