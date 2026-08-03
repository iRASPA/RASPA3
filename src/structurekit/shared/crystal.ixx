module;

export module crystal;

import std;

import double2;
import double3;

export import unit_cell;

import pair_interactions;

// One atom of the structure, as everything here reads it: where it is, what it carries, and which row of the
// interaction table it uses.
export struct CrystalAtom
{
  double3 position{};   // Cartesian, Å
  double charge{0.0};   // e
  std::size_t type{0};  // index into PairInteractions
};

// The contents of one unit cell, and the object the diagram and energy routines are given.
//
// It is deliberately not a framework in the sense the simulation engine means. A framework there is a
// component of a system: it has molecules, groups, connectivity, intramolecular potentials, fractional
// molecules and a Widom average, all of which exist because something is going to move in it. Nothing moves
// here. A pore diameter is a property of a fixed set of spheres in a fixed cell, and the fixed set of spheres
// in a fixed cell is all of what this holds, which is what lets this library be compiled and tested without
// the engine at all.
//
// The name and the space group are not used in any calculation. They are what the reports are headed with,
// and they travel with the structure so that a routine can write its own file without knowing where the
// structure came from. `SampledStructure` is the same idea taken one step further, for the samplers that need
// only a contact radius per atom and not a type or a charge.
export struct Crystal
{
  std::string name;
  std::size_t spaceGroupHallNumber{1};

  UnitCell unitCell;

  double mass{0.0};  // g/mol, one unit cell

  std::vector<CrystalAtom> atoms;

  // The same atoms in the same order, in fractional coordinates wrapped into [0, 1). Carried rather than
  // derived because most of the grid code works in fractional coordinates throughout and would otherwise
  // convert the whole cell again for every grid it builds.
  std::vector<double3> fractionalPositions;

  std::size_t size() const { return atoms.size(); }

  std::vector<double3> cartesianPositions() const;

  // Each atom's own size and strength, unmixed, in the order the atoms are stored. What wants this is the
  // device code, which cannot index a table of pairs from inside a kernel as cheaply as it can carry two
  // numbers per atom.
  std::vector<double2> lennardJonesParameters(const PairInteractions& interactions) const;
};
