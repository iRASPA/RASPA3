module;

export module structure_input;

import std;

import simulationbox;
import framework;
import forcefield;

import unit_cell;
import crystal;
import pair_interactions;

// The one place where the simulation engine's idea of a structure is turned into the structural analysis's.
//
// structurekit is built on the same three kits as raspakit and not on raspakit itself, so nothing there can
// be handed a Framework or a ForceField. That is the point of it: a pore diameter is a property of a set of
// spheres in a cell, and a library that says so can be read, tested and trusted without the engine. The cost
// is this translation, and it is cheap because so little of a framework and a force field is a question about
// a pore: nine fields of the one, and three numbers per pair of types from the other.
//
// This lives above both libraries rather than inside either, so that the arrow still only ever points from
// the engine-aware code into the analysis, and never back.
export namespace StructureInput
{
// The cell, carrying over whether it was found to be rectangular rather than deciding again, so that the
// minimum image convention is taken in exactly the form the engine would have taken it.
UnitCell makeUnitCell(const SimulationBox& simulationBox);

// One unit cell of the framework: where its atoms are, what they carry, and what a report should be headed
// with. Symmetry has already been applied by the time a Framework holds these, so there is nothing to expand
// here.
Crystal makeCrystal(const Framework& framework);

// The force field flattened into a square table of pairs. The mixing rule, the potential shift and the
// tail correction have all been applied by the force field already; what is copied out is the result.
//
// Only the Lennard-Jones size and strength survive the crossing, because that is the only form the analysis
// evaluates. A structure given a Buckingham or a Morse force field will be measured as though the pairs were
// Lennard-Jones with a size and a strength of zero, which is to say not measured at all, and it is better
// that this is stated in one place than discovered in nine grid loops.
PairInteractions makeInteractions(const ForceField& forceField);
}  // namespace StructureInput
