module;

export module energy_shared_energy_backend;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;

import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;
import energy_shared_isosurface;

// The whole of what the energy-based properties need from the machine underneath: four ways of filling in a
// field, and the one sweep over a field that is heavy enough to be worth choosing a machine for.
//
// Everything above this line is arithmetic on a field that already exists, and there is one copy of it. The
// backends are the two ways of producing the field, and they meet here. Choosing between them is a matter of
// handing a different one of these to the same driver, which is why the routes cannot drift apart in
// anything but the fields themselves.
export struct EnergyBackend
{
  // How this backend is named in the files written from its fields, so that running both leaves two sets of
  // results side by side rather than one on top of the other.
  std::string name;

  std::function<ProbeEnergyGrid(const PairInteractions &, const Crystal &, std::string, uint3)> probeEnergyGrid;

  std::function<ElectrostaticPotentialGrid(const PairInteractions &, const Crystal &, uint3, double)>
      electrostaticPotentialGrid;

  std::function<MolecularEnergyGrid(const PairInteractions &, const Crystal &, const LinearProbe &, uint3, std::size_t,
                                    double, const ElectrostaticPotentialGrid *)>
      molecularEnergyGrid;

  // Which framework atom holds each point of the cell, for a molecule rather than a single site: the one
  // whose own contribution to the energy there is the most negative, averaged over how the molecule is
  // turned. The average is Boltzmann-weighted, not plain, and that is the whole of the definition worth
  // arguing over. A molecule in a window can point almost nowhere, and an average over the orientations it
  // cannot take up is an average over configurations that never happen; weighting by exp(-U/kT) asks which
  // atom pulls hardest on the molecule as it actually sits.
  //
  // Only the part of the energy that belongs to one atom decides the label: the dispersion, and the near
  // half of the electrostatic sum, which is a sum over pairs. The far half is a sum over waves and belongs
  // to the framework as a whole, so it takes no part in the division, though it is present in the weights,
  // which are about how the molecule sits and not about who is responsible for it.
  //
  // With one site at the centre and one orientation this reduces, exactly, to the label the single-site
  // field carries, which is the test it is held to.
  std::function<std::vector<std::int32_t>(const PairInteractions &, const Crystal &, const LinearProbe &, uint3,
                                          std::size_t, double, const ElectrostaticPotentialGrid *)>
      molecularStrongestAtom;

  std::function<IsosurfaceArea(const Crystal &, std::span<const float>, uint3, double)> isosurfaceArea;

  // The same surface as triangles rather than as a number, three corners to a triangle, in fractional
  // coordinates. It is what a surface has to be handed back as before it can be divided among the atoms.
  std::function<std::vector<double3>(std::span<const float>, uint3, double)> isosurfaceTriangles;

  // A sphere about every point of the void, each point keeping the largest that reached it, given the room
  // at every point and how far a sphere may reach past itself. It is the heavy half of a pore-size
  // distribution and the only sweep here that is not dwarfed by building the field it runs on.
  std::function<std::vector<float>(uint3, const UnitCell &, std::span<const float>, double)> poreRadiusField;
};

// A molecule's landscape and the framework potential it was built against, which is empty when the molecule
// carries no charge or the electrostatics were turned off.
export struct MolecularField
{
  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  // Whether this came out of the store below rather than being built. The time in `grid` is the time the
  // landscape took when it was made, which is a true statement about it either way.
  bool reused{false};
};

// Builds both, in the order they depend on each other. The framework's potential does not depend on how the
// molecule is turned, so it is built once and every orientation reads the same one; that is what keeps the
// electrostatics from multiplying the cost of the landscape by the number of orientations.
//
// Every molecular property starts here, so none of them can differ in how the field beneath it was made.
//
// The last landscape built is kept, and asking for the same one again returns it rather than building it a
// second time. That matters because the landscape is the whole cost of these properties: on MFI at 128³ with
// a hundred and twenty-eight orientations it is eighteen seconds against hundredths for anything read off
// it, so a run asking for four properties spent three quarters of its time rebuilding a field it already
// had. Keeping it also removes a way for one report to disagree with another, since two properties asked for
// in one run are now guaranteed to be answers about the same field and not merely about fields built from
// the same arguments.
//
// One is kept rather than several. Properties are asked for in a row against the same molecule and grid, so
// one is all that is ever wanted, and a landscape is large enough that keeping several of them is not free:
// a 128³ field is 16 MB and a 512³ field is a gigabyte. What is held is released before the next is built.
export MolecularField buildMolecularField(const EnergyBackend &backend, const PairInteractions &interactions,
                                          const Crystal &framework, const LinearProbe &probe, uint3 gridSize,
                                          std::size_t numberOfOrientations, double temperature,
                                          bool useElectrostatics, double relativePrecision);

// Lets go of the landscape being held, for a caller that is done with the energy properties and wants the
// memory back. Asking for the same one afterwards builds it again.
export void forgetMolecularField();
