module;

export module energy_shared_molecular_energy_grid;

import std;

import int3;
import uint3;
import double3;
import unit_cell;

import energy_shared_linear_probe;

// The energy landscape of a rigid linear molecule in a framework, sampled over the unit cell.
//
// A molecule with a shape does not have one energy at a point, it has an energy for every way it can be
// turned there. Two reductions of that are kept, and they answer different questions.
//
// `minimumEnergy` is the least over orientations, so the landscape a molecule would see if it were always
// turned the best way. It is the zero-temperature limit, and a barrier read from it is a lower bound on the
// true one, since it grants the molecule the best orientation everywhere at no cost.
//
// `freeEnergy` is -kT ln <exp(-U/kT)> over orientations, the landscape left when the turning is averaged out
// rather than optimised. It is the one to use. A molecule threading a narrow window can point almost nowhere,
// and the freedom it gives up to get through is a real cost that no minimum will ever show. The difference
// between the two fields is exactly that cost, so having both is worth more than having either.
//
// They coincide, exactly, when the molecule is a single site at its own centre. That is the case the other
// route computes, and it is the test this one is held to.
//
// The landscape is a value, not a computation. Which arithmetic filled it in is the business of the builder
// that returned it, and nothing downstream of here can tell or needs to.
export struct MolecularEnergyGrid
{
  uint3 gridSize{0, 0, 0};
  UnitCell unitCell;

  // Both indexed by `voxelIndex`, x varying fastest, in the force field's internal units of energy.
  std::vector<float> freeEnergy;
  std::vector<float> minimumEnergy;

  LinearProbe probe;
  std::size_t numberOfOrientations{0};
  bool overHemisphere{true};
  double temperature{0.0};
  double cutOff{0.0};
  int3 numberOfImageShells{0, 0, 0};
  double ceiling{0.0};

  // Whether the molecule's partial charges were acted on. A charged molecule whose electrostatics are left
  // out is short of the term that decides most of what carbon dioxide does, so which of these it is matters
  // more than any other single thing about the landscape.
  bool chargesIncluded{false};
  bool chargesIgnored{false};
  double ewaldAlpha{0.0};
  std::size_t numberOfWaveVectors{0};

  // Which builder filled the field in, as it appears in the names of the files written from it. Nothing
  // computed from the field depends on this; it is here so that running both routes leaves two sets of
  // results side by side rather than one on top of the other.
  std::string backend{"unknown"};

  double seconds{0.0};

  MolecularEnergyGrid();
  ~MolecularEnergyGrid();

  std::size_t numberOfVoxels() const { return this->freeEnergy.size(); }

  std::size_t voxelIndex(std::size_t i, std::size_t j, std::size_t k) const
  {
    return (k * this->gridSize.y + j) * this->gridSize.x + i;
  }

  double3 fractionalPosition(std::size_t i, std::size_t j, std::size_t k) const;
  double3 spacing() const;

  // The deepest point of each landscape.
  double deepestFreeEnergy() const;
  double deepestMinimumEnergy() const;
};
