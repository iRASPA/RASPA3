module;

export module energy_shared_well_field;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;

// Types the well-surface field is made of, shared by the processor and GPU builders. The arithmetic that
// fills the field lives next door, behind EnergyBackend; this module is only the numbers and the description
// of the neighbourhood they were built from.
//
// The constants belong here rather than with the surface report so that a backend can use them without
// importing the surface module, which itself imports the backend.

export inline constexpr double wellContactPrefactor = 1.122462048309373;  // 2^(1/6)
export inline constexpr double wellDummyCoreRadius = 1.0;  // Å, massless charge sites only (1/r floor, not a wall)
export inline constexpr double wellSoftminTau = 0.4;                      // Å
export inline constexpr double wellEnergyScalePerKelvin = 0.001;
export inline constexpr double wellFilamentReliabilityThreshold = 0.45;
export inline constexpr double wellFilamentMinimumArea = 2.0;  // Å²
export inline constexpr double wellRefinementReach = 0.75;     // Å

// What one probe site does to one framework atom: the force field's own mixed parameters, the contact
// radius 2^(1/6) sigma of the pair, and the charge product of the near electrostatic half.
export struct SitePair
{
  double epsilon4{0.0};
  double sigma2{0.0};
  double contact{0.0};
  double shift{0.0};
  double chargeProduct{0.0};  // coulombFactor * q_site * q_atom
};

// Where a site hangs off the centre of mass and what it carries. A site may have dispersion, charge, or
// both; one that carries neither was dropped before this was built.
export struct ProbeSite
{
  double offset{0.0};  // Å along the molecular axis, from the centre of mass
  double charge{0.0};
  bool dispersion{true};  // a bare charge site has no pair energy; it still inherits a Pauli core (sigma)
};

// One periodic image of a framework atom near a point: the Cartesian separation from the atom image to the
// point, and which atom it is, so the mixed parameters can be looked up per site.
export struct NearbyImage
{
  double3 dr;
  std::size_t atom;
};

// erfc(alpha r) / r against r², read off with straight lines. Built once per neighbourhood; the GPU copies
// the table into a float buffer and uses the same bins and scale.
export struct ScreenedCoulomb
{
  static constexpr double smallestSquared = 0.25;  // Å², r = 0.5 Å
  static constexpr std::size_t bins = 8192;

  double alpha{0.0};
  double scale{0.0};  // bins per Å²
  std::vector<double> table;

  static double exactly(double alpha, double rr);
  void build(double alphaValue, double largestSquared);
  double at(double rr) const;
};

// The arguments that rebuild the neighbourhood the field (or the refinement) walks. Both backends take this
// rather than the fat CPU Neighbourhood, so the GPU host can pack its own buffers from the same inputs.
export struct NeighbourhoodParameters
{
  LinearProbe probe;
  std::size_t numberOfOrientations{1};
  double thermalEnergy{0.0};
  double extraReach{0.0};
  double blockedEnergyPerAngstrom{0.0};
  double ceiling{0.0};
  const ElectrostaticPotentialGrid *potential{nullptr};
  double coulombFactor{0.0};
};

// One grid of the well field: energy, contact distance and medial reliability, three separate arrays in the
// same voxel order the energy grids use, x varying fastest.
export struct WellField
{
  uint3 gridSize{0, 0, 0};
  UnitCell unitCell;

  std::vector<float> energy;
  std::vector<float> distance;
  std::vector<float> reliability;
  // Per-orientation energy at each voxel, layout [voxel * numberOfOrientations + o], filled when the
  // product is at most 1e9 floats (~4 GB). Empty for a spherical probe. ρ(r, ω) NLDFT reads this; the
  // Helmholtz average in `energy` is still the well-surface field.
  std::vector<float> orientationEnergy;

  // What the field was built with, carried along so that the surface is refined against the same
  // arithmetic that made it. `numberOfOrientations` is the count actually used (one for a probe with no
  // length); `thermalEnergy` is the kT of the orientational average in the field's own units, zero when
  // the least over orientations was taken instead.
  LinearProbe probe;
  std::size_t numberOfOrientations{1};
  double thermalEnergy{0.0};

  // Whether the probe's partial charges were acted on. A charged probe whose electrostatics are left out is
  // short of the term that decides most of what carbon dioxide does, so which of these it is matters more
  // than any other single thing about the field.
  bool chargesIncluded{false};
  bool chargesIgnored{false};
  double ewaldAlpha{0.0};
  std::size_t numberOfWaveVectors{0};

  std::string probeName;
  double cutOff{0.0};
  double ceiling{0.0};
  double seconds{0.0};

  WellField();
  ~WellField();

  std::size_t numberOfVoxels() const { return this->energy.size(); }

  std::size_t voxelIndex(std::size_t i, std::size_t j, std::size_t k) const
  {
    return (k * this->gridSize.y + j) * this->gridSize.x + i;
  }

  double deepestEnergy() const;
};
