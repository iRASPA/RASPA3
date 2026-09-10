module;

export module energy_shared_well_surface;

import std;

import uint3;
import double3;
import crystal;
import pair_interactions;
import unit_cell;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_energy_backend;
import energy_shared_well_field;

export import energy_shared_well_field;

// The well surface: where an adsorbed molecule sits, as opposed to the iso-surface, which is where it turns
// back.
//
// Topology comes from a distance field and geometry from the energy.
//
//   - `computeWellField` supplies three numbers per grid point: the energy U, the additively weighted
//     (Apollonius) distance d = min over atoms of (|x - a| - 2^(1/6) sigma), whose zero level set is the
//     probe-contact offset surface of the framework, and the medial reliability rel (1 against one wall, 0
//     on the medial axis of a channel where opposing walls cancel). Blocking spheres join d through the min
//     and ramp U, so a pocket closes the surface the way a sphere of framework would.
//   - The field handed to marching cubes is max(-d, s (U - iso)): its zero set bounds the region
//     { d > 0 and U < iso }, the pore trimmed to wells at least iso deep, with smooth isosurface caps at the
//     trim. Being the boundary of a region, it is watertight and single-sheeted --- no interior membranes,
//     domes, or flaps, which were the failure modes of walking down the energy from the zero surface.
//   - Each vertex is then slid along the ray to its nearest atom onto the exact 1D minimum of the analytic
//     energy: the true multi-atom well floor.
//   - Where the channel is narrower than the probe's contact diameter the sheet cannot exist (d < 0 across
//     the whole cross-section), yet those are the deepest wells of all: the transverse minima have merged
//     onto the channel axis, a 1D filament. It is extracted as a thin tube, the zero set of
//     max(rel - r0, d, s (U - iso)). Packing does not use that tube as a cylinder. The ridge is
//     the Gelb-Gubbins inscribed radius of the energy iso (distance to U = iso, the first half of
//     the energy PSD), restricted to the filament: a wrapping tube is a 1-D file (whole channels
//     times whole sites along the period), a slab wrapping two axes is a 2-D midplane. A flattened
//     ribbon that still wraps one axis is a file, not a midplane. Wrapping is read from the medial
//     tight-channel set, not from {U < iso}, so a pinched energy well still packs as a file.
//     A²/(4πV) of the whole mesh is kept only as a check.
//
// A molecule with a shape of its own goes through the same construction with its two scalars generalized.
// Its energy at a centre-of-mass position is the orientational free energy -kT ln <exp(-U/kT)> over a set
// of axis directions (the least of them at kT = 0), and its contact distance is the largest over directions
// of the least (|x_site - a| - 2^(1/6) sigma) over site-atom pairs: zero where the best-fitting way round
// just touches the framework, negative where no way round fits. The dimer that cannot pass a window head-on
// still slips through sideways, and the max over directions is what says so. Sites carrying charge but no
// dispersion (the massless centre site of N2 and O2) have no pair energy and no contact radius, and a field
// with no electrostatics in it owes them nothing. A single site at the centre is the same molecule every
// way round, so one direction stands for the sphere and both reductions collapse, exactly, to the
// single-atom definitions --- through the same code, which is the sharpest test of it.
//
// A charged probe handed the framework's electrostatic potential feels it: each site adds q times the
// tabulated far half of the Ewald sum, and the near half, erfc(alpha r)/r pair by pair, rides along in the
// same walk over neighbours the dispersion is already making. The quadrupole this puts on nitrogen and
// carbon dioxide is most of what they do in a charged framework, and it deepens and turns the wells; it
// takes no part in the contact distance, because charges do not change where the molecule fits, only how it
// sits. The conversion factor of charge² per Ångström into energy is taken as an argument for the same
// reason the temperature is: the measurement is not allowed to know the unit system.
//
// This replaces an earlier construction that walked from every vertex of the zero-energy surface along the
// wall normal until the energy stopped falling. That walk is the right physics on a smooth isolated wall and
// the wrong mesh on a real framework: the crease set of the energy contains one-sided sheets --- window
// membranes, intersection domes, flaps over wall bumps --- and no local quantity separates them from the
// wall sheet. The Apollonius contact surface cannot produce those, because d has exactly one zero crossing
// along any ray into a wall.

export struct SheetPatch
{
  double area{0.0};  // Å²
  std::array<double, 3> energy{0.0, 0.0, 0.0};
};

export struct FilamentVoxel
{
  double volume{0.0};  // Å³ of the tube of centres
  double length{0.0};  // Å along a 1-D file; 0 means this voxel is not packed as a line
  double area{0.0};    // Å² of a 2-D midplane; 0 means this voxel is not packed as a sheet
  double energy{0.0};
};

export struct WellSurface
{
  double area{0.0};  // Å² per unit cell, the refined sheet
  double gravimetricArea{0.0};  // m²/g
  double volumetricArea{0.0};   // m²/cm³

  double weightedArea{0.0};
  double gravimetricWeightedArea{0.0};

  double filamentArea{0.0};
  double filamentVolume{0.0};  // Å³ of voxels that belong to the merged-well filament
  double filamentLength{0.0};  // Å, cylinder A²/(4πV) of the whole overlay; diagnostic only
  double filamentRidgeLength{0.0};  // Å, 1-D packing: graph length of the medial ridge
  double filamentMedialArea{0.0};  // Å², 2-D packing: half the planar overlay (area / 16.2 in BET)

  std::size_t numberOfTriangles{0};
  std::size_t numberOfRejectedTriangles{0};
  std::size_t numberOfFilamentTriangles{0};
  std::size_t numberOfTrimmedVertices{0};

  double meanDepth{0.0};
  double deepestWell{0.0};
  double isoValue{0.0};
  double energyScale{0.0};
  double thermalEnergy{0.0};
  double seconds{0.0};

  std::vector<SheetPatch> patches;
  std::vector<FilamentVoxel> filamentVoxels;

  double enhancement() const { return (this->area > 0.0) ? this->weightedArea / this->area : 0.0; }
};

// `thermalEnergy` is the kT of the orientational free energy, in the force field's internal units; at or
// below zero the least over orientations is taken instead. A probe with no length ignores both it and
// `numberOfOrientations`. A charged probe acts on `potential` when one is supplied, and `coulombFactor` is
// then the conversion of charge² per Ångström into the field's units --- it must be the same one the
// potential was built with, or the two halves of the Ewald sum will not add up to anything.
export WellField computeWellField(const PairInteractions &interactions, const Crystal &framework,
                                  const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                  double thermalEnergy, std::span<const BlockingSphere> blockingSpheres,
                                  double blockedEnergyPerAngstrom, double ceiling,
                                  const ElectrostaticPotentialGrid *potential = nullptr,
                                  double coulombFactor = 0.0);

// The trim level actually used: the isovalue, unless that is below the deepest well on the grid (nothing
// would remain), in which case a quarter of the deepest well.
export float effectiveTrimIsovalue(const WellField &field, double isovalue);

// Maps the well field onto the contact sheet, refines the vertices onto the analytic well floor, extracts
// the merged-well filament, and measures them. The probe, the orientations and the kT of the orientational
// average are read off the field, so the refinement works the same arithmetic that built it. `energyScale`
// is angstrom per internal energy unit (iRASPA's 0.001 Å/K times EnergyToKelvin); the zero set itself does
// not depend on it. `thermalEnergy` is kT in the same units as the field, for the Boltzmann weight of the
// area, or zero to leave that weight out.
export WellSurface wellSurfaceOfField(const Crystal &framework, const PairInteractions &interactions,
                                      const WellField &field, double isoValue, double energyScale,
                                      double thermalEnergy, std::span<const BlockingSphere> blockingSpheres,
                                      const EnergyBackend *backend = nullptr,
                                      const ElectrostaticPotentialGrid *potential = nullptr,
                                      double coulombFactor = 0.0);

export void writeWellSurface(std::ostream &stream, const WellSurface &surface);

export struct MolecularWellSurface
{
  WellSurface surface;
  WellField field;
  double isoValue{0.0};

  MolecularWellSurface();
  ~MolecularWellSurface();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double level, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, std::span<const BlockingSphere> blockingSpheres = {}, double energyScale = 0.0,
           double blockedEnergyPerAngstrom = 0.0, double ceiling = 0.0, bool useElectrostatics = true,
           double relativePrecision = 1e-6);
};
