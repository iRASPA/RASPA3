module;

export module energy_shared_bet_surface_area;

import std;

import uint3;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_well_surface;

// A predicted BET surface area, built from the well surface and the energy-grid Henry coefficient rather
// than from an experimental isotherm.
//
// The isotherm is a Fowler–Guggenheim / Bragg–Williams mean-field lattice gas on the geometric adsorption
// sites the well surface already supplies: each contact-sheet triangle holds area / 16.2 Å² sites at its
// well-floor energy U; a merged well packs as a 1-D file (L / v_L^{1/3}) or a 2-D midplane (area / 16.2),
// again at the filament energy. Molecules on a site attract each other with a mean-field coupling
// W = -2 q_L, fixed so bulk liquid nitrogen coexists with its vapour at P0 --- the heat of liquefaction is
// reproduced by construction. The occupancy of a site is the minimiser of its grand potential
//
//   omega(theta) = (U - mu) theta + (W/2) theta² + kT [theta ln theta + (1-theta) ln(1-theta)],
//
// a Maxwell construction: below the mean-field critical temperature the site condenses, jumping from a
// Boltzmann-dilute occupancy to a filled one when mu - U crosses W/2. Deep wells fill first, shallow
// corners last. A continuum lattice gas on the energy-grid voxels themselves is deliberately not used:
// in an ultramicroporous channel the region of attractive centres has vanishing volume (a filament of a
// few Å³), so packing dV / v_L undercounts by orders of magnitude while GCMC still places a molecule
// every λ along the channel; the geometric capacities keep that packing.
//
// The chemical potential is anchored so the x -> 0 slope equals the exact grid Henry limit
// (P0 V / kT) <exp(-U/kT)>; the anchor absorbs the vibrational prefactor a bare lattice model misses.
// When a sheet and a filament both exist, the filament uses the sheet's C so a deep core cannot steal the
// Henry match and crush the wall.
//
// The BET number is then a fit to that isotherm over the Rouquerol-consistent window. In a micropore the
// fitted n_m can read more than the geometric monolayer (filling already saturated in the window), which
// is the well-known artifact of BET on ultramicroporous materials and is left in, not corrected away.
//
// The independent-site BET/Langmuir model survives as `fromSamples` for unit tests; `run` uses the
// lattice-gas isotherm on the same geometric sites.

export inline constexpr double nitrogenCrossSection = 16.2;          // Å²
export inline constexpr double nitrogenLiquidVolume = 57.7;          // Å³ / molecule
export inline constexpr double nitrogenLinearSpacing = 3.864107;      // Å, cbrt(v_L), 1-D file along a filament
// How close two nitrogen centres are allowed in the packing that reads the capacity off the field. This is
// the first peak of liquid nitrogen's pair distribution, not v_L^{1/3}: that is the spacing of a simple
// cubic lattice at liquid density, and a greedy packing at that distance leaves gaps a thermal liquid
// fills, undercounting the cages (faujasite 80 against WHAM's 120, beta 16 against 25) while a 10-ring
// file is almost insensitive to the difference. Two centres abreast still will not fit in a channel
// narrower than this, which is the point of counting rather than measuring.
export inline constexpr double nitrogenPackingDistance = 3.5;  // Å
// How far a placed molecule attracts its neighbours: the first minimum of liquid nitrogen's pair
// distribution, so the shell is the first coordination shell and nothing else. The number of neighbours
// that shell holds at liquid density follows from its volume, and it is the coordination at which the
// model has to reproduce the heat of liquefaction; a molecule with fewer neighbours than that is bound by
// proportionately less.
export inline constexpr double nitrogenCouplingDistance = 5.4;  // Å
export inline constexpr double nitrogenLiquidCoordination =
    4.18879020478639 * nitrogenCouplingDistance * nitrogenCouplingDistance * nitrogenCouplingDistance /
    nitrogenLiquidVolume;
export inline constexpr double nitrogenHeatOfLiquefactionKJMol = 5.57;
export inline constexpr double nitrogenBETTemperature = 77.355;      // K
export inline constexpr double nitrogenSaturationPressure = 101325.0;  // Pa

export struct IsothermPoint
{
  double relativePressure{0.0};
  double moleculesPerCell{0.0};
};

// One adsorption site for the lattice-gas isotherm: geometric capacity (molecules) and well-floor energy.
export struct LatticeSite
{
  double capacity{0.0};
  double energy{0.0};
};

// The mean-field lattice-gas isotherm over geometric sites. `henryPrefactor` is the factor the Henry
// anchor applied to the bare lattice activity; `saturationCapacity` is n at x -> 1.
export struct LatticeIsotherm
{
  double henryPrefactor{1.0};
  double saturationCapacity{0.0};
  std::vector<IsothermPoint> isotherm;
};

// Pure and unit-testable: energies, thermalEnergy (kT) and heatOfLiquefaction (q_L) in the same internal
// units; henryMoleculesPerCell the exact x -> 0 anchor (P0 V / kT) <exp(-U/kT)>.
//
// With `neighbours` the attraction acts between sites rather than within them: each site carries one
// molecule, the field it feels is q_L / z_L per occupied neighbour, and the occupancies are the solution
// of the Bragg–Williams equations on that graph, chosen between the empty and the full start by grand
// potential. A molecule with a bulk liquid's worth of neighbours then condenses at P0 exactly as before,
// while one in a file of two does not condense at all --- the mean-field transition temperature falls
// with the coordination, and at z = 2 it is 59 K, below the 77 K the isotherm is measured at. Without a
// graph the older per-site form is used: the site is a patch of wall holding many molecules and the
// coupling is internal to it.
export LatticeIsotherm latticeGasIsotherm(std::span<const LatticeSite> sites, double henryMoleculesPerCell,
                                          double thermalEnergy, double heatOfLiquefaction,
                                          std::span<const std::vector<std::uint32_t>> neighbours = {});

// Contact-sheet and filament packing turned into lattice sites. A 1-D file packs L / v_L^{1/3}, a 2-D
// midplane area / 16.2, a blob (no mesh) dV / v_L. When a sheet is present the filament's energy is the
// sheet's, so a deep core cannot steal the Henry match.
export struct GeometricSites
{
  std::vector<LatticeSite> sites;
  // Which sites are close enough to attract each other, one list per site, and empty when the caller has
  // no geometry to say. A mesh site is a patch of wall holding many molecules and has no such answer; a
  // packed molecule does, and it is what stops the whole of it condensing at one pressure.
  std::vector<std::vector<std::uint32_t>> neighbours;
  double sheetCapacity{0.0};
  double filamentFileCapacity{0.0};
  double filamentMidplaneCapacity{0.0};
  double filamentBlobCapacity{0.0};
  double filamentCapacity{0.0};
};

// `filamentBoundaryArea` is the area of the filament overlay's own surface. The contact sheet is part of
// that same boundary, so the share of it the sheet has drawn is a share of the pore already counted, and
// the filament keeps only the rest. Zero means the caller has no overlay and the filament counts in full.
export GeometricSites geometricAdsorptionSites(std::span<const SheetPatch> patches,
                                               std::span<const FilamentVoxel> filament,
                                               double filamentBoundaryArea = 0.0);

// The same sites read off the field instead of off the meshes drawn from it. The mesh route has to decide
// first what shape the pore is --- a wall to spread a monolayer over, a file to thread, a midplane to tile
// --- and then apply the packing rule for that shape, and it is those rules rather than the field that put
// mordenite 43% below what a simulated isotherm holds and beta 28% above. Filling the attractive region
// deepest-first, refusing to place a molecule within a packing distance of one already placed, needs no
// such decision: a channel too narrow for two abreast admits a single file because the geometry says so,
// and a cage packs at bulk density for the same reason. Every molecule placed is one site carrying the
// energy at the point it sits, so the site energies come out spread as finely as the field is sampled.
export GeometricSites fieldAdsorptionSites(const WellField &field,
                                           double packingDistance = nitrogenPackingDistance);

export struct BETSurfaceArea
{
  double gravimetricArea{0.0};  // m²/g, the BET number
  double volumetricArea{0.0};   // m²/cm³
  double monolayerCapacity{0.0};  // n_m, molecules per unit cell
  double cConstant{0.0};
  double windowLow{0.0};
  double windowHigh{0.0};
  double rSquared{0.0};
  double henryAnchor{1.0};  // f, the prefactor the Henry match imposed
  double sheetCapacity{0.0};              // molecules / cell from the contact sheet
  double sheetGravimetricArea{0.0};       // m²/g of that sheet
  double filamentFileCapacity{0.0};       // molecules / cell, 1-D L / v_L^{1/3}
  double filamentMidplaneCapacity{0.0};   // molecules / cell, 2-D midplane / 16.2 Å²
  double filamentBlobCapacity{0.0};       // molecules / cell, dV / v_L when no mesh
  double filamentCapacity{0.0};           // molecules / cell, sum of the three filament parts
  double fieldCapacity{0.0};              // molecules / cell packed into the attractive field, what the isotherm uses
  double fieldCoordination{0.0};          // mean neighbours of a packed molecule, against z_L in the bulk liquid
  double filamentFileGravimetricArea{0.0};
  double filamentMidplaneGravimetricArea{0.0};
  double multilayerRoom{0.0};     // layers the pore volume can hold; 0 means an open surface
  double unsheetedFilamentShare{0.0};  // share of the filament overlay the contact sheet has not drawn
  double saturationCapacity{0.0};  // molecules / cell the lattice isotherm holds at x -> 1

  // Read off the isotherm rather than fitted to it: the Gurvich volume is the loading where the isotherm
  // has stopped rising, taken as liquid nitrogen, and the t-plot is the same isotherm compared against a
  // non-porous silica reference over 3.2 to 4.2 Å of film thickness.
  double saturationLoading{0.0};       // molecules / cell at the top of the isotherm
  double microporeVolume{0.0};         // mL/g, that loading as liquid nitrogen
  double tPlotMicroporeVolume{0.0};    // mL/g, intercept of the t-plot
  double tPlotExternalArea{0.0};       // m²/g, its slope: mesopore plus external surface
  double tPlotRSquared{0.0};
  std::size_t tPlotNumberOfPoints{0};
  // Set when no window in the isotherm was a BET line and it was read as Type I instead: the monolayer is
  // the plateau and C comes from the Henry slope.
  bool plateauReading{false};

  std::vector<IsothermPoint> isotherm;

  WellSurface surface;
  WellField field;
  double seconds{0.0};

  BETSurfaceArea();
  ~BETSurfaceArea();

  // Pure core, unit-testable without a framework. Energies and thermalEnergy are in the same internal units.
  // `multilayerRoom` is how many monolayers the pore volume can hold; zero means an open surface with the
  // classical unlimited stack. `filamentBoundaryArea` is the area of the filament overlay's own surface,
  // of which the contact sheet is a part; zero lets the filament count in full.
  static BETSurfaceArea fromSamples(std::span<const SheetPatch> patches, std::span<const FilamentVoxel> filament,
                                    double henryMoleculesPerCell, double thermalEnergy, double heatOfLiquefaction,
                                    double mass, double cellVolume, double multilayerRoom = 0.0,
                                    double filamentBoundaryArea = 0.0);

  // The Rouquerol window and BET line fitted to an isotherm someone else produced (the lattice model, a
  // GCMC run, a measurement). Only the fit fields are filled. `crossSection` [Å²] and `liquidVolume`
  // [Å³/molecule] come from the adsorbate Component (e.g. N2: 16.2 and 57.7).
  static BETSurfaceArea fromIsotherm(std::vector<IsothermPoint> isotherm, double mass, double cellVolume,
                                     double crossSection, double liquidVolume);

  // Refit slope and intercept of the BET line inside a fixed Rouquerol window (no window search).
  // Used for block jackknife error bars on the area: the window comes from the full-data fit.
  // On failure (too few points, non-physical line) monolayerCapacity and the areas stay zero.
  static BETSurfaceArea fromIsothermFixedWindow(std::vector<IsothermPoint> isotherm, double mass, double cellVolume,
                                                double windowLow, double windowHigh, double crossSection,
                                                double liquidVolume);

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double level, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, std::span<const BlockingSphere> blockingSpheres = {}, double energyScale = 0.0,
           double blockedEnergyPerAngstrom = 0.0, double ceiling = 0.0, bool useElectrostatics = true,
           double relativePrecision = 1e-6);
};

export void writeBETSurfaceArea(std::ostream &stream, const BETSurfaceArea &bet);

// The properties that come off the isotherm rather than off the BET line: the Gurvich micropore volume and
// the t-plot. Every route that produces an isotherm reports them the same way, so the block is shared.
export void writeIsothermDerivedProperties(std::ostream &stream, const BETSurfaceArea &bet);
