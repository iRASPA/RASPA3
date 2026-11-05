module;

export module energy_shared_nldft;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_well_field;
import energy_shared_bet_surface_area;

// A predicted nitrogen isotherm by classical density functional theory, solved on the framework's own
// energy grid rather than on an idealised slit or cylinder.
//
// The kernel approach of experimental porosimetry exists for the inverse problem: an isotherm is measured,
// the structure is unknown, and a library of slit/cylinder pores is deconvolved against it. Here the
// structure is known exactly, so the forward problem is solved instead: the external potential V_ext(r) is
// the framework's dispersion-plus-Ewald field. Two densities are available:
//
//   --nldft                spherical ρ(r) of a single-site probe; V_ext is the (clipped) field
//   --nldft --molecule N2  ρ(r, ω) of the linear guest; V_ext is U(r, ω) in the ideal-rotator term.
//                          White Bear FMT is a fused dumbbell at the TraPPE N–N σ (not the Ravikovitch
//                          COM sphere). Guest–guest attraction is TraPPE N–N 12-6 mean field at the
//                          dispersing sites (±0.55 Å); guest–guest Coulomb/quadrupole is not in F_att.
//                          U is soft TraPPE (Coulomb gated at r < σ). Spherical --nldft still clips
//                          V_ext at 12 kT and uses Ravikovitch WCA.
//
// A 77 K Helmholtz PMF stuffed into spherical ρ is not the dimer: that path stays the clipped
// spherical option. The molecular path recovers the GCMC Henry ⟨e^{-βU(r,ω)}⟩.
//
// The grand potential is the standard NLDFT functional,
//
//   Omega[rho] = F_id[rho] + F_hs[rho] + F_att[rho] + int rho (V_ext - mu),
//
// with the hard-sphere repulsion by fundamental measure theory in the White Bear form (its uniform limit
// is the Carnahan-Starling fluid) and the attraction in mean field: spherical --nldft uses the Ravikovitch
// WCA tail on the COM; --nldft --molecule N2 uses TraPPE N–N 12-6 from r = σ at the dispersing sites
// (the r < σ core is fused-dumbbell FMT at that σ). The Euler-Lagrange equation
//
//   rho(r) = rho_b exp[ -beta ( V_ext + mu_hs[rho](r) - mu_hs,b + int rho u_att - rho_b a ) ]
//
// is solved by descending Omega with L-BFGS over ln rho, every step held to lowering it: over the density at each voxel
// for a spherical density, and for ρ(r, ω) over what each orientation holds at each voxel, so that the orientations are
// those of the density itself and the hard-sphere repulsion is all that keeps either below close packing. Only a
// converged profile is the warm start for the next pressure, and it is all that is carried there. All convolutions (the
// four FMT weighted densities and the attractive tail) are products in Fourier space, so each evaluation is a fixed
// number of FFTs over the periodic cell, for any cell shape. The FMT weights are non-negative kernels laid onto the
// grid in real space, which keeps n3, n2 and n2 - |n2v| from going negative on any grid, and the third term of the free
// energy is the Rosenfeld-Tarazona one, so the free energy density is non-negative and the functional bounded below.
//
// The bulk fluid the pore is in equilibrium with is the same functional at uniform density: its gas-liquid
// coexistence at the run temperature is found by a Maxwell construction once, and the relative pressure
// x = P / P0 is measured against that model P0, exactly as kernel builders do (the model fluid's P0, not
// the experimental table value, is what makes x -> 1 mean bulk condensation inside the model). Layering,
// micropore filling and capillary condensation then emerge from the real geometry; nothing is synthesised
// from sites or packings. For spherical ρ the x -> 0 limit is rho_b exp(-beta U). For ρ(r, ω) it is
// rho_b ⟨exp(-beta U(r,ω))⟩_ω, the GCMC Henry coefficient on the same orientations.
//
// The pressure is swept up, warm-starting each point from the one below, and where that fills a pore it is swept back
// down from the filled profile. Across a hysteresis loop both are metastable in opposite directions, so the isotherm
// reported is the equilibrium one: at each pressure, whichever has the lower grand potential Omega, the same
// construction the bulk coexistence is found by. Below the loop that is the sweep up and above it the sweep down.
//
// The BET number reported alongside is a fit to this isotherm over the Rouquerol-consistent window, by the
// criteria a measured or GCMC isotherm is read by, and off the curve through the pressures solved at rather
// than off those pressures, so that it does not turn on where they fall.

// N2 at 77 K (Ravikovitch, Vishnyakov, Russo and Neimark, Langmuir 16, 2311 (2000)).
export inline constexpr double nldftSigma = 3.575;               // Å
export inline constexpr double nldftEpsilonKelvin = 94.45;       // K
export inline constexpr double nldftHardSphereDiameter = 3.575;  // Å
export inline constexpr double nldftCutoffFactor = 5.0;          // r_c = 5 sigma
export inline constexpr double nldftWellFloorInKT = 12.0;        // clip V_ext below -this * T

export struct NLDFTOptions
{
  double temperature{77.355};                        // K
  double sigma{nldftSigma};                          // Å
  double epsilon{nldftEpsilonKelvin};                // K
  double hardSphereDiameter{nldftHardSphereDiameter};  // Å
  double cutoff{nldftCutoffFactor * nldftSigma};     // Å

  int pressurePoints{80};
  // Defaults match the old slit-kernel habit (x from 1e-5 to 0.35). The CLI / report path overwrites
  // them with nldftApplyHenryPressureWindow so ultramicroporous Type I fills are inside the solve,
  // the same policy WHAM and TMMC use for their pressure grids.
  double xLow{1.0e-5};
  double xHigh{0.35};

  // RMS |rho_EL - rho| against the RMS density, at which loadings have settled to about five figures. For ρ(r, ω)
  // both are over the orientations at each voxel, not over their sum.
  double tolerance{1.0e-4};
  // The descent also stops once the grand potential and the loading have both moved, over a hundred of its steps, by
  // less than this fraction of the loading: the tail of a pore packed close to solid, where the residual is still
  // coming down long after neither of them shows it. Zero is a tenth of the tolerance.
  double steadyTolerance{0.0};
  std::size_t maxIterations{2000};  // per pressure point (warm-started from the last *converged* profile)
  bool seedLiquid{false};           // start from the bulk liquid (empty-box coexistence test); default is Henry

  // Attractive wells deeper than this (in units of kT) make exp(-βV_ext) overflow the packing
  // constraint, so the Euler-Lagrange equation has no solution and Rouquerol sees a step at x = 0.
  // Twelve kT is a physical N2 heat of adsorption; a spherical density cannot occupy a deeper
  // orientational well, which is the extra depth of a 77 K molecular PMF.
  double wellFloorInKT{nldftWellFloorInKT};

  // Orientation-dependent (dumbbell) FMT for ρ(r, ω). Offsets are along the molecular axis in Å,
  // diameter is the dispersing site's σ. Site weights are V_union / (n_sites V_sphere) so the
  // uniform packing is the fused dumbbell, not n_sites independent spheres. Empty offsets keep
  // spherical FMT on the COM density.
  std::vector<double> dumbbellSiteOffsets{};
  double dumbbellSiteDiameter{0.0};
  bool headTailSymmetric{true};

  // Guest–guest attraction on dispersing sites (TraPPE N–N). Offsets along the axis in Å; σ and ε/k
  // in Å and K. Empty offsets keep the Ravikovitch COM WCA tail. The pair is the full 12-6 from
  // contact (r = σ) to the cutoff, not the WCA split.
  std::vector<double> attractiveSiteOffsets{};
  double attractiveSigma{0.0};
  double attractiveEpsilon{0.0};

  // The energy grid is far finer than the density profile needs: FMT's contact peaks are ~0.2 Å wide,
  // while the field is built at ~0.06-0.15 Å. Before solving, each axis is decimated by the largest
  // divisor of its point count that keeps the spacing at or under this, which is worth an order of
  // magnitude in the FFTs. Zero switches the decimation off.
  double targetSpacing{0.25};  // Å
};

// The model fluid's own bulk coexistence at the run temperature, by Maxwell construction on the uniform
// limit of the functional. Densities in molecules / Å³, pressure in Pa. If the uniform limit has no
// van der Waals loop, x is measured against experimental N2 P0 while μ and P stay the same EOS.
export struct NLDFTBulk
{
  double saturationPressure{0.0};
  double gasDensity{0.0};
  double liquidDensity{0.0};
  double meanFieldIntegral{0.0};  // a = int u_att d³r, in K Å³ (negative)
  bool experimentalSaturation{false};
};

export NLDFTBulk nldftBulkCoexistence(const NLDFTOptions &options);

// The gas-branch bulk density in equilibrium with relative pressure x = P / P0(model), molecules / Å³.
export double nldftGasDensity(const NLDFTOptions &options, const NLDFTBulk &bulk, double x);

// Same pressure-window policy as WHAM / TMMC: start where the Henry line holds
// `nldftPressureWindowFillFraction` of the liquid-packing capacity V_cell / v_L, never above
// `nldftPressureWindowFloorPa`, and end at the isotherm's P0 (xHigh = 1).
// `henryMoleculesPerCellAtExperimentalP0` is the grid ⟨e^{-βU}⟩ occupancy at experimental P0
// (NLDFTIsotherm::henryMoleculesPerCell). `isothermSaturationPressurePa` is the P0 used for x.
export inline constexpr double nldftPressureWindowFillFraction = 1.0e-3;
export inline constexpr double nldftPressureWindowFloorPa = 1.0;
export void nldftApplyHenryPressureWindow(NLDFTOptions &options, double henryMoleculesPerCellAtExperimentalP0,
                                          double cellVolumeAngstrom3, double isothermSaturationPressurePa);

// Pure and unit-testable: energies in Kelvin on the periodic grid (x varying fastest), any cell shape.
// One point of a sweep, with what says how stable it is. The grand potential is of the profile at this pressure, in
// units of kT per unit cell, and comparing it between branches at one pressure says which is the stable state. It
// means that only for a profile that solved the equation at this pressure: one that ran out of iterations, or that
// was held from a lower pressure, is not a state at this pressure at all, and `stationary` says which is which.
export struct NLDFTBranchPoint
{
  double relativePressure{0.0};
  double moleculesPerCell{0.0};
  double grandPotential{0.0};      // kT / cell
  double chemicalPotential{0.0};   // mu / kT of the reservoir, on the same zero as the grand potential
  double residual{0.0};            // how far the profile was from the equation, against the density it is measured by
  bool stationary{false};
};

export struct NLDFTGridIsotherm
{
  NLDFTBulk bulk;
  // Molecules / cell against x = P / P0(model): the equilibrium isotherm, which is the branch of lower grand potential
  // at each pressure, and the one the BET fit is read from. Where a pore's two branches differ neither of them is the
  // stable state throughout, the sweep up staying empty too long and the sweep down staying full too long.
  std::vector<IsothermPoint> isotherm;
  // The two branches themselves, at the same pressures as `isotherm`. Each point of the sweep up is warm-started from
  // the one below it and each point of the sweep down from the one above it. Where the sweep up never fills a pore
  // nothing is swept back down, and the two are the same.
  std::vector<NLDFTBranchPoint> adsorptionIsotherm;
  std::vector<NLDFTBranchPoint> desorptionIsotherm;
  // Where the equilibrium isotherm leaves the sweep up for the sweep down: the lowest pressure at which the filled
  // state is the stable one. Zero where the two branches never differ.
  double equilibriumTransition{0.0};
  // Pressures at which neither branch has a profile that solves the equation, where the isotherm carries the sweep up's
  // loading for want of anything better. A count above zero says that what is read off the isotherm there, the BET
  // fit included if its window reaches them, is not a result of the functional.
  std::size_t pointsWithoutASolution{0};
  std::size_t unconverged{0};           // pressure points that ran out of iterations
  bool molecular{false};                // ρ(r, ω) for a linear guest; otherwise spherical ρ(r)
  std::size_t numberOfOrientations{1};
};

export NLDFTGridIsotherm nldftIsothermOnGrid(std::span<const float> energyKelvin, uint3 gridSize,
                                             const UnitCell &cell, const NLDFTOptions &options,
                                             std::span<const float> orientationKelvin = {},
                                             std::size_t nOrientations = 1);

// Attractive wells deeper than -wellFloorInKT * T fill at x = 0 for a spherical density. Clip them so
// a leftover Coulomb hole cannot set the Henry coefficient. Do not use well-field distance < 0 as a
// wall mask: that signed distance is clearance to rmin of the linear probe, which is negative through
// most of a zeolite channel even where a spherical COM still fits. `energyKelvin` is overwritten.
export struct NLDFTExternalFieldStats
{
  std::size_t clippedVoxels{0};
  double wellFloorKelvin{0.0};
};

export NLDFTExternalFieldStats nldftMaskAndClipExternalPotential(std::span<float> energyKelvin,
                                                                double temperature, float ceilingKelvin,
                                                                double wellFloorInKT = nldftWellFloorInKT);

export struct NLDFTIsotherm
{
  NLDFTBulk bulk;
  std::vector<IsothermPoint> isotherm;                // see NLDFTGridIsotherm, as for the next three
  std::vector<NLDFTBranchPoint> adsorptionIsotherm;
  std::vector<NLDFTBranchPoint> desorptionIsotherm;
  double equilibriumTransition{0.0};
  std::size_t pointsWithoutASolution{0};
  std::size_t unconverged{0};

  // The Rouquerol/BET read of the isotherm (only the fit fields of it are filled).
  BETSurfaceArea fit;

  double henryMoleculesPerCell{0.0};  // exact grid limit, the x -> 0 consistency anchor
  WellField field;
  double seconds{0.0};
  bool molecular{false};
  std::size_t numberOfOrientations{1};

  NLDFTIsotherm();
  ~NLDFTIsotherm();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations, double temperature,
           std::span<const BlockingSphere> blockingSpheres = {}, bool useElectrostatics = true,
           double relativePrecision = 1e-6);
};
