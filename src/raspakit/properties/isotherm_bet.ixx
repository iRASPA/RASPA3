module;

export module isotherm_bet;

import std;

import json;
import system;
import component;

export import energy_shared_bet_surface_area;

/// BET (Rouquerol) extraction from a reweighted isotherm. Probe cross-section, liquid volume and
/// saturation pressure P0 come from the adsorbate Component ('CrossSection', 'LiquidVolume',
/// 'SaturationPressure'). Used by WHAM and parallel-TMMC when 'ComputeBET' is set.

/// One (pressure [Pa], loading [molecules / crystallographic unit cell]) point of a simulated isotherm.
export struct SimulatedIsothermPoint
{
  double pressure{0.0};
  double moleculesPerCell{0.0};
};

/// Probe properties for BET area, Gurvich packing and relative pressure x = P / P0.
export struct BETProbeProperties
{
  double crossSection{0.0};        ///< [Å²]
  double liquidVolume{0.0};        ///< [Å³ / molecule]
  double saturationPressure{0.0};  ///< P0 [Pa]
};

/// Read CrossSection, LiquidVolume and SaturationPressure from a component; throws if any is missing or ≤ 0.
export BETProbeProperties requireBETProbeProperties(const Component& component);

/// Convert an isotherm in (P [Pa], n [molecules / cell]) to the relative-pressure form the Rouquerol
/// fit wants, dropping the endpoints x <= 0 and x >= 1 that the BET plot is undefined at.
export std::vector<IsothermPoint> nitrogenBETIsotherm(std::span<const SimulatedIsothermPoint> points,
                                                      double saturationPressure);

/// Fit the BET line, Gurvich volume and t-plot to an isotherm. `mass` is the crystallographic
/// unit-cell mass [g/mol] and `cellVolume` its volume [Å³].
export BETSurfaceArea fitNitrogenBET(std::span<const SimulatedIsothermPoint> points, double mass, double cellVolume,
                                     const BETProbeProperties& probe);

/// Refit slope/intercept inside a fixed Rouquerol window (no window search). Returns zero area on failure.
export BETSurfaceArea fitNitrogenBETFixedWindow(std::span<const SimulatedIsothermPoint> points, double mass,
                                                double cellVolume, double windowLow, double windowHigh,
                                                const BETProbeProperties& probe);

/// Write the Rouquerol/BET numbers (area, n_m, C, window, Gurvich, t-plot) to a report stream.
/// Optional block confidence-interval errors print next to the area / n_m / C when provided.
export void writeNitrogenBETSummary(std::ostream& stream, const BETSurfaceArea& bet, const BETProbeProperties& probe,
                                    std::string_view indent = "",
                                    std::optional<double> gravimetricAreaError = std::nullopt,
                                    std::optional<double> monolayerCapacityError = std::nullopt,
                                    std::optional<double> cConstantError = std::nullopt);

/// Write the BET-plot table (x, n, n(1-x), x/(n(1-x))) that the fit is read from, marking the window.
export void writeNitrogenBETTable(std::ostream& stream, const BETSurfaceArea& bet,
                                  const BETProbeProperties& probe);

/// JSON object of the fitted numbers (no isotherm table).
export nlohmann::json nitrogenBETJson(const BETSurfaceArea& bet, const BETProbeProperties& probe);

/// Number of Widom insertions used to place the nitrogen BET pressure span (order of magnitude of K is enough).
export constexpr std::size_t nitrogenBETHenryInsertions = 200000;

/// Henry-based nitrogen BET pressure placement: sampling ladder (WHAM) and reweighting span (WHAM and TMMC).
export struct NitrogenBETPressurePlan
{
  double henryCoefficientPerCell{0.0};  ///< Molecules per crystallographic unit cell per Pa.
  double packingCapacity{0.0};          ///< Gurvich packing n = V_cell / v_L [molecules / cell].
  double lowestPressure{1.0};           ///< min(1 Pa, 10^{-3} n_Gurvich / K_H).
  double highestPressure{0.0};  ///< Saturation pressure P0 from the component [Pa].
  std::vector<double> samplingPressures{};             ///< Log-spaced WHAM rungs from lowest to P0.
  std::size_t whamReweightingNumberOfPressures{400};   ///< Dense WHAM isotherm grid.
  std::size_t tmmcReweightingNumberOfPressures{100};   ///< TMMC grid: 12 points/decade, at least 100.
  BETProbeProperties probe{};                          ///< Cross-section, liquid volume, P0 from Component.
};

/// Place the nitrogen BET pressure span from a Widom Henry coefficient of the empty framework.
/// `numberOfThreads` sizes the WHAM sampling ladder (8–16 rungs), matching raspa3-cli `--threads`.
export NitrogenBETPressurePlan planNitrogenBETPressures(System system, double temperature,
                                                        std::size_t numberOfThreads);

/// Write the Henry coefficient and chosen pressure span to a report stream.
export void writeNitrogenBETPressurePlan(std::ostream& stream, const NitrogenBETPressurePlan& plan);

/// Re-place a pressure ladder so neighboring rungs have equal Δn.
///
/// `pressures` must be strictly positive and sorted. Occupancies are made non-decreasing, a
/// monotonic PCHIP is fit to n(log P), and new pressures are the inverse at uniform n between
/// the measured min and max. Endpoints are kept. A flat isotherm returns the original ladder.
export std::vector<double> placePressuresByEqualOccupancy(std::span<const double> pressures,
                                                          std::span<const double> occupancies);

/// Neighboring auto-WHAM rungs may not be farther apart in log P than this many
/// equal-log steps: Δlog P ≤ multiplier × (log P0 − log P_Henry) / (N − 1).
export constexpr double nitrogenBETMaxLogSpacingMultiplier = 2.0;
/// Extra rungs inserted into a log-skeleton interval of the auto WHAM ladder.
export constexpr std::size_t nitrogenBETMaxGradientExtrasPerInterval = 2;
/// Langmuir coverages of the short GCMC probes used to fit Langmuir–Freundlich / Toth.
export constexpr std::array<double, 3> nitrogenBETLangmuirProbeCoverages{0.1, 0.5, 0.9};

/// Type I sketch used to place the auto WHAM pressure ladder.
export enum class NitrogenBETIsothermModel
{
  Langmuir,              ///< m = 1, b = K_H / n_sat.
  LangmuirFreundlich,    ///< Sips: θ = (bP)^m / (1 + (bP)^m).
  Toth                   ///< θ = bP / (1 + (bP)^t)^{1/t}, b = K_H / n_sat.
};

/// One short GCMC occupancy scout (P0 ceiling or a Langmuir-coverage probe).
export struct NitrogenBETScoutPoint
{
  double pressure{0.0};          ///< Scout pressure [Pa].
  double meanOccupancy{0.0};     ///< Plateau mean [molecules in the simulation cell].
  double standardDeviation{0.0}; ///< Plateau standard deviation.
  std::size_t peakOccupancy{0};  ///< Highest occupancy seen.
  std::size_t cycles{0};         ///< Scout cycles run.
  bool converged{false};         ///< True when the loading plateaued before the cycle cap.
  std::size_t gurvichCap{0};     ///< Hard upper bound (helium-void Gurvich, else V / v_L).
};

/// Fitted Type I isotherm from K_H, the P0 plateau, and the Langmuir-coverage probes.
export struct NitrogenBETPreIsothermFit
{
  NitrogenBETIsothermModel model{NitrogenBETIsothermModel::Langmuir};
  double saturationLoadingPerCell{0.0};  ///< n_sat [molecules / cell].
  double affinity{0.0};                  ///< Langmuir / Sips / Toth b [1/Pa].
  double exponent{1.0};                  ///< Sips m or Toth t (1 for Langmuir).
  double rSquared{0.0};                  ///< Linearized Sips r² (or Toth n-space r²).
  double langmuirRSquared{0.0};          ///< Linearized Sips r² with m fixed at 1.
  bool langmuirRejected{false};          ///< Nested test rejected m = 1.
  bool hillPlotFlat{true};               ///< Consecutive pairwise Sips slopes agree.
  double hillExponentLow{1.0};           ///< Pairwise m on the lower two valid probes.
  double hillExponentHigh{1.0};          ///< Pairwise m on the upper two valid probes.
  const char* selectionReason{"Langmuir from K_H and n_sat"};
};

/// Unbiased GCMC occupancy at a set pressure. The input system is copied.
export NitrogenBETScoutPoint scoutNitrogenBETOccupancy(System system, double pressurePa);

/// Langmuir P(θ) = (θ / (1 − θ)) / b. Returns 0 if b or θ is not usable.
export double langmuirPressureAtCoverage(double coverage, double affinity);

/// Choose Langmuir vs Sips vs Toth from K_H, n_sat, and the Langmuir-coverage probes.
/// Probe loadings are in molecules per crystallographic unit cell.
export NitrogenBETPreIsothermFit fitNitrogenBETPreIsotherm(double henryCoefficientPerCell,
                                                           double saturationLoadingPerCell,
                                                           double gurvichLoadingPerCell,
                                                           std::span<const SimulatedIsothermPoint> probes);

/// n(P) [molecules / cell] from a pre-isotherm fit.
export double nitrogenBETPreIsothermLoadingPerCell(const NitrogenBETPreIsothermFit& fit, double pressure);

/// Place WHAM rungs from a fitted Type I isotherm: equal Fisher overlap (arcsine in θ
/// for Langmuir/Sips; ∫ dθ/√(θ(1−θ^t)) for Toth) on a log-spaced skeleton (gap ≤
/// `nitrogenBETMaxLogSpacingMultiplier` equal-log steps). Leftover rungs go into the
/// skeleton intervals with the largest Fisher-coordinate jump, at most
/// `nitrogenBETMaxGradientExtrasPerInterval` per interval. Endpoints are kept.
export std::vector<double> placePressuresByPreIsothermFisher(double lowestPressure, double highestPressure,
                                                             std::size_t numberOfRungs,
                                                             const NitrogenBETPreIsothermFit& fit);

/// Place WHAM rungs from a Langmuir sketch of K_H and n_sat (Fisher overlap, log skeleton).
/// A missing Henry coefficient, a non-positive plateau, or a flat window returns a
/// log-spaced ladder.
export std::vector<double> placePressuresByLangmuirEqualOccupancy(double lowestPressure, double highestPressure,
                                                                  std::size_t numberOfRungs,
                                                                  double henryCoefficientPerCell,
                                                                  double saturationLoadingPerCell);

/// Pin T, P and Peng–Robinson fugacity coefficients on a system (WHAM replica re-placement).
export void pinSystemPengRobinsonPressure(System& system, double pressurePa);

/// Unbiased GCMC scout at the component saturation pressure P0: filling ceiling for TMMC N_max / WHAM cycles.
export struct NitrogenBETFillingCeiling
{
  std::size_t maxMacrostate{0};   ///< Ceiling placed just above the plateau loading.
  double meanOccupancy{0.0};      ///< Plateau mean occupancy [molecules in the simulation cell].
  double standardDeviation{0.0};  ///< Plateau standard deviation.
  std::size_t peakOccupancy{0};   ///< Highest occupancy seen (not used to place the ceiling).
  std::size_t cycles{0};          ///< Scout cycles run.
  bool converged{false};          ///< True when the loading plateaued before the cycle cap.
  std::size_t gurvichCap{0};      ///< Hard upper bound (helium-void Gurvich, else V / v_L).
  double saturationPressure{0.0}; ///< P0 used for the scout [Pa].
};

/// Scout occupancy at P0 with the system's own swap moves and return N_max.
/// The input system is copied; TMMC window bounds on the original are not used (the scout must
/// be free to fill past the default maxMacrostate of 100).
export NitrogenBETFillingCeiling scoutNitrogenBETFillingCeiling(System system);

/// Write the P0 scout and chosen filling ceiling to a report stream.
export void writeNitrogenBETFillingCeiling(std::ostream& stream, const NitrogenBETFillingCeiling& ceiling);

/// Liquid packing of a volume [Å³]: ceil(V / liquidVolume). Used as a hard cap on the P0 scout.
export inline std::size_t gurvichOccupancy(double volume, double liquidVolume)
{
  return std::max(1uz, static_cast<std::size_t>(std::ceil(volume / liquidVolume)));
}
