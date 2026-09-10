module;

module isotherm_bet;

import std;

import energy_shared_bet_surface_area;
import json;
import int3;
import double3;
import randomnumbers;
import mc_moves;
import mc_moves_widom;
import property_lambda_probability_histogram;
import equation_of_states;
import component;
import system;
import units;

BETProbeProperties requireBETProbeProperties(const Component& component)
{
  if (!component.crossSection.has_value() || !(component.crossSection.value() > 0.0))
  {
    throw std::runtime_error(
        std::format("[ComputeBET]: component '{}' needs a positive 'CrossSection' [Å²]\n", component.name));
  }
  if (!component.liquidVolume.has_value() || !(component.liquidVolume.value() > 0.0))
  {
    throw std::runtime_error(
        std::format("[ComputeBET]: component '{}' needs a positive 'LiquidVolume' [Å³/molecule]\n", component.name));
  }
  if (!component.saturationPressure.has_value() || !(component.saturationPressure.value() > 0.0))
  {
    throw std::runtime_error(std::format(
        "[ComputeBET]: component '{}' needs a positive 'SaturationPressure' P0 [Pa]\n", component.name));
  }
  return BETProbeProperties{.crossSection = component.crossSection.value(),
                            .liquidVolume = component.liquidVolume.value(),
                            .saturationPressure = component.saturationPressure.value()};
}

std::vector<IsothermPoint> nitrogenBETIsotherm(std::span<const SimulatedIsothermPoint> points,
                                               double saturationPressure)
{
  std::vector<IsothermPoint> isotherm;
  isotherm.reserve(points.size());
  for (const SimulatedIsothermPoint& point : points)
  {
    const double x = point.pressure / saturationPressure;
    if (x > 0.0 && x < 1.0)
    {
      isotherm.push_back(IsothermPoint{x, point.moleculesPerCell});
    }
  }
  return isotherm;
}

BETSurfaceArea fitNitrogenBET(std::span<const SimulatedIsothermPoint> points, double mass, double cellVolume,
                              const BETProbeProperties& probe)
{
  return BETSurfaceArea::fromIsotherm(nitrogenBETIsotherm(points, probe.saturationPressure), mass, cellVolume,
                                      probe.crossSection, probe.liquidVolume);
}

BETSurfaceArea fitNitrogenBETFixedWindow(std::span<const SimulatedIsothermPoint> points, double mass, double cellVolume,
                                         double windowLow, double windowHigh, const BETProbeProperties& probe)
{
  return BETSurfaceArea::fromIsothermFixedWindow(nitrogenBETIsotherm(points, probe.saturationPressure), mass,
                                                 cellVolume, windowLow, windowHigh, probe.crossSection,
                                                 probe.liquidVolume);
}

void writeNitrogenBETSummary(std::ostream& stream, const BETSurfaceArea& bet, const BETProbeProperties& probe,
                             std::string_view indent, std::optional<double> gravimetricAreaError,
                             std::optional<double> monolayerCapacityError, std::optional<double> cConstantError)
{
  std::print(stream, "{}BET (Rouquerol consistency, P0 = {:.0f} Pa, σ = {:.1f} Å²)\n", indent,
             probe.saturationPressure, probe.crossSection);
  if (bet.plateauReading)
  {
    std::print(stream, "{}    no linear BET window: Type I reading (monolayer = plateau, C from Henry slope)\n",
               indent);
  }
  if (gravimetricAreaError.has_value())
  {
    std::print(stream, "{}    BET area:              {:14.4f} ± {:.4f} [m²/g]   {:12.4f} [m²/cm³]\n", indent,
               bet.gravimetricArea, *gravimetricAreaError, bet.volumetricArea);
  }
  else
  {
    std::print(stream, "{}    BET area:              {:14.4f} [m²/g]   {:12.4f} [m²/cm³]\n", indent, bet.gravimetricArea,
               bet.volumetricArea);
  }
  if (monolayerCapacityError.has_value())
  {
    std::print(stream, "{}    monolayer n_m:         {:14.4f} ± {:.4f} [molecules / cell]\n", indent,
               bet.monolayerCapacity, *monolayerCapacityError);
  }
  else
  {
    std::print(stream, "{}    monolayer n_m:         {:14.4f} [molecules / cell]\n", indent, bet.monolayerCapacity);
  }
  if (cConstantError.has_value())
  {
    std::print(stream, "{}    C constant:            {:14.4f} ± {:.4f} [-]\n", indent, bet.cConstant, *cConstantError);
  }
  else
  {
    std::print(stream, "{}    C constant:            {:14.4f} [-]\n", indent, bet.cConstant);
  }
  std::print(stream, "{}    fit window:            {:.4g} -- {:.4g} in P/P0, r² = {:.5f}\n", indent, bet.windowLow,
             bet.windowHigh, bet.rSquared);
  if (bet.microporeVolume > 0.0)
  {
    std::print(stream, "{}    micropore volume:      {:14.4f} [mL/g] (Gurvich, {:.4f} molecules / cell as liquid)\n",
               indent, bet.microporeVolume, bet.saturationLoading);
    if (bet.tPlotNumberOfPoints >= 3)
    {
      std::print(stream, "{}    t-plot micropore vol.: {:14.4f} [mL/g]\n", indent, bet.tPlotMicroporeVolume);
      std::print(stream, "{}    t-plot meso+external:  {:14.4f} [m²/g] ({} points, r² = {:.5f})\n", indent,
                 bet.tPlotExternalArea, bet.tPlotNumberOfPoints, bet.tPlotRSquared);
    }
  }
}

void writeNitrogenBETTable(std::ostream& stream, const BETSurfaceArea& bet, const BETProbeProperties& probe)
{
  std::print(stream, "# BET plot of the simulated isotherm (P0 = {:.0f} Pa, σ = {:.1f} Å²)\n",
             probe.saturationPressure, probe.crossSection);
  if (bet.plateauReading)
  {
    std::print(stream, "# Type I reading: no linear BET window, monolayer = plateau, C from Henry slope\n");
  }
  else
  {
    std::print(stream, "# Rouquerol window: {:.6g} -- {:.6g} in P/P0, r² = {:.5f}\n", bet.windowLow, bet.windowHigh,
               bet.rSquared);
  }
  std::print(stream, "# column 5 is 1 inside the fitted window (the points the BET line is read from)\n");
  std::print(stream, "# {:>12} {:>16} {:>16} {:>16} {:>8}\n", "x = P/P0", "n [molec/cell]", "n(1-x)", "x/(n(1-x))",
             "window");
  for (const IsothermPoint& point : bet.isotherm)
  {
    const double g = point.moleculesPerCell * (1.0 - point.relativePressure);
    const double y = (g > 0.0) ? point.relativePressure / g : 0.0;
    const int inWindow =
        (!bet.plateauReading && point.relativePressure >= bet.windowLow && point.relativePressure <= bet.windowHigh)
            ? 1
            : 0;
    std::print(stream, "  {:>12.6g} {:>16.8g} {:>16.8g} {:>16.8g} {:>8d}\n", point.relativePressure,
               point.moleculesPerCell, g, y, inWindow);
  }
}

nlohmann::json nitrogenBETJson(const BETSurfaceArea& bet, const BETProbeProperties& probe)
{
  return nlohmann::json{{"gravimetricArea", bet.gravimetricArea},
                        {"volumetricArea", bet.volumetricArea},
                        {"monolayerCapacity", bet.monolayerCapacity},
                        {"cConstant", bet.cConstant},
                        {"windowLow", bet.windowLow},
                        {"windowHigh", bet.windowHigh},
                        {"rSquared", bet.rSquared},
                        {"plateauReading", bet.plateauReading},
                        {"saturationLoading", bet.saturationLoading},
                        {"microporeVolume", bet.microporeVolume},
                        {"tPlotMicroporeVolume", bet.tPlotMicroporeVolume},
                        {"tPlotExternalArea", bet.tPlotExternalArea},
                        {"tPlotRSquared", bet.tPlotRSquared},
                        {"tPlotNumberOfPoints", bet.tPlotNumberOfPoints},
                        {"saturationPressure", probe.saturationPressure},
                        {"crossSection", probe.crossSection},
                        {"liquidVolume", probe.liquidVolume}};
}

namespace
{
constexpr double ladderBottomFillFraction = 1.0e-3;
constexpr double ladderSpacing = 5.0;
constexpr std::size_t minimumNumberOfPressures = 8;
constexpr std::size_t maximumNumberOfPressures = 16;
constexpr std::size_t reweightingPressurePointsPerDecade = 12;
constexpr std::size_t minimumNumberOfTMMCReweightingPressures = 100;
constexpr std::size_t numberOfWHAMReweightingPressures = 400;

double henryCoefficientPerCell(System& system, double temperature, double cells)
{
  RandomNumber rng(std::size_t{1});
  system.forceField.initializeAutomaticCutOff(system.simulationBox);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  std::ostringstream discard;
  system.createExternalFieldInterpolationGrid(discard, 0);
  system.createFrameworkInterpolationGrids(discard);
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();

  double sum = 0.0;
  for (std::size_t insertion = 0uz; insertion < nitrogenBETHenryInsertions; ++insertion)
  {
    sum += MC_Moves::WidomMove(rng, system, 0uz);
  }
  const double averageRosenbluthWeight = sum / static_cast<double>(nitrogenBETHenryInsertions);

  return averageRosenbluthWeight * Units::AvogadroConstant * system.simulationBox.volume * Units::LengthUnit *
         Units::LengthUnit * Units::LengthUnit / (Units::MolarGasConstant * temperature * cells);
}

std::vector<double> pressureLadder(double lowestPressure, double highestPressure, std::size_t maximumThreads)
{
  const double logMin = std::log(lowestPressure);
  const double logMax = std::log(highestPressure);
  const std::size_t physicalRungs =
      std::clamp(static_cast<std::size_t>(std::ceil((logMax - logMin) / std::log(ladderSpacing))) + 1uz,
                 minimumNumberOfPressures, maximumNumberOfPressures);
  const std::size_t rungs =
      maximumThreads >= minimumNumberOfPressures
          ? std::clamp(maximumThreads, minimumNumberOfPressures, maximumNumberOfPressures)
          : physicalRungs;

  std::vector<double> pressures(rungs);
  for (std::size_t i = 0; i < rungs; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(rungs - 1);
    pressures[i] = std::exp(logMin + t * (logMax - logMin));
  }
  pressures.back() = highestPressure;
  return pressures;
}
}  // namespace

NitrogenBETPressurePlan planNitrogenBETPressures(System system, double temperature, std::size_t numberOfThreads)
{
  if (!system.framework.has_value())
  {
    throw std::runtime_error(
        "[ComputeBET]: automatic pressure placement needs a framework (unit-cell volume for the "
        "Gurvich packing and a host for the Henry coefficient)\n");
  }
  if (system.components.empty())
  {
    throw std::runtime_error("[ComputeBET]: automatic pressure placement needs an adsorbate component\n");
  }

  const BETProbeProperties probe = requireBETProbeProperties(system.components.front());

  const int3 numberOfUnitCells = system.framework->numberOfUnitCells;
  const double cells =
      static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
  const double unitCellVolume = system.simulationBox.volume / cells;

  std::print("  Computing Henry coefficient ({} Widom insertions) for BET pressure placement...\n",
             nitrogenBETHenryInsertions);
  std::cout << std::flush;

  NitrogenBETPressurePlan plan;
  plan.probe = probe;
  plan.packingCapacity = unitCellVolume / probe.liquidVolume;
  plan.henryCoefficientPerCell = henryCoefficientPerCell(system, temperature, cells);
  plan.highestPressure = probe.saturationPressure;
  plan.lowestPressure = 1.0;
  if (plan.henryCoefficientPerCell > 0.0)
  {
    plan.lowestPressure =
        std::min(1.0, ladderBottomFillFraction * plan.packingCapacity / plan.henryCoefficientPerCell);
  }

  double samplingHighest = std::max(plan.highestPressure, plan.lowestPressure * ladderSpacing);
  plan.samplingPressures = pressureLadder(plan.lowestPressure, samplingHighest, numberOfThreads);
  plan.whamReweightingNumberOfPressures = numberOfWHAMReweightingPressures;
  const double reweightingDecades = std::log10(plan.highestPressure / plan.lowestPressure);
  plan.tmmcReweightingNumberOfPressures = std::max(
      minimumNumberOfTMMCReweightingPressures,
      static_cast<std::size_t>(reweightingPressurePointsPerDecade * std::ceil(reweightingDecades)));

  std::print("  Henry coefficient {:.4e} molecules/cell/Pa, pressures {:.4g} -- {:.4g} Pa\n",
             plan.henryCoefficientPerCell, plan.lowestPressure, plan.highestPressure);
  std::cout << std::flush;
  return plan;
}

void writeNitrogenBETPressurePlan(std::ostream& stream, const NitrogenBETPressurePlan& plan)
{
  std::print(stream, "BET pressure placement\n");
  std::print(stream, "    Henry coefficient:                     {:.6e} [molecules/cell/Pa] ({} Widom insertions)\n",
             plan.henryCoefficientPerCell, nitrogenBETHenryInsertions);
  std::print(stream, "    Cross-section σ:                       {:.2f} [Å²]\n", plan.probe.crossSection);
  std::print(stream, "    Liquid volume v_L:                     {:.2f} [Å³/molecule]\n", plan.probe.liquidVolume);
  std::print(stream, "    Saturation pressure P0:                {:.5e} [Pa]\n", plan.probe.saturationPressure);
  std::print(stream, "    Gurvich packing n_L:                   {:.4f} [molecules / cell]\n", plan.packingCapacity);
  std::print(stream, "    Pressure span:                         {:.5e} - {:.5e} [Pa]\n", plan.lowestPressure,
             plan.highestPressure);
  std::print(stream,
             "    (bottom is min(1 Pa, {:g} of packing / K_H); top is P0. A fixed 1 Pa bottom\n"
             "     puts every replica on the saturated plateau of a strongly binding framework.)\n"
             "    WHAM then scouts occupancy at Langmuir θ = 0.1, 0.5, 0.9, fits Langmuir vs\n"
             "     Langmuir–Freundlich vs Toth, and places rungs at equal Fisher overlap on a\n"
             "     log-spaced skeleton (a point at least every {:g}× equal-log steps, at most\n"
             "     {} extras per interval).\n\n",
             ladderBottomFillFraction, nitrogenBETMaxLogSpacingMultiplier, nitrogenBETMaxGradientExtrasPerInterval);
}

void pinSystemPengRobinsonPressure(System& system, double pressurePa)
{
  const double T = system.temperature;
  system.input_pressure = pressurePa;
  system.pressure = pressurePa / Units::PressureConversionFactor;
  system.input_pressureTensorDiagonal = double3(pressurePa, pressurePa, pressurePa);
  system.pressureTensorDiagonal = double3(system.pressure, system.pressure, system.pressure);
  for (Component& component : system.components)
  {
    component.fugacityCoefficient = std::nullopt;
  }
  system.equationOfState =
      EquationOfState(EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals, T, pressurePa,
                      system.simulationBox, system.heliumVoidFraction, system.components);
}

namespace
{
struct MonotonePCHIP
{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> d;

  double operator()(double query) const
  {
    if (query <= x.front())
    {
      return y.front();
    }
    if (query >= x.back())
    {
      return y.back();
    }
    const auto upper = std::ranges::upper_bound(x, query);
    const std::size_t i = static_cast<std::size_t>(std::distance(x.begin(), upper) - 1);
    const double h = x[i + 1] - x[i];
    const double t = (query - x[i]) / h;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;
    return h00 * y[i] + h10 * h * d[i] + h01 * y[i + 1] + h11 * h * d[i + 1];
  }

  double invert(double target) const
  {
    const double loY = y.front();
    const double hiY = y.back();
    if (target <= loY)
    {
      return x.front();
    }
    if (target >= hiY)
    {
      return x.back();
    }
    double lo = x.front();
    double hi = x.back();
    for (std::size_t iteration = 0; iteration < 80uz; ++iteration)
    {
      const double mid = 0.5 * (lo + hi);
      if ((*this)(mid) < target)
      {
        lo = mid;
      }
      else
      {
        hi = mid;
      }
    }
    return 0.5 * (lo + hi);
  }
};

MonotonePCHIP makeMonotonePCHIP(std::vector<double> x, std::vector<double> y)
{
  const std::size_t n = x.size();
  std::vector<double> d(n, 0.0);
  if (n == 2uz)
  {
    const double slope = (y[1] - y[0]) / (x[1] - x[0]);
    d[0] = slope;
    d[1] = slope;
    return {std::move(x), std::move(y), std::move(d)};
  }

  std::vector<double> h(n - 1uz);
  std::vector<double> delta(n - 1uz);
  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    h[i] = x[i + 1] - x[i];
    delta[i] = (y[i + 1] - y[i]) / h[i];
  }

  for (std::size_t i = 1; i + 1 < n; ++i)
  {
    if (delta[i - 1] * delta[i] <= 0.0)
    {
      d[i] = 0.0;
    }
    else
    {
      const double w1 = 2.0 * h[i] + h[i - 1];
      const double w2 = h[i] + 2.0 * h[i - 1];
      d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
    }
  }

  d[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);
  if (d[0] * delta[0] < 0.0)
  {
    d[0] = 0.0;
  }
  else if (delta[0] * delta[1] < 0.0 && std::abs(d[0]) > std::abs(3.0 * delta[0]))
  {
    d[0] = 3.0 * delta[0];
  }

  d[n - 1] = ((2.0 * h[n - 2] + h[n - 3]) * delta[n - 2] - h[n - 2] * delta[n - 3]) / (h[n - 2] + h[n - 3]);
  if (d[n - 1] * delta[n - 2] < 0.0)
  {
    d[n - 1] = 0.0;
  }
  else if (delta[n - 2] * delta[n - 3] < 0.0 && std::abs(d[n - 1]) > std::abs(3.0 * delta[n - 2]))
  {
    d[n - 1] = 3.0 * delta[n - 2];
  }

  return {std::move(x), std::move(y), std::move(d)};
}
}  // namespace

std::vector<double> placePressuresByEqualOccupancy(std::span<const double> pressures,
                                                   std::span<const double> occupancies)
{
  if (pressures.size() < 2uz || occupancies.size() != pressures.size())
  {
    return {pressures.begin(), pressures.end()};
  }

  std::vector<double> logP(pressures.size());
  std::vector<double> n(occupancies.begin(), occupancies.end());
  for (std::size_t i = 0; i < pressures.size(); ++i)
  {
    if (!(pressures[i] > 0.0))
    {
      return {pressures.begin(), pressures.end()};
    }
    logP[i] = std::log(pressures[i]);
    if (i > 0uz)
    {
      n[i] = std::max(n[i], n[i - 1]);
      if (!(logP[i] > logP[i - 1]))
      {
        return {pressures.begin(), pressures.end()};
      }
    }
  }

  const double nMin = n.front();
  const double nMax = n.back();
  if (!(nMax > nMin + 1.0e-12))
  {
    return {pressures.begin(), pressures.end()};
  }

  const MonotonePCHIP spline = makeMonotonePCHIP(logP, n);
  std::vector<double> placed(pressures.size());
  const std::size_t last = pressures.size() - 1uz;
  placed.front() = pressures.front();
  placed.back() = pressures.back();
  for (std::size_t k = 1; k < last; ++k)
  {
    const double nTarget = nMin + static_cast<double>(k) / static_cast<double>(last) * (nMax - nMin);
    placed[k] = std::exp(spline.invert(nTarget));
  }
  return placed;
}

double langmuirPressureAtCoverage(double coverage, double affinity)
{
  if (!(affinity > 0.0) || !(coverage > 0.0) || !(coverage < 1.0))
  {
    return 0.0;
  }
  return (coverage / (1.0 - coverage)) / affinity;
}

namespace
{
constexpr double sipsCoverageLow = 0.05;
constexpr double sipsCoverageHigh = 0.95;
constexpr double nSatScanFactor = 1.2;
constexpr std::size_t nSatScanSteps = 25;
constexpr double langmuirExponentTolerance = 0.15;
constexpr double hillRelativeTolerance = 0.4;
constexpr double tothLangmuirExponentLow = 0.85;
constexpr double tothLangmuirExponentHigh = 1.15;
constexpr std::size_t tothFisherTableSize = 1024;

std::vector<double> logSpacedPressures(double lo, double hi, std::size_t count)
{
  std::vector<double> pressures(count);
  const double logMin = std::log(lo);
  const double logMax = std::log(hi);
  for (std::size_t i = 0; i < count; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(count - 1);
    pressures[i] = std::exp(logMin + t * (logMax - logMin));
  }
  pressures.front() = lo;
  pressures.back() = hi;
  return pressures;
}

double sipsCoverage(double pressure, double affinity, double exponent)
{
  if (!(pressure > 0.0) || !(affinity > 0.0) || !(exponent > 0.0))
  {
    return 0.0;
  }
  const double x = std::pow(affinity * pressure, exponent);
  if (!std::isfinite(x) || x > 1.0e300)
  {
    return 1.0;
  }
  return x / (1.0 + x);
}

double tothCoverage(double pressure, double affinity, double exponent)
{
  if (!(pressure > 0.0) || !(affinity > 0.0) || !(exponent > 0.0))
  {
    return 0.0;
  }
  const double x = affinity * pressure;
  const double xt = std::pow(x, exponent);
  if (!std::isfinite(xt) || xt > 1.0e300)
  {
    return 1.0;
  }
  return x / std::pow(1.0 + xt, 1.0 / exponent);
}

double modelCoverage(const NitrogenBETPreIsothermFit& fit, double pressure)
{
  switch (fit.model)
  {
    case NitrogenBETIsothermModel::Toth:
      return tothCoverage(pressure, fit.affinity, fit.exponent);
    case NitrogenBETIsothermModel::LangmuirFreundlich:
      return sipsCoverage(pressure, fit.affinity, fit.exponent);
    case NitrogenBETIsothermModel::Langmuir:
      return sipsCoverage(pressure, fit.affinity, 1.0);
  }
  return sipsCoverage(pressure, fit.affinity, 1.0);
}

double sipsPressureFromCoverage(double coverage, double affinity, double exponent)
{
  const double theta = std::min(std::max(coverage, 1.0e-15), 1.0 - 1.0e-15);
  return std::pow(theta / (1.0 - theta), 1.0 / exponent) / affinity;
}

double tothPressureFromCoverage(double coverage, double affinity, double exponent)
{
  const double theta = std::min(std::max(coverage, 1.0e-15), 1.0 - 1.0e-15);
  const double thetaT = std::pow(theta, exponent);
  if (!(1.0 - thetaT > 0.0))
  {
    return 0.0;
  }
  return std::pow(thetaT / (1.0 - thetaT), 1.0 / exponent) / affinity;
}

double modelPressureFromCoverage(const NitrogenBETPreIsothermFit& fit, double coverage)
{
  switch (fit.model)
  {
    case NitrogenBETIsothermModel::Toth:
      return tothPressureFromCoverage(coverage, fit.affinity, fit.exponent);
    case NitrogenBETIsothermModel::LangmuirFreundlich:
      return sipsPressureFromCoverage(coverage, fit.affinity, fit.exponent);
    case NitrogenBETIsothermModel::Langmuir:
      return sipsPressureFromCoverage(coverage, fit.affinity, 1.0);
  }
  return sipsPressureFromCoverage(coverage, fit.affinity, 1.0);
}

double sipsFisherCoordinate(double coverage)
{
  const double theta = std::min(std::max(coverage, 0.0), 1.0);
  return std::asin(std::sqrt(theta));
}

double tothIntegrand(double coverage, double exponent)
{
  const double theta = std::max(coverage, 1.0e-18);
  const double inner = theta * (1.0 - std::pow(theta, exponent));
  if (!(inner > 0.0))
  {
    return 0.0;
  }
  return 1.0 / std::sqrt(inner);
}

struct FisherTable
{
  std::vector<double> coverage{};
  std::vector<double> coordinate{};
};

FisherTable makeTothFisherTable(double thetaLo, double thetaHi, double exponent)
{
  FisherTable table;
  const double lo = std::min(std::max(thetaLo, 1.0e-12), 1.0 - 1.0e-12);
  const double hi = std::min(std::max(thetaHi, lo + 1.0e-12), 1.0 - 1.0e-12);
  table.coverage.resize(tothFisherTableSize);
  table.coordinate.resize(tothFisherTableSize);
  table.coverage.front() = lo;
  table.coordinate.front() = 0.0;
  const double step = (hi - lo) / static_cast<double>(tothFisherTableSize - 1);
  for (std::size_t i = 1; i < tothFisherTableSize; ++i)
  {
    const double previous = lo + static_cast<double>(i - 1) * step;
    const double current = (i + 1 == tothFisherTableSize) ? hi : lo + static_cast<double>(i) * step;
    table.coverage[i] = current;
    table.coordinate[i] =
        table.coordinate[i - 1] + 0.5 * (tothIntegrand(previous, exponent) + tothIntegrand(current, exponent)) *
                                      (current - previous);
  }
  return table;
}

double interpolateMonotone(const std::vector<double>& x, const std::vector<double>& y, double query)
{
  if (query <= x.front())
  {
    return y.front();
  }
  if (query >= x.back())
  {
    return y.back();
  }
  const auto upper = std::ranges::upper_bound(x, query);
  const std::size_t i = static_cast<std::size_t>(std::distance(x.begin(), upper));
  const double t = (query - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] + t * (y[i] - y[i - 1]);
}

struct LinearSipsFit
{
  double exponent{1.0};
  double intercept{0.0};
  double rSquared{0.0};
  double residualSum{0.0};
  double slopeError{0.0};
  std::size_t numberOfPoints{0};
};

LinearSipsFit fitLinearSips(std::span<const SimulatedIsothermPoint> probes, double nSat, double forcedExponent)
{
  LinearSipsFit fit;
  std::vector<double> x;
  std::vector<double> y;
  x.reserve(probes.size());
  y.reserve(probes.size());
  for (const SimulatedIsothermPoint& point : probes)
  {
    if (!(point.pressure > 0.0) || !(point.moleculesPerCell > sipsCoverageLow * nSat) ||
        !(point.moleculesPerCell < sipsCoverageHigh * nSat))
    {
      continue;
    }
    x.push_back(std::log(point.pressure));
    y.push_back(std::log(point.moleculesPerCell / (nSat - point.moleculesPerCell)));
  }
  fit.numberOfPoints = x.size();
  if (fit.numberOfPoints < 2uz)
  {
    return fit;
  }

  double sumX = 0.0;
  double sumY = 0.0;
  double sumXX = 0.0;
  double sumXY = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    sumX += x[i];
    sumY += y[i];
    sumXX += x[i] * x[i];
    sumXY += x[i] * y[i];
  }
  const double n = static_cast<double>(x.size());
  const double det = n * sumXX - sumX * sumX;
  if (forcedExponent > 0.0)
  {
    fit.exponent = forcedExponent;
    fit.intercept = (sumY - fit.exponent * sumX) / n;
  }
  else if (std::abs(det) < 1.0e-18)
  {
    return fit;
  }
  else
  {
    fit.exponent = (n * sumXY - sumX * sumY) / det;
    fit.intercept = (sumY * sumXX - sumX * sumXY) / det;
  }
  if (!(std::isfinite(fit.exponent) && std::isfinite(fit.intercept)))
  {
    fit.numberOfPoints = 0;
    return fit;
  }

  double ssRes = 0.0;
  double ssTot = 0.0;
  const double meanY = sumY / n;
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    const double predicted = fit.exponent * x[i] + fit.intercept;
    ssRes += (y[i] - predicted) * (y[i] - predicted);
    ssTot += (y[i] - meanY) * (y[i] - meanY);
  }
  fit.residualSum = ssRes;
  fit.rSquared = (ssTot > 0.0) ? 1.0 - ssRes / ssTot : 1.0;
  if (forcedExponent <= 0.0 && x.size() > 2uz && det > 0.0)
  {
    const double mse = ssRes / static_cast<double>(x.size() - 2uz);
    const double sxx = sumXX - sumX * sumX / n;
    if (sxx > 0.0)
    {
      fit.slopeError = std::sqrt(mse / sxx);
    }
  }
  return fit;
}

double tothLoadingSumOfSquares(std::span<const SimulatedIsothermPoint> probes, double nSat, double affinity,
                               double exponent)
{
  double ss = 0.0;
  for (const SimulatedIsothermPoint& point : probes)
  {
    if (!(point.pressure > 0.0))
    {
      continue;
    }
    const double predicted = nSat * tothCoverage(point.pressure, affinity, exponent);
    const double residual = predicted - point.moleculesPerCell;
    ss += residual * residual;
  }
  return ss;
}

double fitTothExponent(std::span<const SimulatedIsothermPoint> probes, double nSat, double affinity)
{
  double bestT = 1.0;
  double bestSS = tothLoadingSumOfSquares(probes, nSat, affinity, 1.0);
  for (int i = 4; i <= 40; ++i)
  {
    const double t = 0.05 * static_cast<double>(i);
    const double ss = tothLoadingSumOfSquares(probes, nSat, affinity, t);
    if (ss < bestSS)
    {
      bestSS = ss;
      bestT = t;
    }
  }
  double lo = std::max(0.15, bestT - 0.05);
  double hi = std::min(2.0, bestT + 0.05);
  for (std::size_t refine = 0; refine < 20uz; ++refine)
  {
    const double left = lo + (hi - lo) / 3.0;
    const double right = hi - (hi - lo) / 3.0;
    if (tothLoadingSumOfSquares(probes, nSat, affinity, left) < tothLoadingSumOfSquares(probes, nSat, affinity, right))
    {
      hi = right;
    }
    else
    {
      lo = left;
    }
  }
  return 0.5 * (lo + hi);
}

double tothRSquared(std::span<const SimulatedIsothermPoint> probes, double nSat, double affinity, double exponent)
{
  double mean = 0.0;
  std::size_t count = 0uz;
  for (const SimulatedIsothermPoint& point : probes)
  {
    if (!(point.pressure > 0.0))
    {
      continue;
    }
    mean += point.moleculesPerCell;
    ++count;
  }
  if (count < 2uz)
  {
    return 0.0;
  }
  mean /= static_cast<double>(count);
  double ssTot = 0.0;
  for (const SimulatedIsothermPoint& point : probes)
  {
    if (!(point.pressure > 0.0))
    {
      continue;
    }
    ssTot += (point.moleculesPerCell - mean) * (point.moleculesPerCell - mean);
  }
  const double ssRes = tothLoadingSumOfSquares(probes, nSat, affinity, exponent);
  return (ssTot > 0.0) ? 1.0 - ssRes / ssTot : 1.0;
}

void hillExponents(std::span<const SimulatedIsothermPoint> probes, double nSat, double& low, double& high, bool& flat)
{
  std::vector<SimulatedIsothermPoint> valid;
  for (const SimulatedIsothermPoint& point : probes)
  {
    if (point.pressure > 0.0 && point.moleculesPerCell > sipsCoverageLow * nSat &&
        point.moleculesPerCell < sipsCoverageHigh * nSat)
    {
      valid.push_back(point);
    }
  }
  std::ranges::sort(valid, {}, &SimulatedIsothermPoint::pressure);
  flat = true;
  low = 1.0;
  high = 1.0;
  if (valid.size() < 3uz)
  {
    return;
  }
  auto yOf = [&](const SimulatedIsothermPoint& point)
  { return std::log(point.moleculesPerCell / (nSat - point.moleculesPerCell)); };
  const double m01 =
      (yOf(valid[1]) - yOf(valid[0])) / (std::log(valid[1].pressure) - std::log(valid[0].pressure));
  const double m12 =
      (yOf(valid[2]) - yOf(valid[1])) / (std::log(valid[2].pressure) - std::log(valid[1].pressure));
  low = m01;
  high = m12;
  const double scale = std::max({std::abs(m01), std::abs(m12), 0.5});
  flat = std::abs(m01 - m12) <= hillRelativeTolerance * scale;
}

bool rejectLangmuirExponent(const LinearSipsFit& sips, const LinearSipsFit& langmuir)
{
  if (sips.numberOfPoints < 2uz || !(sips.rSquared >= 0.85) || !(sips.rSquared + 1.0e-6 >= langmuir.rSquared))
  {
    return false;
  }
  const double delta = std::abs(sips.exponent - 1.0);
  if (sips.numberOfPoints == 2uz)
  {
    return delta > langmuirExponentTolerance;
  }
  if (sips.slopeError > 0.0)
  {
    return delta > langmuirExponentTolerance && delta > 2.0 * sips.slopeError;
  }
  return delta > langmuirExponentTolerance;
}
}  // namespace

NitrogenBETPreIsothermFit fitNitrogenBETPreIsotherm(double henryCoefficientPerCell, double saturationLoadingPerCell,
                                                    double gurvichLoadingPerCell,
                                                    std::span<const SimulatedIsothermPoint> probes)
{
  NitrogenBETPreIsothermFit fit;
  fit.saturationLoadingPerCell = saturationLoadingPerCell;
  if (henryCoefficientPerCell > 0.0 && saturationLoadingPerCell > 0.0)
  {
    fit.affinity = henryCoefficientPerCell / saturationLoadingPerCell;
  }
  fit.exponent = 1.0;
  fit.model = NitrogenBETIsothermModel::Langmuir;
  fit.selectionReason = "Langmuir from K_H and n_sat";

  if (!(saturationLoadingPerCell > 0.0) || probes.empty())
  {
    return fit;
  }

  double maxProbe = 0.0;
  for (const SimulatedIsothermPoint& point : probes)
  {
    maxProbe = std::max(maxProbe, point.moleculesPerCell);
  }
  const double nSatMin = std::max(saturationLoadingPerCell, maxProbe * 1.02);
  double nSatMax = nSatMin;
  if (gurvichLoadingPerCell > nSatMin)
  {
    nSatMax = std::min(nSatScanFactor * saturationLoadingPerCell, gurvichLoadingPerCell);
    nSatMax = std::max(nSatMax, nSatMin);
  }

  LinearSipsFit bestSips;
  double bestNSat = saturationLoadingPerCell;
  for (std::size_t step = 0; step <= nSatScanSteps; ++step)
  {
    const double t = static_cast<double>(step) / static_cast<double>(nSatScanSteps);
    const double nSat = nSatMin + t * (nSatMax - nSatMin);
    const LinearSipsFit candidate = fitLinearSips(probes, nSat, 0.0);
    if (candidate.numberOfPoints >= 2uz &&
        (bestSips.numberOfPoints < 2uz || candidate.rSquared > bestSips.rSquared))
    {
      bestSips = candidate;
      bestNSat = nSat;
    }
  }
  if (bestSips.numberOfPoints < 2uz)
  {
    bestSips = fitLinearSips(probes, saturationLoadingPerCell, 0.0);
    bestNSat = saturationLoadingPerCell;
  }

  const LinearSipsFit langmuirLine = fitLinearSips(probes, bestNSat, 1.0);
  fit.langmuirRSquared = langmuirLine.rSquared;
  fit.rSquared = bestSips.rSquared;
  fit.langmuirRejected = rejectLangmuirExponent(bestSips, langmuirLine);

  bool hillFlat = true;
  hillExponents(probes, bestNSat, fit.hillExponentLow, fit.hillExponentHigh, hillFlat);
  fit.hillPlotFlat = hillFlat;

  if (!hillFlat && fit.affinity > 0.0)
  {
    const double t = fitTothExponent(probes, saturationLoadingPerCell, fit.affinity);
    fit.model = NitrogenBETIsothermModel::Toth;
    fit.saturationLoadingPerCell = saturationLoadingPerCell;
    fit.exponent = t;
    fit.rSquared = tothRSquared(probes, saturationLoadingPerCell, fit.affinity, t);
    if (t >= tothLangmuirExponentLow && t <= tothLangmuirExponentHigh && !fit.langmuirRejected)
    {
      fit.model = NitrogenBETIsothermModel::Langmuir;
      fit.exponent = 1.0;
      fit.rSquared = fit.langmuirRSquared;
      fit.selectionReason = "Toth exponent consistent with Langmuir";
    }
    else
    {
      fit.selectionReason = "Hill plot not flat; Toth with b = K_H / n_sat";
    }
    return fit;
  }

  if (fit.langmuirRejected && bestSips.exponent > 0.05 && bestSips.exponent < 4.0)
  {
    fit.model = NitrogenBETIsothermModel::LangmuirFreundlich;
    fit.saturationLoadingPerCell = bestNSat;
    fit.exponent = bestSips.exponent;
    fit.affinity = std::exp(bestSips.intercept / bestSips.exponent);
    fit.selectionReason = "Linearized Langmuir–Freundlich rejected m = 1";
    return fit;
  }

  return fit;
}

std::vector<double> placePressuresByPreIsothermFisher(double lowestPressure, double highestPressure,
                                                      std::size_t numberOfRungs, const NitrogenBETPreIsothermFit& fit)
{
  if (numberOfRungs < 2uz || !(lowestPressure > 0.0) || !(highestPressure > lowestPressure))
  {
    if (numberOfRungs == 0uz)
    {
      return {};
    }
    if (numberOfRungs == 1uz)
    {
      return {lowestPressure};
    }
    return logSpacedPressures(std::max(lowestPressure, 1.0e-30), std::max(highestPressure, lowestPressure * 1.000001),
                              numberOfRungs);
  }

  if (!(fit.affinity > 0.0) || !(fit.saturationLoadingPerCell > 0.0) || !(fit.exponent > 0.0))
  {
    return logSpacedPressures(lowestPressure, highestPressure, numberOfRungs);
  }

  const double thetaLo = modelCoverage(fit, lowestPressure);
  const double thetaHi = modelCoverage(fit, highestPressure);
  if (!(thetaHi > thetaLo + 1.0e-12))
  {
    return logSpacedPressures(lowestPressure, highestPressure, numberOfRungs);
  }

  FisherTable tothTable;
  const bool useToth = fit.model == NitrogenBETIsothermModel::Toth;
  if (useToth)
  {
    tothTable = makeTothFisherTable(thetaLo, thetaHi, fit.exponent);
  }

  auto fisherOfCoverage = [&](double theta)
  {
    if (useToth)
    {
      return interpolateMonotone(tothTable.coverage, tothTable.coordinate, theta);
    }
    return sipsFisherCoordinate(theta);
  };
  auto coverageOfFisher = [&](double coordinate)
  {
    if (useToth)
    {
      return interpolateMonotone(tothTable.coordinate, tothTable.coverage, coordinate);
    }
    const double chi = std::min(std::max(coordinate, 0.0), 0.5 * std::numbers::pi);
    return std::sin(chi) * std::sin(chi);
  };
  auto fisherOfPressure = [&](double pressure) { return fisherOfCoverage(modelCoverage(fit, pressure)); };
  auto invertFisher = [&](double coordinate)
  {
    const double theta = coverageOfFisher(coordinate);
    double pressure = modelPressureFromCoverage(fit, theta);
    if (!(pressure > 0.0) || !std::isfinite(pressure))
    {
      pressure = lowestPressure;
    }
    return pressure;
  };

  const std::size_t last = numberOfRungs - 1uz;
  const double logRange = std::log(highestPressure / lowestPressure);
  const double maxLogGap = nitrogenBETMaxLogSpacingMultiplier * logRange / static_cast<double>(last);
  const std::size_t skeletonIntervals =
      std::max(1uz, static_cast<std::size_t>(std::ceil(logRange / maxLogGap - 1.0e-12)));
  const std::size_t nSkeleton = std::min(numberOfRungs, skeletonIntervals + 1uz);
  std::vector<double> skeleton = logSpacedPressures(lowestPressure, highestPressure, nSkeleton);

  const std::size_t nIntervals = nSkeleton - 1uz;
  std::vector<double> fisherJump(nIntervals, 0.0);
  for (std::size_t i = 0; i < nIntervals; ++i)
  {
    fisherJump[i] = fisherOfPressure(skeleton[i + 1]) - fisherOfPressure(skeleton[i]);
  }

  std::vector<std::size_t> extra(nIntervals, 0uz);
  const std::size_t extrasToPlace = numberOfRungs - nSkeleton;
  for (std::size_t placedExtra = 0; placedExtra < extrasToPlace; ++placedExtra)
  {
    std::size_t best = nIntervals;
    double bestScore = -1.0;
    for (std::size_t i = 0; i < nIntervals; ++i)
    {
      if (extra[i] >= nitrogenBETMaxGradientExtrasPerInterval)
      {
        continue;
      }
      const double score = fisherJump[i] / static_cast<double>(extra[i] + 1uz);
      if (score > bestScore)
      {
        bestScore = score;
        best = i;
      }
    }
    if (best == nIntervals)
    {
      break;
    }
    ++extra[best];
  }

  std::vector<double> placed;
  placed.reserve(numberOfRungs);
  placed.push_back(skeleton.front());
  for (std::size_t i = 0; i < nIntervals; ++i)
  {
    const double lo = skeleton[i];
    const double hi = skeleton[i + 1];
    const std::size_t k = extra[i];
    const double fLo = fisherOfPressure(lo);
    const double fHi = fisherOfPressure(hi);
    for (std::size_t j = 1; j <= k; ++j)
    {
      double pressure = 0.0;
      if (fHi > fLo + 1.0e-12)
      {
        const double target = fLo + static_cast<double>(j) / static_cast<double>(k + 1uz) * (fHi - fLo);
        pressure = invertFisher(target);
      }
      else
      {
        const double t = static_cast<double>(j) / static_cast<double>(k + 1uz);
        pressure = std::exp(std::log(lo) + t * std::log(hi / lo));
      }
      pressure = std::min(hi, std::max(lo, pressure));
      if (!(pressure > placed.back()))
      {
        pressure = std::nextafter(placed.back(), hi);
      }
      if (pressure >= hi)
      {
        pressure = std::nextafter(hi, lo);
      }
      placed.push_back(pressure);
    }
    placed.push_back(hi);
  }
  while (placed.size() < numberOfRungs)
  {
    std::size_t gap = 1;
    double largestLogGap = 0.0;
    for (std::size_t i = 1; i < placed.size(); ++i)
    {
      const double logGap = std::log(placed[i] / placed[i - 1]);
      if (logGap > largestLogGap)
      {
        largestLogGap = logGap;
        gap = i;
      }
    }
    placed.insert(placed.begin() + static_cast<std::ptrdiff_t>(gap),
                  std::exp(0.5 * (std::log(placed[gap - 1]) + std::log(placed[gap]))));
  }
  return placed;
}

std::vector<double> placePressuresByLangmuirEqualOccupancy(double lowestPressure, double highestPressure,
                                                          std::size_t numberOfRungs, double henryCoefficientPerCell,
                                                          double saturationLoadingPerCell)
{
  NitrogenBETPreIsothermFit fit;
  fit.model = NitrogenBETIsothermModel::Langmuir;
  fit.saturationLoadingPerCell = saturationLoadingPerCell;
  if (henryCoefficientPerCell > 0.0 && saturationLoadingPerCell > 0.0)
  {
    fit.affinity = henryCoefficientPerCell / saturationLoadingPerCell;
  }
  fit.exponent = 1.0;
  return placePressuresByPreIsothermFisher(lowestPressure, highestPressure, numberOfRungs, fit);
}

double nitrogenBETPreIsothermLoadingPerCell(const NitrogenBETPreIsothermFit& fit, double pressure)
{
  return fit.saturationLoadingPerCell * modelCoverage(fit, pressure);
}

namespace
{
constexpr std::size_t scoutBlockCycles = 1000;
constexpr std::size_t scoutMaximumCycles = 15000;
constexpr std::size_t scoutStableBlocks = 3;
constexpr double scoutPlateauTolerance = 0.01;
constexpr double scoutUnconvergedMargin = 1.5;
constexpr double scoutSigmaMargin = 5.0;
constexpr double scoutMinimumMargin = 1.15;
constexpr double scoutGurvichMargin = 1.2;

void pinNitrogenBETPressure(System& system, double pressurePa)
{
  pinSystemPengRobinsonPressure(system, pressurePa);
}

std::size_t gurvichCapForSystem(const System& system, double liquidVolume)
{
  std::size_t gurvichCap = gurvichOccupancy(system.simulationBox.volume, liquidVolume);
  if (system.heliumVoidFraction > 0.0 && system.heliumVoidFraction < 1.0)
  {
    const double packed =
        scoutGurvichMargin * system.heliumVoidFraction * system.simulationBox.volume / liquidVolume;
    gurvichCap = std::max(1uz, static_cast<std::size_t>(std::ceil(packed)));
  }
  return gurvichCap;
}
}  // namespace

NitrogenBETScoutPoint scoutNitrogenBETOccupancy(System system, double pressurePa)
{
  if (system.components.empty())
  {
    throw std::runtime_error("[ComputeBET]: occupancy scout needs an adsorbate component\n");
  }
  if (!(pressurePa > 0.0))
  {
    throw std::runtime_error("[ComputeBET]: occupancy scout needs a positive pressure\n");
  }

  pinNitrogenBETPressure(system, pressurePa);
  system.tmmc.doTMMC = false;
  system.tmmc.rejectOutOfBound = false;

  const std::size_t seed =
      1uz + static_cast<std::size_t>(std::llround(std::abs(std::log(std::max(pressurePa, 1.0e-30))) * 1000.0)) %
                32768uz;
  RandomNumber rng(seed);
  system.forceField.initializeAutomaticCutOff(system.simulationBox);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  std::ostringstream discard;
  system.createExternalFieldInterpolationGrid(discard, 0);
  system.createFrameworkInterpolationGrids(discard);
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();

  const BETProbeProperties probe = requireBETProbeProperties(system.components.front());

  std::size_t fractionalMoleculeSystem = 0uz;
  NitrogenBETScoutPoint point;
  point.pressure = pressurePa;
  point.gurvichCap = gurvichCapForSystem(system, probe.liquidVolume);

  double previousBlockMean = 0.0;
  std::size_t stableBlocks = 0uz;
  double plateauSum = 0.0;
  double plateauSumOfSquares = 0.0;
  std::size_t plateauCycles = 0uz;

  for (Component& component : system.components)
  {
    component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                           system.containsTheFractionalMolecule);
  }

  std::print("  Scouting occupancy at {:.5e} Pa...\n", pressurePa);
  std::cout << std::flush;

  while (point.cycles < scoutMaximumCycles)
  {
    double blockSum = 0.0;
    double blockSumOfSquares = 0.0;
    std::size_t blockPeak = 0uz;
    for (std::size_t cycle = 0uz; cycle < scoutBlockCycles; ++cycle)
    {
      const std::size_t steps = std::max(system.numberOfMolecules(), 20uz) * system.numerOfAdsorbateComponents();
      for (std::size_t step = 0uz; step < steps; ++step)
      {
        const std::size_t selectedComponent = system.randomComponent(rng);
        MC_Moves::performRandomMoveInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        system.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
            system.lambdaWangLandauIsActive(selectedComponent));
      }
      const double occupancy = static_cast<double>(system.numberOfIntegerMoleculesPerComponent[0]);
      blockSum += occupancy;
      blockSumOfSquares += occupancy * occupancy;
      blockPeak = std::max(blockPeak, system.numberOfIntegerMoleculesPerComponent[0]);
    }
    const double blockMean = blockSum / static_cast<double>(scoutBlockCycles);
    const bool isFirstBlock = point.cycles == 0uz;
    point.cycles += scoutBlockCycles;

    if (!isFirstBlock)
    {
      point.peakOccupancy = std::max(point.peakOccupancy, blockPeak);
      const double growth = blockMean - previousBlockMean;
      if (growth <= scoutPlateauTolerance * std::max(1.0, previousBlockMean))
      {
        ++stableBlocks;
        plateauSum += blockSum;
        plateauSumOfSquares += blockSumOfSquares;
        plateauCycles += scoutBlockCycles;
      }
      else
      {
        stableBlocks = 0uz;
        plateauSum = 0.0;
        plateauSumOfSquares = 0.0;
        plateauCycles = 0uz;
      }
    }
    previousBlockMean = blockMean;
    std::print("    {} cycles, mean loading {:.1f}, peak {}\n", point.cycles, blockMean, point.peakOccupancy);
    std::cout << std::flush;

    if (stableBlocks >= scoutStableBlocks)
    {
      point.converged = true;
      break;
    }
  }

  if (plateauCycles > 0uz)
  {
    const double inverseCycles = 1.0 / static_cast<double>(plateauCycles);
    point.meanOccupancy = plateauSum * inverseCycles;
    point.standardDeviation =
        std::sqrt(std::max(0.0, plateauSumOfSquares * inverseCycles - point.meanOccupancy * point.meanOccupancy));
  }
  else
  {
    point.meanOccupancy = previousBlockMean;
  }
  return point;
}

NitrogenBETFillingCeiling scoutNitrogenBETFillingCeiling(System system)
{
  if (system.components.empty())
  {
    throw std::runtime_error(
        "[ComputeBET]: automatic filling ceiling needs an adsorbate component to scout occupancy at P0\n");
  }

  const BETProbeProperties probe = requireBETProbeProperties(system.components.front());

  std::print("  Scouting occupancy at P0 ({:.0f} Pa) for the filling ceiling N_max...\n",
             probe.saturationPressure);
  std::cout << std::flush;

  const NitrogenBETScoutPoint point = scoutNitrogenBETOccupancy(std::move(system), probe.saturationPressure);
  NitrogenBETFillingCeiling ceiling;
  ceiling.meanOccupancy = point.meanOccupancy;
  ceiling.standardDeviation = point.standardDeviation;
  ceiling.peakOccupancy = point.peakOccupancy;
  ceiling.cycles = point.cycles;
  ceiling.converged = point.converged;
  ceiling.gurvichCap = point.gurvichCap;
  ceiling.saturationPressure = probe.saturationPressure;

  const double occupancyCeiling = std::max(ceiling.meanOccupancy + scoutSigmaMargin * ceiling.standardDeviation,
                                           scoutMinimumMargin * ceiling.meanOccupancy);
  const std::size_t scouted = static_cast<std::size_t>(
      std::ceil(ceiling.converged ? occupancyCeiling : scoutUnconvergedMargin * occupancyCeiling));
  ceiling.maxMacrostate = ceiling.gurvichCap;
  if (scouted > 0uz)
  {
    ceiling.maxMacrostate = std::min(ceiling.gurvichCap, scouted + std::max(3uz, (scouted + 15uz) / 16uz));
  }
  ceiling.maxMacrostate = std::max(1uz, ceiling.maxMacrostate);

  std::print("  Filling ceiling N_max = {} (plateau mean {:.1f} ± {:.1f}, Gurvich cap {})\n", ceiling.maxMacrostate,
             ceiling.meanOccupancy, ceiling.standardDeviation, ceiling.gurvichCap);
  std::cout << std::flush;
  return ceiling;
}

void writeNitrogenBETFillingCeiling(std::ostream& stream, const NitrogenBETFillingCeiling& ceiling)
{
  std::print(stream, "BET filling ceiling (unbiased GCMC scout at P0 = {:.0f} Pa)\n", ceiling.saturationPressure);
  std::print(stream, "    Plateau mean occupancy:                {:.2f} ± {:.2f} [molecules]\n", ceiling.meanOccupancy,
             ceiling.standardDeviation);
  std::print(stream, "    Peak occupancy (not used for N_max):   {}\n", ceiling.peakOccupancy);
  std::print(stream, "    Scout cycles:                          {}{}\n", ceiling.cycles,
             ceiling.converged ? "" : " (did not plateau; 1.5× margin applied)");
  std::print(stream, "    Gurvich cap:                           {}\n", ceiling.gurvichCap);
  std::print(stream, "    Macrostate / cycle-length ceiling:     {} molecules\n\n", ceiling.maxMacrostate);
}
