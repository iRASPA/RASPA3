module;

module nitrogen_bet_setup;

import std;

import system;
import framework;
import int3;
import isotherm_bet;
import json;

NitrogenBETAutoSetupResult applyNitrogenBETAutoSetup(System& templateSystem, double temperature,
                                                     std::size_t numberOfThreads,
                                                     std::size_t scoutMaximumCycles, bool planPressures,
                                                     bool scoutFillingCeiling)
{
  NitrogenBETAutoSetupResult result;

  if (planPressures)
  {
    result.pressurePlan = planNitrogenBETPressures(templateSystem, temperature, numberOfThreads);
  }

  if (scoutFillingCeiling)
  {
    result.fillingCeiling = scoutNitrogenBETFillingCeiling(templateSystem, scoutMaximumCycles);
    templateSystem.tmmc.maxMacrostate =
        std::max(templateSystem.tmmc.minMacrostate + 1uz, result.fillingCeiling->maxMacrostate);
  }

  return result;
}

void placeNitrogenBETPressureLadderFromHenryAndSaturation(
    std::span<System> systems, std::span<const double> temperatures, std::vector<double>& pressures,
    const NitrogenBETPressurePlan& plan, const std::optional<NitrogenBETFillingCeiling>& fillingCeiling,
    std::size_t scoutMaximumCycles, std::ostream& stream, nlohmann::json& outputJson)
{
  if (systems.empty() || pressures.size() < 3uz || temperatures.empty())
  {
    return;
  }

  const std::size_t numberOfPressures = pressures.size();
  const std::size_t numberOfTemperatures = temperatures.size();

  double cells = 1.0;
  if (systems.front().framework.has_value())
  {
    const int3 numberOfUnitCells = systems.front().framework->numberOfUnitCells;
    cells = static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
  }

  const double henry = plan.henryCoefficientPerCell;
  double nSatPerCell = plan.packingCapacity;
  const char* nSatSource = "Gurvich packing";
  if (fillingCeiling.has_value() && fillingCeiling->meanOccupancy > 0.0 && cells > 0.0)
  {
    nSatPerCell = fillingCeiling->meanOccupancy / cells;
    nSatSource = "P0-scout plateau";
  }
  const double gurvichPerCell = plan.packingCapacity;

  const std::vector<double> logSpacedPressures = pressures;
  const double affinity = (nSatPerCell > 0.0) ? henry / nSatPerCell : 0.0;

  std::vector<SimulatedIsothermPoint> probes;
  if (affinity > 0.0 && nSatPerCell > 0.0)
  {
    std::print("  Scouting Langmuir-coverage probes (θ = 0.1, 0.5, 0.9) for the WHAM pressure ladder...\n");
    std::cout << std::flush;
    for (double coverage : nitrogenBETLangmuirProbeCoverages)
    {
      const double probePressure = langmuirPressureAtCoverage(coverage, affinity);
      if (!(probePressure > logSpacedPressures.front() * 1.0001) ||
          !(probePressure < logSpacedPressures.back() * 0.9999))
      {
        std::print("    skip θ = {:g}: Langmuir P = {:.5e} Pa is outside the Henry–P0 span\n", coverage,
                   probePressure);
        std::cout << std::flush;
        continue;
      }
      const NitrogenBETScoutPoint scout =
          scoutNitrogenBETOccupancy(systems.front(), probePressure, scoutMaximumCycles);
      SimulatedIsothermPoint point;
      point.pressure = probePressure;
      point.moleculesPerCell = (cells > 0.0) ? scout.meanOccupancy / cells : scout.meanOccupancy;
      probes.push_back(point);
      std::print("    θ_L = {:g}: P = {:.5e} Pa, n = {:.4f} / cell (Langmuir {:.4f})\n", coverage, probePressure,
                 point.moleculesPerCell, nSatPerCell * coverage);
      std::cout << std::flush;
    }
  }

  const NitrogenBETPreIsothermFit fit =
      fitNitrogenBETPreIsotherm(henry, nSatPerCell, gurvichPerCell, probes);
  pressures = placePressuresByPreIsothermFisher(logSpacedPressures.front(), logSpacedPressures.back(),
                                                numberOfPressures, fit);

  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      const std::size_t replicaId = temperatureIndex * numberOfPressures + pressureIndex;
      pinSystemPengRobinsonPressure(systems[replicaId], pressures[pressureIndex]);
    }
  }

  const char* modelName = "Langmuir";
  switch (fit.model)
  {
    case NitrogenBETIsothermModel::LangmuirFreundlich:
      modelName = "Langmuir-Freundlich";
      break;
    case NitrogenBETIsothermModel::Toth:
      modelName = "Toth";
      break;
    case NitrogenBETIsothermModel::Langmuir:
      modelName = "Langmuir";
      break;
  }

  std::print(stream, "Pressure-ladder re-placement (pre-isotherm Type I fit, Fisher overlap)\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    Henry coefficient:                     {:.6e} [molecules/cell/Pa]\n", henry);
  std::print(stream, "    n_sat:                                 {:.4f} [molecules / cell] ({})\n", nSatPerCell,
             nSatSource);
  std::print(stream, "    Fitted n_sat:                          {:.4f} [molecules / cell]\n",
             fit.saturationLoadingPerCell);
  std::print(stream, "    Model:                                 {} (b = {:.6e} 1/Pa, exponent = {:.4f})\n", modelName,
             fit.affinity, fit.exponent);
  std::print(stream, "    Selection:                             {}\n", fit.selectionReason);
  std::print(stream, "    Langmuir rejected (m ≠ 1):             {}\n", fit.langmuirRejected ? "yes" : "no");
  std::print(stream, "    Hill plot flat:                        {} (m = {:.3f}, {:.3f})\n",
             fit.hillPlotFlat ? "yes" : "no", fit.hillExponentLow, fit.hillExponentHigh);
  std::print(stream, "    Fit r²:                                {:.5f} (Langmuir line {:.5f})\n", fit.rSquared,
             fit.langmuirRSquared);
  if (affinity > 0.0)
  {
    std::print(stream, "    Langmuir 1/b (half filling):           {:.5e} [Pa]\n", 1.0 / affinity);
  }
  std::print(stream, "    Log-skeleton max gap:                  {:g} × equal-log step\n",
             nitrogenBETMaxLogSpacingMultiplier);
  std::print(stream, "    Fisher extras per skeleton interval:   at most {}\n",
             nitrogenBETMaxGradientExtrasPerInterval);
  if (!probes.empty())
  {
    std::print(stream, "    Langmuir-coverage probes:\n");
    std::print(stream, "      θ_L          P [Pa]          n_scout / cell    n_model / cell\n");
    for (const SimulatedIsothermPoint& point : probes)
    {
      const double langmuirTheta = (affinity * point.pressure) / (1.0 + affinity * point.pressure);
      std::print(stream, "      {:<8.2f}  {:13.5e}    {:14.4f}    {:14.4f}\n", langmuirTheta, point.pressure,
                 point.moleculesPerCell, nitrogenBETPreIsothermLoadingPerCell(fit, point.pressure));
    }
  }
  std::print(stream, "    replica    log-spaced P [Pa]   model n         production P [Pa]   model n\n");
  std::print(stream, "    --------------------------------------------------------------------------------\n");
  for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
  {
    std::print(stream, "    {:7d}    {:13.5e}    {:14.4f}    {:13.5e}    {:14.4f}\n", pressureIndex,
               logSpacedPressures[pressureIndex],
               nitrogenBETPreIsothermLoadingPerCell(fit, logSpacedPressures[pressureIndex]),
               pressures[pressureIndex], nitrogenBETPreIsothermLoadingPerCell(fit, pressures[pressureIndex]));
  }
  std::print(stream, "\n    Production pressure ladder:                ");
  for (double P : pressures)
  {
    std::print(stream, " {}", P);
  }
  std::print(stream, " [Pa]\n");
  std::print(stream, "    Replica grid after re-placement:\n");
  std::print(stream, "    replica    temperature [K]    pressure [Pa]\n");
  std::print(stream, "    ------------------------------------------------\n");
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      const std::size_t replicaId = temperatureIndex * numberOfPressures + pressureIndex;
      std::print(stream, "    {:7d}    {:15.4f}    {:13.5e}\n", replicaId, temperatures[temperatureIndex],
                 pressures[pressureIndex]);
    }
  }
  std::print(stream, "\n");
  std::flush(stream);

  std::print("  Re-placed {} WHAM rungs ({} Fisher ladder, log skeleton at {:g}× equal-log, n_sat = {})\n",
             numberOfPressures, modelName, nitrogenBETMaxLogSpacingMultiplier, nSatSource);
  std::cout << std::flush;

  outputJson["initialization"]["langmuirSaturationLoadingPerCell"] = nSatPerCell;
  outputJson["initialization"]["langmuirSaturationSource"] = nSatSource;
  outputJson["initialization"]["preIsothermModel"] = modelName;
  outputJson["initialization"]["preIsothermAffinity"] = fit.affinity;
  outputJson["initialization"]["preIsothermExponent"] = fit.exponent;
  outputJson["initialization"]["preIsothermSaturationLoadingPerCell"] = fit.saturationLoadingPerCell;
  outputJson["initialization"]["preIsothermLangmuirRejected"] = fit.langmuirRejected;
  outputJson["initialization"]["preIsothermHillPlotFlat"] = fit.hillPlotFlat;
  outputJson["initialization"]["preIsothermRSquared"] = fit.rSquared;

  std::filesystem::create_directories("wham");
  std::ofstream ladderFile("wham/isotherm_ladder.reweighted_histogram.txt", std::ios::trunc);
  std::print(ladderFile, "# WHAM pressure ladder: pre-isotherm {} fit, Fisher overlap on a log-spaced\n", modelName);
  std::print(ladderFile, "# skeleton (gap ≤ {:g}× equal-log) plus at most {} extras per interval.\n",
             nitrogenBETMaxLogSpacingMultiplier, nitrogenBETMaxGradientExtrasPerInterval);
  std::print(ladderFile, "# n_sat source: {}\n", nSatSource);
  std::print(ladderFile, "# selection: {}\n", fit.selectionReason);
  std::print(ladderFile, "# b = {:.8e} 1/Pa, exponent = {:.6f}, fitted n_sat = {:.6f} / cell\n", fit.affinity,
             fit.exponent, fit.saturationLoadingPerCell);
  std::print(ladderFile, "# column 1: log-spaced Henry-to-P0 pressure [Pa]\n");
  std::print(ladderFile, "# column 2: model n at that pressure [molecules / cell]\n");
  std::print(ladderFile, "# column 3: production pressure [Pa]\n");
  std::print(ladderFile, "# column 4: model n at the production pressure [molecules / cell]\n\n");
  for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
  {
    std::print(ladderFile, "{: .6e}   {: .6e}   {: .6e}   {: .6e}\n", logSpacedPressures[pressureIndex],
               nitrogenBETPreIsothermLoadingPerCell(fit, logSpacedPressures[pressureIndex]),
               pressures[pressureIndex], nitrogenBETPreIsothermLoadingPerCell(fit, pressures[pressureIndex]));
  }
}
