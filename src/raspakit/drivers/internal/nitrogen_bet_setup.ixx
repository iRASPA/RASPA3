module;

export module nitrogen_bet_setup;

import std;

import json;
import system;
import isotherm_bet;

/// Shared driver bootstrap for nitrogen BET auto setup (WHAM and parallel TMMC).
/// Henry / packing / scout physics live in 'isotherm_bet'; this module only orchestrates them.

/// Result of optional Henry pressure planning and P0 filling-ceiling scout.
export struct NitrogenBETAutoSetupResult
{
  std::optional<NitrogenBETPressurePlan> pressurePlan;
  std::optional<NitrogenBETFillingCeiling> fillingCeiling;
};

/// Run Henry-based pressure planning and/or a P0 filling-ceiling scout on the template system.
/// When a ceiling is computed, raises 'templateSystem.tmmc.maxMacrostate' to at least that ceiling.
export NitrogenBETAutoSetupResult applyNitrogenBETAutoSetup(System& templateSystem, double temperature,
                                                            std::size_t numberOfThreads,
                                                            std::size_t scoutMaximumCycles, bool planPressures,
                                                            bool scoutFillingCeiling);

/// Re-pin the auto WHAM pressure ladder from a pre-isotherm Type I fit (Langmuir,
/// Langmuir–Freundlich, or Toth): Fisher overlap on a 2× equal-log skeleton.
///
/// 'pressures' enters as the log-spaced Henry–P0 skeleton and leaves as the production ladder.
/// 'systems' is temperature-major: index = temperatureIndex * pressures.size() + pressureIndex.
/// Writes the ladder report to 'stream', stdout, 'outputJson', and wham/isotherm_ladder.*.txt.
export void placeNitrogenBETPressureLadderFromHenryAndSaturation(
    std::span<System> systems, std::span<const double> temperatures, std::vector<double>& pressures,
    const NitrogenBETPressurePlan& plan, const std::optional<NitrogenBETFillingCeiling>& fillingCeiling,
    std::size_t scoutMaximumCycles, std::ostream& stream, nlohmann::json& outputJson);
