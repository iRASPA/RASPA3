module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <fstream>
#include <map>
#include <string>
#endif

export module mc_moves_cputime;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import archive;
import json;
import mc_moves_move_types;

/**
 * \brief Stores CPU timing statistics for various Monte Carlo moves.
 *
 * The MCMoveCpuTime struct accumulates the duration of different Monte Carlo moves
 * and their components in a simulation. It provides methods to reset the timings,
 * write the timing statistics to strings, and output them in JSON format.
 */
export struct MCMoveCpuTime
{
  /**
   * \brief Default constructor.
   *
   * Initializes all timing statistics to zero.
   */
  MCMoveCpuTime();

  std::uint64_t versionNumber{2};  ///< Version number for serialization purposes.

  std::map<MoveTypes, std::map<std::string, std::chrono::duration<double>>> timingMap;

  std::chrono::duration<double> propertySampling{0.0};            ///< Time spent on property sampling.
  std::chrono::duration<double> energyPressureComputation{0.0};   ///< Time spent on energy and pressure computations.
  std::chrono::duration<double> pressureFrameworkTime{0.0};       ///< Time spent on energy and pressure computations.
  std::chrono::duration<double> pressureIntermolecularTime{0.0};  ///< Time spent on energy and pressure computations.
  std::chrono::duration<double> pressureEwaldTime{0.0};           ///< Time spent on energy and pressure computations.
  std::chrono::duration<double> pressureTailTime{0.0};            ///< Time spent on energy and pressure computations.
  std::chrono::duration<double> pressureRestTime{0.0};            ///< Time spent on energy and pressure computations.

  /**
   * \brief Calculates the total CPU time spent on all recorded Monte Carlo moves.
   *
   * \return The total duration of all moves.
   */
  inline std::chrono::duration<double> total() const
  {
    std::chrono::duration<double> total{0.0};
    for (auto& [moveType, moveTimings] : timingMap)
    {
      total += moveTimings.at("Total");
    }
    total += propertySampling;
    total += energyPressureComputation;
    return total;
  }

  /**
   * \brief Resets all timing statistics to zero.
   */
  void clearTimingStatistics();

  /**
   * \brief Writes the CPU time statistics to a string.
   *
   * \return A string containing the timing statistics.
   */
  const std::string writeMCMoveCPUTimeStatistics() const;

  /**
   * \brief Writes the CPU time statistics for a specific component.
   *
   * \param componentId The identifier of the component.
   * \param componentName The name of the component.
   * \return A string containing the timing statistics for the component.
   */
  const std::string writeMCMoveCPUTimeStatistics(std::size_t componentId, const std::string& componentName) const;

  /**
   * \brief Writes the overall CPU time statistics.
   *
   * \param total The total simulation time.
   * \return A string containing the overall timing statistics.
   */
  const std::string writeMCMoveCPUTimeStatistics(std::chrono::duration<double> total) const;

  /**
   * \brief Returns the system-level CPU time statistics in JSON format.
   *
   * \return A JSON object containing the system-level timing statistics.
   */
  const nlohmann::json jsonSystemMCMoveCPUTimeStatistics() const;

  /**
   * \brief Returns the component-level CPU time statistics in JSON format.
   *
   * \return A JSON object containing the component-level timing statistics.
   */
  const nlohmann::json jsonComponentMCMoveCPUTimeStatistics() const;

  /**
   * \brief Returns the overall CPU time statistics in JSON format.
   *
   * \param total The total simulation time.
   * \return A JSON object containing the overall timing statistics.
   */
  const nlohmann::json jsonOverallMCMoveCPUTimeStatistics(std::chrono::duration<double> total) const;

  MCMoveCpuTime(const MCMoveCpuTime&) = default;

  std::map<std::string, std::chrono::duration<double>>& operator[](const MoveTypes& move) { return timingMap[move]; }

  inline MCMoveCpuTime& operator=(const MCMoveCpuTime& b)
  {
    propertySampling = b.propertySampling;
    energyPressureComputation = b.energyPressureComputation;
    timingMap = b.timingMap;
    return *this;
  }

  inline MCMoveCpuTime& operator+=(const MCMoveCpuTime& b)
  {
    propertySampling += b.propertySampling;
    energyPressureComputation += b.energyPressureComputation;

    for (auto& [moveType, moveTimings] : timingMap)
    {
      for (auto& [timingName, time] : moveTimings)
      {
        time += b.timingMap.at(moveType).at(timingName);
      }
    }
    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MCMoveCpuTime& t);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MCMoveCpuTime& t);
};

export inline MCMoveCpuTime operator+(const MCMoveCpuTime& a, const MCMoveCpuTime& b)
{
  MCMoveCpuTime m;

  m.propertySampling = a.propertySampling + b.propertySampling;
  m.energyPressureComputation = a.energyPressureComputation + b.energyPressureComputation;

  for (auto& [moveType, moveTimings] : m.timingMap)
  {
    for (auto& [timingName, time] : moveTimings)
    {
      time = a.timingMap.at(moveType).at(timingName) + b.timingMap.at(moveType).at(timingName);
    }
  }

  return m;
}
