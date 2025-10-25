module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
#endif

module mc_moves_cputime;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import stringutils;
import archive;
import json;
import mc_moves_move_types;

MCMoveCpuTime::MCMoveCpuTime()
    : timingMap{{MoveTypes::Translation,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                     {"Framework-Molecule", std::chrono::duration<double>::zero()},
                     {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::RandomTranslation,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                     {"Framework-Molecule", std::chrono::duration<double>::zero()},
                     {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::Rotation,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                     {"Framework-Molecule", std::chrono::duration<double>::zero()},
                     {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::RandomRotation,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                     {"Framework-Molecule", std::chrono::duration<double>::zero()},
                     {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::ReinsertionCBMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::PartialReinsertionCBMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::Swap,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"Insertion-Total", std::chrono::duration<double>::zero()},
                     {"Deletion-Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Tail", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::SwapCBMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"Insertion-Total", std::chrono::duration<double>::zero()},
                     {"Deletion-Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Tail", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::SwapCFCMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"Insertion-ExternalField", std::chrono::duration<double>::zero()},
                     {"Insertion-Framework", std::chrono::duration<double>::zero()},
                     {"Insertion-Molecule", std::chrono::duration<double>::zero()},
                     {"Insertion-Ewald", std::chrono::duration<double>::zero()},
                     {"Insertion-Tail", std::chrono::duration<double>::zero()},
                     {"Deletion-ExternalField", std::chrono::duration<double>::zero()},
                     {"Deletion-Framework", std::chrono::duration<double>::zero()},
                     {"Deletion-Molecule", std::chrono::duration<double>::zero()},
                     {"Deletion-Ewald", std::chrono::duration<double>::zero()},
                     {"Deletion-Tail", std::chrono::duration<double>::zero()},
                     {"Lambda-ExternalField", std::chrono::duration<double>::zero()},
                     {"Lambda-Framework", std::chrono::duration<double>::zero()},
                     {"Lambda-Molecule", std::chrono::duration<double>::zero()},
                     {"Lambda-Ewald", std::chrono::duration<double>::zero()},
                     {"Lambda-Tail", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::SwapCBCFCMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"Insertion-ExternalField", std::chrono::duration<double>::zero()},
                     {"Insertion-Framework", std::chrono::duration<double>::zero()},
                     {"Insertion-Molecule", std::chrono::duration<double>::zero()},
                     {"Insertion-NonEwald", std::chrono::duration<double>::zero()},
                     {"Insertion-Ewald", std::chrono::duration<double>::zero()},
                     {"Insertion-Tail", std::chrono::duration<double>::zero()},
                     {"Deletion-ExternalField", std::chrono::duration<double>::zero()},
                     {"Deletion-Framework", std::chrono::duration<double>::zero()},
                     {"Deletion-Molecule", std::chrono::duration<double>::zero()},
                     {"Deletion-NonEwald", std::chrono::duration<double>::zero()},
                     {"Deletion-Ewald", std::chrono::duration<double>::zero()},
                     {"Deletion-Tail", std::chrono::duration<double>::zero()},
                     {"Lambda-ExternalField", std::chrono::duration<double>::zero()},
                     {"Lambda-Framework", std::chrono::duration<double>::zero()},
                     {"Lambda-Molecule", std::chrono::duration<double>::zero()},
                     {"Lambda-NonEwald", std::chrono::duration<double>::zero()},
                     {"Lambda-Ewald", std::chrono::duration<double>::zero()},
                     {"Lambda-Tail", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::GibbsSwapCBMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Tail", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::GibbsSwapCFCMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"LambdaInterchange-NonEwald", std::chrono::duration<double>::zero()},
                     {"LambdaInterchange-Ewald", std::chrono::duration<double>::zero()},
                     {"LambdaInterchange-Tail", std::chrono::duration<double>::zero()},
                     {"LambdaChange-NonEwald", std::chrono::duration<double>::zero()},
                     {"LambdaChange-Ewald", std::chrono::duration<double>::zero()},
                     {"LambdaChange-Tail", std::chrono::duration<double>::zero()},
                     {"LambdaShuffle-NonEwald", std::chrono::duration<double>::zero()},
                     {"LambdaShuffle-Ewald", std::chrono::duration<double>::zero()},
                     {"LambdaShuffle-Tail", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::Widom,
                 {{"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}}},
                {MoveTypes::WidomCFCMC,
                 {{"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField", std::chrono::duration<double>::zero()},
                  {"Molecule", std::chrono::duration<double>::zero()},
                  {"Framework", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()}}},
                {MoveTypes::WidomCBCFCMC,
                 {{"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField", std::chrono::duration<double>::zero()},
                  {"Molecule", std::chrono::duration<double>::zero()},
                  {"Framework", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()}}},
                {MoveTypes::VolumeChange,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Tail", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::GibbsVolume,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"NonEwald", std::chrono::duration<double>::zero()},
                     {"Tail", std::chrono::duration<double>::zero()},
                     {"Ewald", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::ParallelTempering,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"Energy", std::chrono::duration<double>::zero()},
                     {"Fugacity", std::chrono::duration<double>::zero()},
                 }},
                {MoveTypes::HybridMC,
                 {
                     {"Total", std::chrono::duration<double>::zero()},
                     {"Integration", std::chrono::duration<double>::zero()},
                 }}}
{
}

void MCMoveCpuTime::clearTimingStatistics()
{
  for (auto& [moveType, moveTimings] : timingMap)
  {
    for (auto& [timingName, time] : moveTimings)
    {
      time = std::chrono::duration<double>::zero();
    }
  }
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics() const
{
  std::ostringstream stream;

  for (const MoveTypes& moveType : systemMoves)
  {
    if (timingMap.find(moveType) == timingMap.end()) continue;
    auto& moveTimings = timingMap.at(moveType);
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      std::print(stream, "\n");
      std::print(stream, "{:<32} {:14f} [s]\n", moveNames[moveType], moveTimings.at("Total").count());
      for (auto& [timingName, time] : moveTimings)
      {
        if (timingName != "Total")
        {
          std::print(stream, "    {:<28s} {:14f} [s]\n", timingName, time.count());
        }
      }
    }
  }

  for (const MoveTypes& moveType : crossSystemMoves)
  {
    if (timingMap.find(moveType) == timingMap.end()) continue;
    auto& moveTimings = timingMap.at(moveType);
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", moveNames[moveType], moveTimings.at("Total").count());
      for (auto& [timingName, time] : moveTimings)
      {
        if (timingName != "Total")
        {
          std::print(stream, "    {:<27s} {:14f} [s]\n", timingName, time.count());
        }
      }
    }
  }

  std::print(stream, "\n");
  std::print(stream, "Property sampling               {:14f} [s]\n", propertySampling.count());
  std::print(stream, "Energy/pressure sampling:       {:14f} [s]\n", energyPressureComputation.count());
  std::print(stream, "\n\n");

  return stream.str();
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics(std::size_t componentId,
                                                              const std::string& componentName) const
{
  std::ostringstream stream;
  std::print(stream, "Component {} {}\n", componentId, componentName);
  for (const MoveTypes& moveType : componentMoves)
  {
    if (timingMap.find(moveType) == timingMap.end()) continue;
    auto& moveTimings = timingMap.at(moveType);
    double total = moveTimings.at("Total").count();
    if (total > 0.0)
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", moveNames[moveType], total);
      for (auto& [timingName, time] : moveTimings)
      {
        if (timingName != "Total")
        {
          std::print(stream, "    {:<27s} {:14f} [s]\n", timingName, time.count());

          // skip subtracting keys "...-Total" for overhead (they are summed qts)
          if (timingName.find("Total") == std::string::npos)
          {
            total -= time.count();
          }
        }
      }
      std::print(stream, "    {:<27s} {:14f} [s]\n", "Overhead", total);
    }
  }
  std::print(stream, "\n\n");
  return stream.str();
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics(std::chrono::duration<double> totalSimulation) const
{
  std::ostringstream stream;
  for (auto& [moveType, moveTimings] : timingMap)
  {
    double total = moveTimings.at("Total").count();
    if (total > 0.0)
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", moveNames[moveType], total);
      for (auto& [timingName, time] : moveTimings)
      {
        if (timingName != "Total")
        {
          std::print(stream, "    {:<27s} {:14f} [s]\n", timingName, time.count());

          // skip subtracting keys "...-Total" for overhead (they are summed qts)
          if (timingName.find("Total") == std::string::npos)
          {
            total -= time.count();
          }
        }
      }
      std::print(stream, "    {:<27s} {:14f} [s]\n", "Overhead", total);
    }
  }

  std::print(stream, "\n");
  std::print(stream, "Property sampling:              {:14f} [s]\n", propertySampling.count());
  std::print(stream, "Energy/pressure sampling:       {:14f} [s]\n\n", energyPressureComputation.count());

  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalSimulation.count());
  std::print(stream, "                All summed:     {:14f} [s]\n", total().count());
  std::print(stream, "                  Overhead:     {:14f} [s]\n", totalSimulation.count() - total().count());
  std::print(stream, "\n\n");

  return stream.str();
}

const nlohmann::json MCMoveCpuTime::jsonSystemMCMoveCPUTimeStatistics() const
{
  nlohmann::json status;

  for (const MoveTypes& moveType : systemMoves)
  {
    if (timingMap.find(moveType) == timingMap.end()) continue;
    auto& moveTimings = timingMap.at(moveType);
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      for (auto& [timingName, time] : moveTimings)
      {
        status[moveNames[moveType]][timingName] = time.count();
      }
    }
  }
  for (const MoveTypes& moveType : crossSystemMoves)
  {
    if (timingMap.find(moveType) == timingMap.end()) continue;
    auto& moveTimings = timingMap.at(moveType);
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      for (auto& [timingName, time] : moveTimings)
      {
        status[moveNames[moveType]][timingName] = time.count();
      }
    }
  }

  status["propertySampling"] = propertySampling.count();
  status["energyPressureSampling"] = energyPressureComputation.count();
  return status;
}

const nlohmann::json MCMoveCpuTime::jsonComponentMCMoveCPUTimeStatistics() const
{
  nlohmann::json status;

  for (auto& [moveType, moveTimings] : timingMap)
  {
    if (moveType == MoveTypes::VolumeChange || moveType == MoveTypes::GibbsVolume || moveType == MoveTypes::HybridMC)
    {
      continue;
    }

    double total = moveTimings.at("Total").count();
    if (total > 0.0)
    {
      status[moveNames[moveType]]["Total"] = total;
      for (auto& [timingName, time] : moveTimings)
      {
        if (timingName != "Total")
        {
          status[moveNames[moveType]][timingName] = time.count();

          // skip subtracting keys "...-Total" for overhead (they are summed qts)
          if (timingName.find("Total") == std::string::npos)
          {
            total -= time.count();
          }
        }
      }
      status[moveNames[moveType]]["Overhead"] = total;
    }
  }
  return status;
}

const nlohmann::json MCMoveCpuTime::jsonOverallMCMoveCPUTimeStatistics(
    std::chrono::duration<double> totalSimulation) const
{
  std::ostringstream stream;
  nlohmann::json status = jsonComponentMCMoveCPUTimeStatistics();
  status.merge_patch(jsonSystemMCMoveCPUTimeStatistics());

  auto& moveTimings = timingMap.at(MoveTypes::ParallelTempering);
  if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
  {
    for (auto& [timingName, time] : moveTimings)
    {
      status[moveNames[MoveTypes::ParallelTempering]][timingName] = time.count();
    }
  }

  status["productionSimulationTime"] = totalSimulation.count();
  status["allSummed"] = total().count();
  status["overhead"] = totalSimulation.count() - total().count();
  return status;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MCMoveCpuTime& t)
{
  archive << t.versionNumber;

  archive << t.propertySampling;
  archive << t.energyPressureComputation;
  archive << t.timingMap;

  archive << t.propertySampling;
  archive << t.energyPressureComputation;
  archive << t.pressureFrameworkTime;
  archive << t.pressureIntermolecularTime;
  archive << t.pressureEwaldTime;
  archive << t.pressureTailTime;
  archive << t.pressureRestTime;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MCMoveCpuTime& t)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > t.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveCpuTime' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> t.propertySampling;
  archive >> t.energyPressureComputation;
  archive >> t.timingMap;

  archive >> t.propertySampling;
  archive >> t.energyPressureComputation;
  archive >> t.pressureFrameworkTime;
  archive >> t.pressureIntermolecularTime;
  archive >> t.pressureEwaldTime;
  archive >> t.pressureTailTime;
  archive >> t.pressureRestTime;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("MCMoveCpuTime: Error in binary restart\n"));
  }
#endif

  return archive;
}
