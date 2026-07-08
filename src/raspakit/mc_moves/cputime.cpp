module;

module mc_moves_cputime;

import std;

import double3;
import stringutils;
import archive;
import json;
import mc_moves_move_types;

MCMoveCpuTime::MCMoveCpuTime() {}

void MCMoveCpuTime::clearTimingStatistics()
{
  for (TimingRow& moveTimings : timingMap)
  {
    moveTimings.fill(std::chrono::duration<double>::zero());
  }
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics() const
{
  std::ostringstream stream;

  for (const Move::Types& moveType : systemMoves)
  {
    const TimingRow& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings[std::to_underlying(Move::Timing::Total)] > std::chrono::duration<double>::zero())
    {
      std::print(stream, "\n");
      std::print(stream, "{:<32} {:14f} [s]\n", Move::moveNames[std::to_underlying(moveType)],
                 moveTimings[std::to_underlying(Move::Timing::Total)].count());
      for (const Move::Timing& timing : Move::timingKeys[std::to_underlying(moveType)])
      {
        std::print(stream, "    {:<28s} {:14f} [s]\n", Move::timingNames[std::to_underlying(timing)],
                   moveTimings[std::to_underlying(timing)].count());
      }
    }
  }

  for (const Move::Types& moveType : crossSystemMoves)
  {
    const TimingRow& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings[std::to_underlying(Move::Timing::Total)] > std::chrono::duration<double>::zero())
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", Move::moveNames[std::to_underlying(moveType)],
                 moveTimings[std::to_underlying(Move::Timing::Total)].count());
      for (const Move::Timing& timing : Move::timingKeys[std::to_underlying(moveType)])
      {
        std::print(stream, "    {:<27s} {:14f} [s]\n", Move::timingNames[std::to_underlying(timing)],
                   moveTimings[std::to_underlying(timing)].count());
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
  for (const Move::Types& moveType : componentMoves)
  {
    const TimingRow& moveTimings = timingMap.at(std::to_underlying(moveType));
    double total = moveTimings[std::to_underlying(Move::Timing::Total)].count();
    if (total > 0.0)
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", Move::moveNames[std::to_underlying(moveType)], total);
      for (const Move::Timing& timing : Move::timingKeys[std::to_underlying(moveType)])
      {
        const std::string& timingName = Move::timingNames[std::to_underlying(timing)];
        std::print(stream, "    {:<27s} {:14f} [s]\n", timingName, moveTimings[std::to_underlying(timing)].count());

        // skip subtracting keys "...-Total" for overhead (they are summed qts)
        if (timingName.find("Total") == std::string::npos)
        {
          total -= moveTimings[std::to_underlying(timing)].count();
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
  for (std::size_t i = 0; i != timingMap.size(); ++i)
  {
    double total = timingMap[i][std::to_underlying(Move::Timing::Total)].count();
    if (total > 0.0)
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", Move::moveNames[i], total);
      for (const Move::Timing& timing : Move::timingKeys[i])
      {
        const std::string& timingName = Move::timingNames[std::to_underlying(timing)];
        std::print(stream, "    {:<27s} {:14f} [s]\n", timingName, timingMap[i][std::to_underlying(timing)].count());

        // skip subtracting keys "...-Total" for overhead (they are summed qts)
        if (timingName.find("Total") == std::string::npos)
        {
          total -= timingMap[i][std::to_underlying(timing)].count();
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

  auto writeMove = [&](const Move::Types& moveType)
  {
    const TimingRow& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings[std::to_underlying(Move::Timing::Total)] > std::chrono::duration<double>::zero())
    {
      const std::string& moveName = Move::moveNames[std::to_underlying(moveType)];
      status[moveName]["Total"] = moveTimings[std::to_underlying(Move::Timing::Total)].count();
      for (const Move::Timing& timing : Move::timingKeys[std::to_underlying(moveType)])
      {
        status[moveName][Move::timingNames[std::to_underlying(timing)]] =
            moveTimings[std::to_underlying(timing)].count();
      }
    }
  };

  for (const Move::Types& moveType : systemMoves)
  {
    writeMove(moveType);
  }
  for (const Move::Types& moveType : crossSystemMoves)
  {
    writeMove(moveType);
  }

  status["propertySampling"] = propertySampling.count();
  status["energyPressureSampling"] = energyPressureComputation.count();
  return status;
}

const nlohmann::json MCMoveCpuTime::jsonComponentMCMoveCPUTimeStatistics() const
{
  nlohmann::json status;

  for (std::size_t i = 0; i != timingMap.size(); ++i)
  {
    if (i == std::to_underlying(Move::Types::VolumeChange) ||
        i == std::to_underlying(Move::Types::AnisotropicVolumeChange) ||
        i == std::to_underlying(Move::Types::GibbsVolume) ||
        i == std::to_underlying(Move::Types::HybridMC))
    {
      continue;
    }

    double total = timingMap[i][std::to_underlying(Move::Timing::Total)].count();
    if (total > 0.0)
    {
      status[Move::moveNames[i]]["Total"] = total;
      for (const Move::Timing& timing : Move::timingKeys[i])
      {
        const std::string& timingName = Move::timingNames[std::to_underlying(timing)];
        status[Move::moveNames[i]][timingName] = timingMap[i][std::to_underlying(timing)].count();

        // skip subtracting keys "...-Total" for overhead (they are summed qts)
        if (timingName.find("Total") == std::string::npos)
        {
          total -= timingMap[i][std::to_underlying(timing)].count();
        }
      }
      status[Move::moveNames[i]]["Overhead"] = total;
    }
  }
  return status;
}

const nlohmann::json MCMoveCpuTime::jsonOverallMCMoveCPUTimeStatistics(
    std::chrono::duration<double> totalSimulation) const
{
  nlohmann::json status = jsonComponentMCMoveCPUTimeStatistics();
  status.merge_patch(jsonSystemMCMoveCPUTimeStatistics());

  const Move::Types parallelTempering = Move::Types::ParallelTempering;
  const TimingRow& moveTimings = timingMap.at(std::to_underlying(parallelTempering));
  if (moveTimings[std::to_underlying(Move::Timing::Total)] > std::chrono::duration<double>::zero())
  {
    const std::string& moveName = Move::moveNames[std::to_underlying(parallelTempering)];
    status[moveName]["Total"] = moveTimings[std::to_underlying(Move::Timing::Total)].count();
    for (const Move::Timing& timing : Move::timingKeys[std::to_underlying(parallelTempering)])
    {
      status[moveName][Move::timingNames[std::to_underlying(timing)]] =
          moveTimings[std::to_underlying(timing)].count();
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
