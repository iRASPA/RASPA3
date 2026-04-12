module;

module mc_moves_cputime;

import std;

import double3;
import stringutils;
import archive;
import json;
import mc_moves_move_types;

MCMoveCpuTime::MCMoveCpuTime()
    : timingMap{std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::Translation [0],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                  {"Framework-Molecule", std::chrono::duration<double>::zero()},
                  {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{
                  //Move::Types::RandomTranslation [1],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                  {"Framework-Molecule", std::chrono::duration<double>::zero()},
                  {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::Rotation [2],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                  {"Framework-Molecule", std::chrono::duration<double>::zero()},
                  {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::RandomRotation [3],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField-Molecule", std::chrono::duration<double>::zero()},
                  {"Framework-Molecule", std::chrono::duration<double>::zero()},
                  {"Molecule-Molecule", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::VolumeChange [4],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::ReinsertionCBMC [5],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::PartialReinsertionCBMC [6],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{
                  //Move::Types::IdentityChangeCBMC [7],
                  {"Total", std::chrono::duration<double>::zero()},
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::Swap [8],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"Insertion-Total", std::chrono::duration<double>::zero()},
                  {"Deletion-Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::SwapCBMC [9],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"Insertion-Total", std::chrono::duration<double>::zero()},
                  {"Deletion-Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()},
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::SwapCFCMC [10],
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
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::SwapCBCFCMC [11],
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
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::GibbsVolume [12],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::GibbsSwapCBMC [13],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()},
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::GibbsSwapCFCMC [14],
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
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::Widom [15],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::WidomCFCMC [16],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField", std::chrono::duration<double>::zero()},
                  {"Molecule", std::chrono::duration<double>::zero()},
                  {"Framework", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::WidomCBCFCMC [17],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"ExternalField", std::chrono::duration<double>::zero()},
                  {"Molecule", std::chrono::duration<double>::zero()},
                  {"Framework", std::chrono::duration<double>::zero()},
                  {"Ewald", std::chrono::duration<double>::zero()},
                  {"NonEwald", std::chrono::duration<double>::zero()},
                  {"Tail", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::ParallelTempering [18],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"Energy", std::chrono::duration<double>::zero()},
                  {"Fugacity", std::chrono::duration<double>::zero()}
                },
                std::map<std::string, std::chrono::duration<double>>{ 
                  //Move::Types::HybridMC [19],
                  {"Total", std::chrono::duration<double>::zero()},
                  {"Integration", std::chrono::duration<double>::zero()},
                }
              }
{
}

void MCMoveCpuTime::clearTimingStatistics()
{
  for (auto& moveTimings : timingMap)
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

  for (const Move::Types& moveType : systemMoves)
  {
    auto& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      std::print(stream, "\n");
      std::print(stream, "{:<32} {:14f} [s]\n", Move::moveNames[std::to_underlying(moveType)], moveTimings.at("Total").count());
      for (auto& [timingName, time] : moveTimings)
      {
        if (timingName != "Total")
        {
          std::print(stream, "    {:<28s} {:14f} [s]\n", timingName, time.count());
        }
      }
    }
  }

  for (const Move::Types& moveType : crossSystemMoves)
  {
    auto& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", Move::moveNames[std::to_underlying(moveType)], moveTimings.at("Total").count());
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
  for (const Move::Types& moveType : componentMoves)
  {
    auto& moveTimings = timingMap.at(std::to_underlying(moveType));
    double total = moveTimings.at("Total").count();
    if (total > 0.0)
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", Move::moveNames[std::to_underlying(moveType)], total);
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
  for(std::size_t i = 0; i != timingMap.size(); ++i)
  {
    double total = timingMap[i].at("Total").count();
    if (total > 0.0)
    {
      std::print(stream, "\n");
      std::print(stream, "{:<31} {:14f} [s]\n", Move::moveNames[i], total);
      for (auto& [timingName, time] : timingMap[i])
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

  for (const Move::Types& moveType : systemMoves)
  {
    auto& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      for (auto& [timingName, time] : moveTimings)
      {
        status[Move::moveNames[std::to_underlying(moveType)]][timingName] = time.count();
      }
    }
  }
  for (const Move::Types& moveType : crossSystemMoves)
  {
    auto& moveTimings = timingMap.at(std::to_underlying(moveType));
    if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
    {
      for (auto& [timingName, time] : moveTimings)
      {
        status[Move::moveNames[std::to_underlying(moveType)]][timingName] = time.count();
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

  for(std::size_t i = 0; i != timingMap.size(); ++i)
  {
    if (i == std::to_underlying(Move::Types::VolumeChange) || 
        i == std::to_underlying(Move::Types::GibbsVolume) || 
        i == std::to_underlying(Move::Types::HybridMC))
    {
      continue;
    }

    double total = timingMap[i].at("Total").count();
    if (total > 0.0)
    {
      status[Move::moveNames[i]]["Total"] = total;
      for (auto& [timingName, time] : timingMap[i])
      {
        if (timingName != "Total")
        {
          status[Move::moveNames[i]][timingName] = time.count();

          // skip subtracting keys "...-Total" for overhead (they are summed qts)
          if (timingName.find("Total") == std::string::npos)
          {
            total -= time.count();
          }
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
  std::ostringstream stream;
  nlohmann::json status = jsonComponentMCMoveCPUTimeStatistics();
  status.merge_patch(jsonSystemMCMoveCPUTimeStatistics());

  auto& moveTimings = timingMap.at(std::to_underlying(Move::Types::ParallelTempering));
  if (moveTimings.at("Total") > std::chrono::duration<double>::zero())
  {
    for (auto& [timingName, time] : moveTimings)
    {
      status[Move::moveNames[std::to_underlying(Move::Types::ParallelTempering)]][timingName] = time.count();
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
