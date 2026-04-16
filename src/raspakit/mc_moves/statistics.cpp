module;

module mc_moves_statistics;

import std;

import archive;
import mc_moves_move_types;

void MCMoveStatistics::clearMoveStatistics()
{
  for(auto& stat: stats)
  {
    std::visit([](auto&& s){ s.clear(); }, stat);
  }
}

void MCMoveStatistics::optimizeMCMoves()
{
  for(auto& stat: stats)
  {
    std::visit([](auto&& s){ s.optimizeAcceptance(); }, stat);
  }
}

void MCMoveStatistics::addAllCounts(const Move::Types& move)
{
  std::visit([](auto&& s){ s.allCounts += 1; }, stats[std::to_underlying(move)]);
}

void MCMoveStatistics::addTrial(const Move::Types& move)
{
  std::visit([](auto&& s){ s.counts += 1; s.totalCounts += 1;}, stats[std::to_underlying(move)]);
}

void MCMoveStatistics::addTrial(const Move::Types& move, std::size_t direction)
{
  if(std::holds_alternative<MoveStatistics<double3>>(stats[std::to_underlying(move)]))
  {
    std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).counts[direction] += 1.0;
    std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).totalCounts[direction] += 1.0;
  }
}

void MCMoveStatistics::addConstructed(const Move::Types& move)
{
  std::visit([](auto&& s){ s.constructed += 1; s.totalConstructed += 1; }, stats[std::to_underlying(move)]);
}

void MCMoveStatistics::addConstructed(const Move::Types& move, std::size_t direction)
{
  if(std::holds_alternative<MoveStatistics<double3>>(stats[std::to_underlying(move)]))
  {
    std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).constructed[direction] += 1.0;
    std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).totalConstructed[direction] += 1.0;
  }
}

void MCMoveStatistics::addAccepted(const Move::Types& move)
{
  std::visit([](auto&& s){ s.accepted += 1; s.totalAccepted += 1; }, stats[std::to_underlying(move)]);
}

void MCMoveStatistics::addAccepted(const Move::Types& move, std::size_t direction)
{
  if(std::holds_alternative<MoveStatistics<double3>>(stats[std::to_underlying(move)]))
  {
    std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).accepted[direction] += 1.0;
    std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).totalAccepted[direction] += 1.0;
  }
}

static std::string formatStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {:20} all:          {:10}\n", name, move.allCounts);
  std::print(stream, "    {:20} total:        {:10}\n", name, move.totalCounts);
  std::print(stream, "    {:20} constructed:  {:10}\n", name, move.totalConstructed);
  std::print(stream, "    {:20} accepted:     {:10}\n", name, move.totalAccepted);
  std::print(stream, "    {:20} fraction:     {:10f}\n", name,
             move.totalAccepted / std::max(1.0, double(move.totalCounts)));
  std::print(stream, "    {:20} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

static std::string formatStatistics(const std::string name, const MoveStatistics<double3>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {:20} all:          {:10}\n", name, move.allCounts);
  std::print(stream, "    {:20} total:        {:10} {:10} {:10}\n", name, move.totalCounts.x, move.totalCounts.y,
             move.totalCounts.z);
  std::print(stream, "    {:20} constructed:  {:10} {:10} {:10}\n", name, move.totalConstructed.x,
             move.totalConstructed.y, move.totalConstructed.z);
  std::print(stream, "    {:20} accepted:     {:10} {:10} {:10}\n", name, move.totalAccepted.x, move.totalAccepted.y,
             move.totalAccepted.z);
  std::print(stream, "    {:20} fraction:     {:10f} {:10f} {:10f}\n", name,
             move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
             move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)),
             move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z)));
  std::print(stream, "    {:20} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y,
             move.maxChange.z);
  return stream.str();
}

static nlohmann::json jsonStatistics(const MoveStatistics<double>& move)
{
  nlohmann::json status;
  status["all"] = move.allCounts;
  status["total"] = move.totalCounts;
  status["constructed"] = move.totalConstructed;
  status["accepted"] = move.totalAccepted;
  status["fraction"] = move.totalAccepted / std::max(1.0, double(move.totalCounts));
  status["maxChange"] = move.maxChange;
  return status;
}

static nlohmann::json jsonStatistics(const MoveStatistics<double3>& move)
{
  nlohmann::json status;
  status["all"] = move.allCounts;
  status["total"] = {move.totalCounts.x, move.totalCounts.y, move.totalCounts.z};
  status["constructed"] = {move.totalConstructed.x, move.totalConstructed.y, move.totalConstructed.z};
  status["accepted"] = {move.totalAccepted.x, move.totalAccepted.y, move.totalAccepted.z};
  status["fraction"] = {move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
                        move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)),
                        move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z))};
  status["maxChange"] = {move.maxChange.x, move.maxChange.y, move.maxChange.z};
  return status;
}

const std::string MCMoveStatistics::writeMCMoveStatistics() const
{
  std::ostringstream stream;
  for(std::size_t i = 0; i != stats.size(); ++i)
  {
    std::visit([&stream, i](auto&& s){
        if(s.allCounts > 0)
        {
          std::print(stream, "{}", formatStatistics(Move::moveNames[i], s)); 
        }
    }, stats[i]);
  }
  return stream.str();
}

const std::string MCMoveStatistics::writeMCMoveStatistics(std::size_t countTotal) const
{
  std::ostringstream stream;

  std::size_t summed = 0;
  for (std::size_t i = 0; i != stats.size(); ++i)
  {
    std::visit([&stream, i, &summed, countTotal](auto&& s){  
        std::size_t moveCount = s.allCounts;
        if (moveCount > 0)
        {
          std::print(stream, "{:<29}{:14} ({:<6.4f} [%])\n", Move::moveNames[i], moveCount,
                     100.0 * static_cast<double>(moveCount) / static_cast<double>(countTotal));
          summed += s.allCounts;
        }
      }, stats[i]);
  }

  std::print(stream, "\n");
  std::print(stream, "Production count MC-steps:   {:14d} [-]\n", countTotal);
  std::print(stream, "               All summed:   {:14d} [-]\n", summed);
  std::print(stream, "               difference:   {:14d} [-]\n", countTotal - summed);

  return stream.str();
}

const nlohmann::json MCMoveStatistics::jsonMCMoveStatistics() const
{
  nlohmann::json status;
  for (std::size_t i = 0; i != stats.size(); ++i)
  {
    std::visit([i, &status](auto&& s){  
       status[Move::moveNames[i]] = jsonStatistics(s);
      }, stats[i]);
  }
  return status;
}


std::string MCMoveStatistics::repr() const
{
  std::ostringstream stream;
  for(std::size_t i = 0; i != stats.size(); ++i)
  {
    std::visit([&stream, i](auto&& s){
        if(s.allCounts > 0)
        {
          std::print(stream, "{}", formatStatistics(Move::moveNames[i], s)); 
        }
    }, stats[i]);
  }
  return stream.str();
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MCMoveStatistics& p)
{
  archive << p.versionNumber;
  archive << p.stats;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MCMoveStatistics& p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.stats;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("MCMoveStatistics: Error in binary restart\n"));
  }
#endif

  return archive;
}
