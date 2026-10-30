module;

module cbmc_statistics;

import std;

import json;

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

const std::string CBMC::InternalMoveStatistics::writeMCMoveStatistics() const
{
  std::ostringstream stream;

  if (ringDisplacementChange.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Ring-displacement", ringDisplacementChange));
  }

  if (ringRotationChange.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Ring-rotation", ringRotationChange));
  }

  if (ringCrankshaftMove.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Ring-crankshaft", ringCrankshaftMove));
  }

  if (rigidTiltRotationChange.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Rigid-tilt rotation", rigidTiltRotationChange));
  }

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

nlohmann::json CBMC::InternalMoveStatistics::jsonMCMoveStatistics() const
{
  nlohmann::json status;

  if (ringDisplacementChange.totalCounts > 0.0)
  {
    status["CBMC Ring-displacement"] = jsonStatistics(ringDisplacementChange);
  }
  if (ringRotationChange.totalCounts > 0.0)
  {
    status["CBMC Ring-rotation"] = jsonStatistics(ringRotationChange);
  }
  if (ringCrankshaftMove.totalCounts > 0.0)
  {
    status["CBMC Ring-crankshaft"] = jsonStatistics(ringCrankshaftMove);
  }
  if (rigidTiltRotationChange.totalCounts > 0.0)
  {
    status["CBMC Rigid-tilt rotation"] = jsonStatistics(rigidTiltRotationChange);
  }

  return status;
}

const std::string CBMC::RecoilGrowthStatistics::writeStatistics(std::size_t numberOfTrialDirections) const
{
  std::ostringstream stream;
  if (grows <= 0.0) return stream.str();

  std::print(stream, "    {:20} grows:        {:10}\n", "CBMC Recoil growth", grows);
  std::print(stream, "    {:20} completed:    {:10} ({:.4f})\n", "CBMC Recoil growth", completed, completed / grows);
  std::print(stream, "    {:20} dead ends:    {:10} ({:.4f})\n", "CBMC Recoil growth", deadEnds, deadEnds / grows);
  std::print(stream, "    {:20} discarded:    {:10} ({:.4f})\n", "CBMC Recoil growth", discarded, discarded / grows);
  std::print(stream, "    {:20} mean m_i/k:   {:10.4f} (k = {})\n\n", "CBMC Recoil growth",
             availableDirectionsSum / std::max(1.0, growSteps) / static_cast<double>(numberOfTrialDirections),
             numberOfTrialDirections);
  return stream.str();
}

nlohmann::json CBMC::RecoilGrowthStatistics::jsonStatistics(std::size_t numberOfTrialDirections) const
{
  nlohmann::json status;
  if (grows <= 0.0) return status;

  status["grows"] = grows;
  status["completed"] = completed;
  status["completedFraction"] = completed / grows;
  status["deadEnds"] = deadEnds;
  status["deadEndsFraction"] = deadEnds / grows;
  status["discarded"] = discarded;
  status["discardedFraction"] = discarded / grows;
  status["meanAvailableDirectionsFraction"] =
      availableDirectionsSum / std::max(1.0, growSteps) / static_cast<double>(numberOfTrialDirections);
  status["numberOfTrialDirections"] = numberOfTrialDirections;
  return status;
}

Archive<std::ofstream>& CBMC::operator<<(Archive<std::ofstream>& archive, const CBMC::InternalMoveStatistics& p)
{
  archive << p.versionNumber;

  archive << p.ringDisplacementChange;
  archive << p.ringRotationChange;
  archive << p.ringCrankshaftMove;
  archive << p.rigidTiltRotationChange;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& CBMC::operator>>(Archive<std::ifstream>& archive, CBMC::InternalMoveStatistics& p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'CBMC::InternalMoveStatistics' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> p.ringDisplacementChange;
  archive >> p.ringRotationChange;
  archive >> p.ringCrankshaftMove;
  archive >> p.rigidTiltRotationChange;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("CBMC::InternalMoveStatistics: Error in binary restart\n"));
  }
#endif

  return archive;
}
