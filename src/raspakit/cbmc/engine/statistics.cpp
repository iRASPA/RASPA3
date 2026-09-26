module;

module cbmc_statistics;

import std;

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

  if (versionNumber < 3)
  {
    // Versions 1 and 2 stored the step sizes of the former internal flexible-bead Monte-Carlo
    // (bond-length, bend-angle, cone-position); read and discard them.
    MoveStatistics<double> legacy{};
    archive >> legacy;
    archive >> legacy;
    archive >> legacy;
  }
  archive >> p.ringDisplacementChange;
  archive >> p.ringRotationChange;
  archive >> p.ringCrankshaftMove;
  if (versionNumber >= 2)
  {
    archive >> p.rigidTiltRotationChange;
  }

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
