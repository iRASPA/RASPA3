module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <numbers>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#endif

module cbmc_move_statistics;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

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

const std::string CBMCMoveStatistics::writeMCMoveStatistics() const
{
  std::ostringstream stream;

  if (bondLengthChange.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Bond-length move", bondLengthChange));
  }

  if (bendAngleChange.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Bend-angle change", bendAngleChange));
  }

  if (conePositionChange.totalCounts > 0.0)
  {
    std::print(stream, "{}", formatStatistics("CBMC Cone-position change", conePositionChange));
  }

  return stream.str();
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const CBMCMoveStatistics& p)
{
  archive << p.versionNumber;

  archive << p.bondLengthChange;
  archive << p.bendAngleChange;
  archive << p.conePositionChange;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, CBMCMoveStatistics& p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'CBMCMoveProbabilitiesSystem' at line {} in file {}\n", location.line(),
                    location.file_name()));
  }

  archive >> p.bondLengthChange;
  archive >> p.bendAngleChange;
  archive >> p.conePositionChange;

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
