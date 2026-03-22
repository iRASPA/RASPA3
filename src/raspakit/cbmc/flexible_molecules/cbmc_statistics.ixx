module;

export module cbmc_move_statistics;

import std;

import archive;
import move_statistics;

export struct CBMCMoveStatistics
{
  std::uint64_t versionNumber{1};

  MoveStatistics<double> bondLengthChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  MoveStatistics<double> bendAngleChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  MoveStatistics<double> conePositionChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};

  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const CBMCMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, CBMCMoveStatistics& p);
};
