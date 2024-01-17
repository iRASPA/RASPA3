module mc_moves_probabilities_system;

import <string>;
import <sstream>;
import <fstream>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

import archive;
import double3;
import move_statistics;
import stringutils;


void MCMoveProbabilitiesSystem::clear()
{
  statistics_VolumeMove.clear();
  statistics_GibbsVolumeMove.clear();
}

void MCMoveProbabilitiesSystem::optimizeAcceptance()
{
  statistics_VolumeMove.optimizeAcceptance(0.01, 1.5);
  statistics_GibbsVolumeMove.optimizeAcceptance(0.01, 1.5);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesSystem &p)
{
  archive << p.versionNumber;

  archive << p.probabilityVolumeMove;
  archive << p.probabilityGibbsVolumeMove;

  archive << p.statistics_VolumeMove;
  archive << p.statistics_GibbsVolumeMove;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesSystem &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.probabilityVolumeMove;
  archive >> p.probabilityGibbsVolumeMove;

  archive >> p.statistics_VolumeMove;
  archive >> p.statistics_GibbsVolumeMove;

  return archive;
}
