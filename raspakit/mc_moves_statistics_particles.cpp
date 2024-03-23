module mc_moves_statistics_particles;

import <string>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#else
  import print;
#endif


import archive;
import double3;
import move_statistics;
import stringutils;


void MCMoveStatisticsParticles::clearMoveStatistics()
{
  translationMove.clear();
  randomTranslationMove.clear();
  rotationMove.clear();
  randomRotationMove.clear();
  reinsertionMove_CBMC.clear();
  identityChangeMove_CBMC.clear();
  swapInsertionMove_CBMC.clear();
  swapDeletionMove_CBMC.clear();
  swapMove_CFCMC.clear();
  swapMove_CFCMC_CBMC.clear();
  WidomMove_CBMC.clear();
  WidomMove_CFCMC.clear();
  WidomMove_CFCMC_CBMC.clear();

  GibbsSwapMove_CBMC.clear();
  GibbsSwapMove_CFCMC.clear();
}

void MCMoveStatisticsParticles::optimizeMCMoves()
{
  translationMove.optimizeAcceptance(0.01, 1.5);
  rotationMove.optimizeAcceptance(0.01, 1.5);

  swapMove_CFCMC_CBMC.optimizeAcceptance(0.0, 1.0);
  GibbsSwapMove_CFCMC.optimizeAcceptance(0.0, 1.0);
}

std::string formatStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10}\n", name, move.totalCounts);
  std::print(stream, "    {} constructed:  {:10}\n", name, move.totalConstructed);
  std::print(stream, "    {} accepted:     {:10}\n", name, move.totalAccepted);
  std::print(stream, "    {} fraction:     {:10f}\n", name, move.totalAccepted / std::max(1.0, double(move.totalCounts)));
  std::print(stream, "    {} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

std::string formatStatistics(const std::string name, const MoveStatistics<double3>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10} {:10} {:10}\n", 
                     name, move.totalCounts.x, move.totalCounts.y, move.totalCounts.z);
  std::print(stream, "    {} constructed:  {:10} {:10} {:10}\n", 
                     name, move.totalConstructed.x, move.totalConstructed.y, move.totalConstructed.z);
  std::print(stream, "    {} accepted:     {:10} {:10} {:10}\n", 
                     name, move.totalAccepted.x, move.totalAccepted.y, move.totalAccepted.z);
  std::print(stream, "    {} fraction:     {:10f} {:10f} {:10f}\n", 
                     name, move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
                     move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)), 
                     move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z)));
  std::print(stream, "    {} max-change:   {:10f} {:10f} {:10f}\n\n", 
                     name, move.maxChange.x, move.maxChange.y, move.maxChange.z);
  return stream.str();
}

const std::string MCMoveStatisticsParticles::writeMCMoveStatistics() const
{
  std::ostringstream stream;
  if (translationMove.totalCounts.x  > 0) 
  {
    std::print(stream, "{}", formatStatistics("Translation", translationMove));
  }
  if (randomTranslationMove.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Random translation", randomTranslationMove));
  }
  if (rotationMove.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Rotation", rotationMove));
  }
  if (randomRotationMove.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Random rotation", randomRotationMove));
  }
  if (reinsertionMove_CBMC.totalCounts > 0) 
  {
    std::print(stream, "{}", formatStatistics("Reinsertion(CBMC)", reinsertionMove_CBMC));
  }
  if (identityChangeMove_CBMC.totalCounts > 0) 
  {
    std::print(stream, "{}", formatStatistics("Identity Swap (CBMC)", identityChangeMove_CBMC));
  }
  if (swapInsertionMove_CBMC.totalCounts > 0) 
  {
    std::print(stream, "{}", formatStatistics("Swap Insertion (CBMC)", swapInsertionMove_CBMC));
  }
  if (swapDeletionMove_CBMC.totalCounts > 0) 
  {
    std::print(stream, "{}", formatStatistics("Swap Deletion (CBMC)", swapDeletionMove_CBMC));
  }
  if (swapMove_CFCMC.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Swap (CFCMC)", swapMove_CFCMC));
  }
  if (swapMove_CFCMC_CBMC.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Swap (CB/CFCMC)", swapMove_CFCMC_CBMC));
  }
  if (WidomMove_CBMC.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Widom (CBMC)", WidomMove_CBMC));
  }
  if (WidomMove_CFCMC.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Widom (CFCMC)", WidomMove_CFCMC));
  }
  if (WidomMove_CFCMC_CBMC.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Widom (CB/CFCMC)", WidomMove_CFCMC_CBMC));
  }

  if (GibbsSwapMove_CBMC.totalCounts > 0) 
  {
    std::print(stream, "{}", formatStatistics("Gibbs Swap (CBMC)", GibbsSwapMove_CBMC));
  }
  if (GibbsSwapMove_CFCMC.totalCounts.x > 0) 
  {
    std::print(stream, "{}", formatStatistics("Gibbs Swap (CFCMC)", GibbsSwapMove_CFCMC));
  }

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsParticles &p)
{
  archive << p.versionNumber;

  archive << p.translationMove;
  archive << p.randomTranslationMove;
  archive << p.rotationMove;
  archive << p.randomRotationMove;
  archive << p.reinsertionMove_CBMC;
  archive << p.identityChangeMove_CBMC;
  archive << p.swapInsertionMove_CBMC;
  archive << p.swapDeletionMove_CBMC;
  archive << p.swapMove_CFCMC;
  archive << p.swapMove_CFCMC_CBMC;
  archive << p.WidomMove_CBMC;
  archive << p.WidomMove_CFCMC;
  archive << p.WidomMove_CFCMC_CBMC;

  archive << p.GibbsSwapMove_CBMC;
  archive << p.GibbsSwapMove_CFCMC;
  archive << p.GibbsSwapMove_CFCMC_CBMC;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsParticles &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveStatisticsParticles' "
                                         "at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.translationMove;
  archive >> p.randomTranslationMove;
  archive >> p.rotationMove;
  archive >> p.randomRotationMove;
  archive >> p.reinsertionMove_CBMC;
  archive >> p.identityChangeMove_CBMC;
  archive >> p.swapInsertionMove_CBMC;
  archive >> p.swapDeletionMove_CBMC;
  archive >> p.swapMove_CFCMC;
  archive >> p.swapMove_CFCMC_CBMC;
  archive >> p.WidomMove_CBMC;
  archive >> p.WidomMove_CFCMC;
  archive >> p.WidomMove_CFCMC_CBMC;

  archive >> p.GibbsSwapMove_CBMC;
  archive >> p.GibbsSwapMove_CFCMC;
  archive >> p.GibbsSwapMove_CFCMC_CBMC;

  return archive;
}

