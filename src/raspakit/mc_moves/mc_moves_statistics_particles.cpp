module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <sstream>
#include <string>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <map>
#include <source_location>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module mc_moves_statistics_particles;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import archive;
import double3;
import move_statistics;
import stringutils;
import hdf5;

void MCMoveStatisticsParticles::clearMoveStatistics()
{
  translationMove.clear();
  randomTranslationMove.clear();
  rotationMove.clear();
  randomRotationMove.clear();
  reinsertionMove_CBMC.clear();
  identityChangeMove_CBMC.clear();
  swapInsertionMove.clear();
  swapDeletionMove.clear();
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

  swapMove_CFCMC.optimizeAcceptance(0.0, 1.0);
  swapMove_CFCMC_CBMC.optimizeAcceptance(0.0, 1.0);
  GibbsSwapMove_CFCMC.optimizeAcceptance(0.0, 1.0);
}

std::string formatStatistics(const std::string name, const MoveStatistics<double> &move)
{
  std::ostringstream stream;
  std::print(stream, "    {} all:          {:10}\n", name, move.allCounts);
  std::print(stream, "    {} total:        {:10}\n", name, move.totalCounts);
  std::print(stream, "    {} constructed:  {:10}\n", name, move.totalConstructed);
  std::print(stream, "    {} accepted:     {:10}\n", name, move.totalAccepted);
  std::print(stream, "    {} fraction:     {:10f}\n", name,
             move.totalAccepted / std::max(1.0, double(move.totalCounts)));
  std::print(stream, "    {} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

std::string formatStatistics(const std::string name, const MoveStatistics<double3> &move)
{
  std::ostringstream stream;
  std::print(stream, "    {} all:          {:10}\n", name, move.allCounts);
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
  std::print(stream, "    {} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y,
             move.maxChange.z);
  return stream.str();
}

std::vector<double> asVector(const MoveStatistics<double> &move)
{
  std::vector<double> v{move.totalCounts, move.totalConstructed, move.totalAccepted,
                        move.totalAccepted / std::max(1.0, double(move.totalCounts)), move.maxChange};
  return v;
}

std::vector<double> asVector(const MoveStatistics<double3> &move)
{
  std::vector<double> v{move.totalCounts.x,
                        move.totalCounts.y,
                        move.totalCounts.z,
                        move.totalConstructed.x,
                        move.totalConstructed.y,
                        move.totalConstructed.z,
                        move.totalAccepted.x,
                        move.totalAccepted.y,
                        move.totalAccepted.z,
                        move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
                        move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)),
                        move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z)),
                        move.maxChange.x,
                        move.maxChange.y,
                        move.maxChange.z};
  return v;
}

void logStatistics(HDF5Writer &hdf5, const std::string &dataset, const MoveStatistics<double> &move)
{
  std::vector data = asVector(move);
}
void logStatistics(HDF5Writer &hdf5, const std::string &dataset, const MoveStatistics<double3> &move)
{
  std::vector data = asVector(move);
}

const std::string MCMoveStatisticsParticles::writeMCMoveStatistics() const
{
  std::ostringstream stream;
  if (translationMove.totalCounts.x > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Translation", translationMove));
  }
  if (randomTranslationMove.totalCounts.x > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Random translation", randomTranslationMove));
  }
  if (rotationMove.totalCounts.x > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Rotation", rotationMove));
  }
  if (randomRotationMove.totalCounts.x > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Random rotation", randomRotationMove));
  }
  if (reinsertionMove_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Reinsertion(CBMC)", reinsertionMove_CBMC));
  }
  if (identityChangeMove_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Identity Swap (CBMC)", identityChangeMove_CBMC));
  }
  if (swapInsertionMove.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Swap Insertion", swapInsertionMove));
  }
  if (swapDeletionMove.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Swap Deletion", swapDeletionMove));
  }
  if (swapInsertionMove_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Swap Insertion (CBMC)", swapInsertionMove_CBMC));
  }
  if (swapDeletionMove_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Swap Deletion (CBMC)", swapDeletionMove_CBMC));
  }
  if (swapMove_CFCMC.totalCounts.x > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Swap (CFCMC)", swapMove_CFCMC));
  }
  if (swapMove_CFCMC_CBMC.totalCounts.x > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Swap (CB/CFCMC)", swapMove_CFCMC_CBMC));
  }
  if (WidomMove_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Widom (CBMC)", WidomMove_CBMC));
  }
  if (WidomMove_CFCMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Widom (CFCMC)", WidomMove_CFCMC));
  }
  if (WidomMove_CFCMC_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Widom (CB/CFCMC)", WidomMove_CFCMC_CBMC));
  }

  if (GibbsSwapMove_CBMC.totalCounts > 0.0) 
  {
    std::print(stream, "{}", formatStatistics("Gibbs Swap (CBMC)", GibbsSwapMove_CBMC));
  }
  if (GibbsSwapMove_CFCMC.totalCounts.x > 0.0)
  {
    std::print(stream, "{}", formatStatistics("Gibbs Swap (CFCMC)", GibbsSwapMove_CFCMC));
  }

  return stream.str();
}

const std::string MCMoveStatisticsParticles::writeMCMoveStatistics(size_t countTotal, size_t componentId,
                                                                   const std::string &componentName) const
{
  std::ostringstream stream;

  if(translationMove.allCounts || randomTranslationMove.allCounts || 
     rotationMove.allCounts || randomRotationMove.allCounts ||
     reinsertionMove_CBMC.allCounts || 
     swapInsertionMove.allCounts || swapDeletionMove.allCounts ||
     swapInsertionMove_CBMC.allCounts || swapDeletionMove_CBMC.allCounts ||
     swapMove_CFCMC.allCounts || swapMove_CFCMC_CBMC.allCounts || 
     GibbsSwapMove_CBMC.allCounts || GibbsSwapMove_CFCMC.allCounts || GibbsSwapMove_CFCMC_CBMC.allCounts ||
     WidomMove_CBMC.allCounts || WidomMove_CFCMC.allCounts || WidomMove_CFCMC_CBMC.allCounts)
  {
    std::print(stream, "Component {} {}\n", componentId, componentName);
    std::print(stream, "\n");

    if (translationMove.allCounts)
    {
      std::print(stream, "    Translation:                 {:14f} [%]\n", 
                 100.0 * static_cast<double>(translationMove.allCounts) / static_cast<double>(countTotal));
    }
    if(randomTranslationMove.allCounts)
    {
      std::print(stream, "    Random Translation:          {:14f} [%]\n",
                 100.0 * static_cast<double>(randomTranslationMove.allCounts) / static_cast<double>(countTotal));
    }
    if(rotationMove.allCounts)
    {
      std::print(stream, "    Rotation:                    {:14f} [%]\n",
                 100.0 * static_cast<double>(rotationMove.allCounts) / static_cast<double>(countTotal));
    }
    if(randomRotationMove.allCounts)
    {
      std::print(stream, "    Random Rotation:             {:14f} [%]\n",
                 100.0 * static_cast<double>(randomRotationMove.allCounts) / static_cast<double>(countTotal));
    }
    if(reinsertionMove_CBMC.allCounts)
    {
      std::print(stream, "    Reinsertion CBMC:            {:14f} [%]\n",
                 100.0 * static_cast<double>(reinsertionMove_CBMC.allCounts) / static_cast<double>(countTotal));
    }
    if(swapInsertionMove.allCounts)
    {
      std::print(stream, "    Insertion:                   {:14f} [%]\n",
                 100.0 * static_cast<double>(swapInsertionMove.allCounts) / static_cast<double>(countTotal));
    }
    if(swapDeletionMove.allCounts)
    {
      std::print(stream, "    Deletion:                    {:14f} [%]\n",
                 100.0 * static_cast<double>(swapDeletionMove.allCounts) / static_cast<double>(countTotal));
    }
    if(swapInsertionMove_CBMC.allCounts)
    {
      std::print(stream, "    Insertion CBMC:              {:14f} [%]\n",
                 100.0 * static_cast<double>(swapInsertionMove_CBMC.allCounts) / static_cast<double>(countTotal));
    }
    if(swapDeletionMove_CBMC.allCounts)
    {
      std::print(stream, "    Deletion CBMC:               {:14f} [%]\n",
                 100.0 * static_cast<double>(swapDeletionMove_CBMC.allCounts) / static_cast<double>(countTotal));
    }
    if(swapMove_CFCMC.allCounts)
    {
      std::print(stream, "    Insertion/deletion CFCMC:    {:14f} [%]\n",
                 100.0 * static_cast<double>(swapMove_CFCMC.allCounts) / static_cast<double>(countTotal));
    }
    if(swapMove_CFCMC_CBMC.allCounts)
    {
      std::print(stream, "    Insertion/deletion CB/CFCMC: {:14f} [%]\n",
                 100.0 * static_cast<double>(swapMove_CFCMC_CBMC.allCounts) / static_cast<double>(countTotal));
    }
    if(GibbsSwapMove_CBMC.allCounts)
    {
      std::print(stream, "    Gibbs-swap CBMC:             {:14f} [%]\n",
                 100.0 * static_cast<double>(GibbsSwapMove_CBMC.allCounts) / static_cast<double>(countTotal));
    }
    if(GibbsSwapMove_CFCMC.allCounts)
    {
      std::print(stream, "    Gibbs-swap CFCMC:            {:14f} [%]\n",
                 100.0 * static_cast<double>(GibbsSwapMove_CFCMC.allCounts) / static_cast<double>(countTotal));
    }
    if(WidomMove_CBMC.allCounts)
    {
      std::print(stream, "    Widom CBMC:                  {:14f} [%]\n",
                 100.0 * static_cast<double>(WidomMove_CBMC.allCounts) / static_cast<double>(countTotal));
    }
    if(WidomMove_CFCMC.allCounts)
    {
      std::print(stream, "    Widom CFCMC:                 {:14f} [%]\n",
                 100.0 * static_cast<double>(WidomMove_CFCMC.allCounts) / static_cast<double>(countTotal));
    }
    if(WidomMove_CFCMC_CBMC.allCounts)
    {
      std::print(stream, "    Widom CB/CFCMC:              {:14f} [%]\n",
                 100.0 * static_cast<double>(WidomMove_CFCMC_CBMC.allCounts) / static_cast<double>(countTotal));
    }



    std::print(stream, "\n\n");
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
  archive << p.swapInsertionMove;
  archive << p.swapDeletionMove;
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
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'MCMoveStatisticsParticles' "
                    "at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> p.translationMove;
  archive >> p.randomTranslationMove;
  archive >> p.rotationMove;
  archive >> p.randomRotationMove;
  archive >> p.reinsertionMove_CBMC;
  archive >> p.identityChangeMove_CBMC;
  archive >> p.swapInsertionMove;
  archive >> p.swapDeletionMove;
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
