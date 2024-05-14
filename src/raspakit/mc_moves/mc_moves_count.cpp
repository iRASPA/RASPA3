module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <string>
#include <sstream>
#include <exception>
#include <source_location>
#include <fstream>
#include <complex>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module mc_moves_count;

#ifndef USE_LEGACY_HEADERS
import <chrono>;
import <string>;
import <sstream>;
import <exception>;
import <source_location>;
import <fstream>;
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

import double3;
import stringutils;
import archive;


void MCMoveCount::clearCountStatistics()
{
  translationMove = size_t{0};
  randomTranslationMove = size_t{0};
  rotationMove = size_t{0};
  randomRotationMove = size_t{0};
  reinsertionMoveCBMC = size_t{0};
  swapInsertionMove = size_t{0};
  swapDeletionMove = size_t{0};
  swapInsertionMoveCBMC = size_t{0};
  swapDeletionMoveCBMC = size_t{0};
  swapLambdaMoveCFCMC = size_t{0};
  swapLambdaMoveCBCFCMC = size_t{0};
  GibbsSwapMoveCBMC = size_t{0};
  GibbsSwapLambdaMoveCFCMC = size_t{0};
  WidomMoveCBMC = size_t{0};
  WidomMoveCFCMC = size_t{0};
  WidomMoveCBCFCMC = size_t{0};
  volumeMove = size_t{0};
  GibbsVolumeMove = size_t{0};
  ParallelTemperingSwap = size_t{0};
}

const std::string MCMoveCount::writeSystemStatistics(size_t countTotal) const
{
  std::ostringstream stream;

  if(volumeMove || GibbsVolumeMove || ParallelTemperingSwap)
  {
    std::print(stream, "\n");
    if(volumeMove)
    {
      std::print(stream, "Volume:                          {:14f} [-]\n", 
                 100.0 * static_cast<double>(volumeMove) / static_cast<double>(countTotal));
    }
    if(GibbsVolumeMove)
    {
      std::print(stream, "Gibbs Volume:                    {:14f} [-]\n", 
                 100.0 * static_cast<double>(GibbsVolumeMove) / static_cast<double>(countTotal));
    }
    if(ParallelTemperingSwap)
    {
      std::print(stream, "Parallel Tempering Swap:         {:14f} [-]\n",
                 100.0 * static_cast<double>(ParallelTemperingSwap) / static_cast<double>(countTotal));
    }

    std::print(stream, "\n\n");
  }

  return stream.str();
}

const std::string MCMoveCount::writeComponentStatistics(size_t countTotal, size_t componentId, 
                                                        const std::string &componentName) const
{
  std::ostringstream stream;

  if(translationMove || randomTranslationMove || rotationMove || randomRotationMove ||
     reinsertionMoveCBMC || swapInsertionMove || swapDeletionMove || swapInsertionMoveCBMC || swapDeletionMoveCBMC ||
     swapLambdaMoveCFCMC || swapLambdaMoveCBCFCMC || GibbsSwapLambdaMoveCFCMC ||
     WidomMoveCBMC || WidomMoveCFCMC || WidomMoveCBCFCMC || 
     volumeMove || GibbsVolumeMove || ParallelTemperingSwap)
  {
    std::print(stream, "Component {} {}\n", componentId, componentName);

    std::print(stream, "\n");
    if(translationMove)
    {
      std::print(stream, "    Translation:                 {:14f} [%]\n", 
                 100.0 * static_cast<double>(translationMove) / static_cast<double>(countTotal));
    }
    if(randomTranslationMove)
    {
      std::print(stream, "    Random Translation:          {:14f} [%]\n", 
                 100.0 * static_cast<double>(randomTranslationMove) / static_cast<double>(countTotal));
    }
    if(rotationMove)
    {
      std::print(stream, "    Rotation:                    {:14f} [%]\n", 
                 100.0 * static_cast<double>(rotationMove) / static_cast<double>(countTotal));
    }
    if(randomRotationMove)
    {
      std::print(stream, "    Random Rotation:             {:14f} [%]\n", 
                 100.0 * static_cast<double>(randomRotationMove) / static_cast<double>(countTotal));
    }
    if(reinsertionMoveCBMC)
    {
      std::print(stream, "    Reinsertion CBMC:            {:14f} [%]\n", 
                 100.0 * static_cast<double>(reinsertionMoveCBMC) / static_cast<double>(countTotal));
    }
    if(swapInsertionMove)
    {
      std::print(stream, "    Insertion:                   {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapInsertionMove) / static_cast<double>(countTotal));
    }
    if(swapDeletionMove)
    {
      std::print(stream, "    Deletion:                    {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapDeletionMove) / static_cast<double>(countTotal));
    }
    if(swapInsertionMoveCBMC)
    {
      std::print(stream, "    Insertion CBMC:              {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapInsertionMoveCBMC) / static_cast<double>(countTotal));
    }
    if(swapDeletionMoveCBMC)
    {
      std::print(stream, "    Deletion CBMC:               {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapDeletionMoveCBMC) / static_cast<double>(countTotal));
    }
    if(swapLambdaMoveCFCMC)
    {
      std::print(stream, "    Insertion/deletion CFCMC:    {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapLambdaMoveCFCMC) / static_cast<double>(countTotal));
    }
    if(swapLambdaMoveCBCFCMC)
    {
      std::print(stream, "    Insertion/deletion CB/CFCMC: {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapLambdaMoveCBCFCMC) / static_cast<double>(countTotal));
    }
    if(GibbsSwapMoveCBMC)
    {
      std::print(stream, "    Gibbs-swap CBMC:             {:14f} [%]\n", 
                 100.0 * static_cast<double>(GibbsSwapMoveCBMC) / static_cast<double>(countTotal));
    }
    if(GibbsSwapLambdaMoveCFCMC)
    {
      std::print(stream, "    Gibbs-swap CFCMC:            {:14f} [%]\n", 
                 100.0 * static_cast<double>(GibbsSwapLambdaMoveCFCMC) / static_cast<double>(countTotal));
    }
    if(WidomMoveCBMC)
    {
      std::print(stream, "    Widom CBMC:                  {:14f} [%]\n", 
                 100.0 * static_cast<double>(WidomMoveCBMC) / static_cast<double>(countTotal));
    }
    if(WidomMoveCFCMC)
    {
      std::print(stream, "    Widom CFCMC:                 {:14f} [%]\n", 
                 100.0 * static_cast<double>(WidomMoveCFCMC) / static_cast<double>(countTotal));
    }
    if(WidomMoveCBCFCMC)
    {
      std::print(stream, "    Widom CB/CFCMC:              {:14f} [%]\n", 
                 100.0 * static_cast<double>(WidomMoveCBCFCMC) / static_cast<double>(countTotal));
    }

    std::print(stream, "\n\n");
  }

  return stream.str();
}

const std::string MCMoveCount::writeAllSystemStatistics(size_t countTotal) const
{
  std::ostringstream stream;

  if(translationMove || randomTranslationMove || rotationMove || randomRotationMove ||
     reinsertionMoveCBMC || swapInsertionMove || swapDeletionMove || swapInsertionMoveCBMC || swapDeletionMoveCBMC ||
     swapLambdaMoveCFCMC || swapLambdaMoveCBCFCMC || GibbsSwapLambdaMoveCFCMC ||
     WidomMoveCBMC || WidomMoveCFCMC || WidomMoveCBCFCMC || 
     volumeMove || GibbsVolumeMove || ParallelTemperingSwap)
  {
    if(translationMove)
    {
      std::print(stream, "Translation:                 {:14f} [%]\n", 
                 100.0 * static_cast<double>(translationMove) / static_cast<double>(countTotal));
    }
    if(randomTranslationMove)
    {
      std::print(stream, "Random Translation:          {:14f} [%]\n", 
                 100.0 * static_cast<double>(randomTranslationMove) / static_cast<double>(countTotal));
    }
    if(rotationMove)
    {
      std::print(stream, "Rotation:                    {:14f} [%]\n", 
                 100.0 * static_cast<double>(rotationMove) / static_cast<double>(countTotal));
    }
    if(randomRotationMove)
    {
      std::print(stream, "Random rotation:             {:14f} [%]\n", 
                 100.0 * static_cast<double>(randomRotationMove) / static_cast<double>(countTotal));
    }
    if(reinsertionMoveCBMC)
    {
      std::print(stream, "Reinsertion CBMC:            {:14f} [%]\n", 
                 100.0 * static_cast<double>(reinsertionMoveCBMC) / static_cast<double>(countTotal));
    }
    if(swapInsertionMove)
    {
      std::print(stream, "Insertion:                   {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapInsertionMove) / static_cast<double>(countTotal));
    }
    if(swapDeletionMove)
    {
      std::print(stream, "Deletion:                    {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapDeletionMove) / static_cast<double>(countTotal));
    }
    if(swapInsertionMoveCBMC)
    {
      std::print(stream, "Insertion CBMC:              {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapInsertionMoveCBMC) / static_cast<double>(countTotal));
    }
    if(swapDeletionMoveCBMC)
    {
      std::print(stream, "Deletion CBMC:               {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapDeletionMoveCBMC) / static_cast<double>(countTotal));
    }
    if(swapLambdaMoveCFCMC)
    {
      std::print(stream, "Insertion/deletion CFCMC:    {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapLambdaMoveCFCMC) / static_cast<double>(countTotal));
    }
    if(swapLambdaMoveCBCFCMC)
    {
      std::print(stream, "Insertion/deletion CB/CFCMC: {:14f} [%]\n", 
                 100.0 * static_cast<double>(swapLambdaMoveCBCFCMC) / static_cast<double>(countTotal));
    }
    if(GibbsSwapMoveCBMC)
    {
      std::print(stream, "Gibbs-swap CBMC:             {:14f} [%]\n", 
                 100.0 * static_cast<double>(GibbsSwapMoveCBMC) / static_cast<double>(countTotal));
    }
    if(GibbsSwapLambdaMoveCFCMC)
    {
      std::print(stream, "Gibbs-swap CFCMC:            {:14f} [%]\n", 
                 100.0 * static_cast<double>(GibbsSwapLambdaMoveCFCMC) / static_cast<double>(countTotal));
    }
    if(volumeMove)
    {
      std::print(stream, "Volume:                      {:14f} [%]\n", 
                 100.0 * static_cast<double>(volumeMove) / static_cast<double>(countTotal));
    }
    if(WidomMoveCBMC)
    {
      std::print(stream, "Widom CBMC:                  {:14f} [%]\n", 
                 100.0 * static_cast<double>(WidomMoveCBMC) / static_cast<double>(countTotal));
    }
    if(WidomMoveCFCMC)
    { 
      std::print(stream, "Widom CFCMC:                 {:14f} [%]\n", 
                 100.0 * static_cast<double>(WidomMoveCFCMC) / static_cast<double>(countTotal));
    }
    if(WidomMoveCBCFCMC)
    {
      std::print(stream, "Widom CB/CFCMC:              {:14f} [%]\n", 
                 100.0 * static_cast<double>(WidomMoveCBCFCMC) / static_cast<double>(countTotal));
    }
    if(GibbsVolumeMove)
    {
      std::print(stream, "Gibbs Volume:                {:14f} [%]\n", 
                 100.0 * static_cast<double>(GibbsVolumeMove) / static_cast<double>(countTotal));
    }
    if(ParallelTemperingSwap)
    {
      std::print(stream, "Parallel Tempering Swap:     {:14f} [%]\n",
                 100.0 * static_cast<double>(ParallelTemperingSwap) / static_cast<double>(countTotal));
    }

    std::print(stream, "\n");
    std::print(stream, "Production mc-move count:    {:14d} [-]\n", countTotal);
    std::print(stream, "                All summed:  {:14d} [-]\n", total());
    std::print(stream, "                difference:  {:14d} [-]\n", countTotal - total());

    std::print(stream, "\n\n");
  }

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveCount &c)
{
  archive << c.versionNumber;

  archive << c.translationMove;
  archive << c.randomTranslationMove;
  archive << c.rotationMove;
  archive << c.randomRotationMove;
  archive << c.reinsertionMoveCBMC;
  archive << c.swapInsertionMove;
  archive << c.swapDeletionMove;
  archive << c.swapInsertionMoveCBMC;
  archive << c.swapDeletionMoveCBMC;
  archive << c.swapLambdaDeletionMove;
  archive << c.swapLambdaMoveCFCMC;
  archive << c.swapLambdaMoveCBCFCMC;
  archive << c.GibbsSwapMoveCBMC;
  archive << c.GibbsSwapLambdaMoveCFCMC;
  archive << c.WidomMoveCBMC;
  archive << c.WidomMoveCFCMC;
  archive << c.WidomMoveCBCFCMC;
  archive << c.volumeMove;
  archive << c.GibbsVolumeMove;
  archive << c.ParallelTemperingSwap;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveCount &c)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > c.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveCpuTime' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> c.translationMove;
  archive >> c.randomTranslationMove;
  archive >> c.rotationMove;
  archive >> c.randomRotationMove;
  archive >> c.reinsertionMoveCBMC;
  archive >> c.swapInsertionMove;
  archive >> c.swapDeletionMove;
  archive >> c.swapInsertionMoveCBMC;
  archive >> c.swapDeletionMoveCBMC;
  archive >> c.swapLambdaDeletionMove;
  archive >> c.swapLambdaMoveCFCMC;
  archive >> c.swapLambdaMoveCBCFCMC;
  archive >> c.GibbsSwapMoveCBMC;
  archive >> c.GibbsSwapLambdaMoveCFCMC;
  archive >> c.WidomMoveCBMC;
  archive >> c.WidomMoveCFCMC;
  archive >> c.WidomMoveCBCFCMC;
  archive >> c.volumeMove;
  archive >> c.GibbsVolumeMove;
  archive >> c.ParallelTemperingSwap;

  return archive;
}
