module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
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
import <print>;
#endif

import double3;
import stringutils;
import archive;
import json;

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

/*
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
*/

const std::string MCMoveCount::writeAllSystemStatistics(size_t countTotal) const
{
  std::ostringstream stream;

  if (translationMove || randomTranslationMove || rotationMove || randomRotationMove || reinsertionMoveCBMC ||
      swapInsertionMove || swapDeletionMove || swapInsertionMoveCBMC || swapDeletionMoveCBMC || swapLambdaMoveCFCMC ||
      swapLambdaMoveCBCFCMC || GibbsSwapLambdaMoveCFCMC || WidomMoveCBMC || WidomMoveCFCMC || WidomMoveCBCFCMC ||
      volumeMove || GibbsVolumeMove || ParallelTemperingSwap)
  {
    if (translationMove)
    {
      std::print(stream, "Translation:                 {:14f} [%]\n",
                 100.0 * static_cast<double>(translationMove) / static_cast<double>(countTotal));
    }
    if (randomTranslationMove)
    {
      std::print(stream, "Random Translation:          {:14f} [%]\n",
                 100.0 * static_cast<double>(randomTranslationMove) / static_cast<double>(countTotal));
    }
    if (rotationMove)
    {
      std::print(stream, "Rotation:                    {:14f} [%]\n",
                 100.0 * static_cast<double>(rotationMove) / static_cast<double>(countTotal));
    }
    if (randomRotationMove)
    {
      std::print(stream, "Random rotation:             {:14f} [%]\n",
                 100.0 * static_cast<double>(randomRotationMove) / static_cast<double>(countTotal));
    }
    if (reinsertionMoveCBMC)
    {
      std::print(stream, "Reinsertion CBMC:            {:14f} [%]\n",
                 100.0 * static_cast<double>(reinsertionMoveCBMC) / static_cast<double>(countTotal));
    }
    if (swapInsertionMove)
    {
      std::print(stream, "Insertion:                   {:14f} [%]\n",
                 100.0 * static_cast<double>(swapInsertionMove) / static_cast<double>(countTotal));
    }
    if (swapDeletionMove)
    {
      std::print(stream, "Deletion:                    {:14f} [%]\n",
                 100.0 * static_cast<double>(swapDeletionMove) / static_cast<double>(countTotal));
    }
    if (swapInsertionMoveCBMC)
    {
      std::print(stream, "Insertion CBMC:              {:14f} [%]\n",
                 100.0 * static_cast<double>(swapInsertionMoveCBMC) / static_cast<double>(countTotal));
    }
    if (swapDeletionMoveCBMC)
    {
      std::print(stream, "Deletion CBMC:               {:14f} [%]\n",
                 100.0 * static_cast<double>(swapDeletionMoveCBMC) / static_cast<double>(countTotal));
    }
    if (swapLambdaMoveCFCMC)
    {
      std::print(stream, "Insertion/deletion CFCMC:    {:14f} [%]\n",
                 100.0 * static_cast<double>(swapLambdaMoveCFCMC) / static_cast<double>(countTotal));
    }
    if (swapLambdaMoveCBCFCMC)
    {
      std::print(stream, "Insertion/deletion CB/CFCMC: {:14f} [%]\n",
                 100.0 * static_cast<double>(swapLambdaMoveCBCFCMC) / static_cast<double>(countTotal));
    }
    if (GibbsSwapMoveCBMC)
    {
      std::print(stream, "Gibbs-swap CBMC:             {:14f} [%]\n",
                 100.0 * static_cast<double>(GibbsSwapMoveCBMC) / static_cast<double>(countTotal));
    }
    if (GibbsSwapLambdaMoveCFCMC)
    {
      std::print(stream, "Gibbs-swap CFCMC:            {:14f} [%]\n",
                 100.0 * static_cast<double>(GibbsSwapLambdaMoveCFCMC) / static_cast<double>(countTotal));
    }
    if (volumeMove)
    {
      std::print(stream, "Volume:                      {:14f} [%]\n",
                 100.0 * static_cast<double>(volumeMove) / static_cast<double>(countTotal));
    }
    if (WidomMoveCBMC)
    {
      std::print(stream, "Widom CBMC:                  {:14f} [%]\n",
                 100.0 * static_cast<double>(WidomMoveCBMC) / static_cast<double>(countTotal));
    }
    if (WidomMoveCFCMC)
    {
      std::print(stream, "Widom CFCMC:                 {:14f} [%]\n",
                 100.0 * static_cast<double>(WidomMoveCFCMC) / static_cast<double>(countTotal));
    }
    if (WidomMoveCBCFCMC)
    {
      std::print(stream, "Widom CB/CFCMC:              {:14f} [%]\n",
                 100.0 * static_cast<double>(WidomMoveCBCFCMC) / static_cast<double>(countTotal));
    }
    if (GibbsVolumeMove)
    {
      std::print(stream, "Gibbs Volume:                {:14f} [%]\n",
                 100.0 * static_cast<double>(GibbsVolumeMove) / static_cast<double>(countTotal));
    }
    if (ParallelTemperingSwap)
    {
      std::print(stream, "Parallel Tempering Swap:     {:14f} [%]\n",
                 100.0 * static_cast<double>(ParallelTemperingSwap) / static_cast<double>(countTotal));
    }

    std::print(stream, "\n");
    std::print(stream, "Production count MC-steps:   {:14d} [-]\n\n", countTotal);

    std::print(stream, "Translation:                 {:14d} [-]\n", translationMove);
    std::print(stream, "Random translation:          {:14d} [-]\n", randomTranslationMove);
    std::print(stream, "Rotation:                    {:14d} [-]\n", rotationMove);
    std::print(stream, "Random rotation:             {:14d} [-]\n", randomRotationMove);
    std::print(stream, "Reinsertion CBMC:            {:14d} [-]\n", reinsertionMoveCBMC);
    std::print(stream, "Insertion:                   {:14d} [-]\n", swapInsertionMove);
    std::print(stream, "Deletion:                    {:14d} [-]\n", swapDeletionMove);
    std::print(stream, "Insertion CBMC:              {:14d} [-]\n", swapInsertionMoveCBMC);
    std::print(stream, "Deletion CBMC:               {:14d} [-]\n", swapDeletionMoveCBMC);
    std::print(stream, "Insertion/deletion CFCMC:    {:14d} [-]\n", swapLambdaMoveCFCMC);
    std::print(stream, "Insertion/deletion CB/CFCMC: {:14d} [-]\n", swapLambdaMoveCBCFCMC);
    std::print(stream, "Gibbs-swap CBMC:             {:14d} [-]\n", GibbsSwapMoveCBMC);
    std::print(stream, "Gibbs-swap CFCMC:            {:14d} [-]\n", GibbsSwapLambdaMoveCFCMC);
    std::print(stream, "Widom CBMC:                  {:14d} [-]\n", WidomMoveCBMC);
    std::print(stream, "Widom CFCMC:                 {:14d} [-]\n", WidomMoveCFCMC);
    std::print(stream, "Widom CB/CFCMC:              {:14d} [-]\n", WidomMoveCBCFCMC);
    std::print(stream, "Volume:                      {:14d} [-]\n", volumeMove);
    std::print(stream, "Gibbs Volume:                {:14d} [-]\n", GibbsVolumeMove);
    std::print(stream, "Parallel Tempering Swap:     {:14d} [-]\n", ParallelTemperingSwap);
    std::print(stream, "                All summed:  {:14d} [-]\n", total());
    std::print(stream, "                difference:  {:14d} [-]\n", countTotal - total());

    std::print(stream, "\n\n");
  }

  return stream.str();
}

const nlohmann::json MCMoveCount::jsonAllSystemStatistics(size_t countTotal) const
{
  nlohmann::json percentage;
  nlohmann::json count;
  nlohmann::json status;

  if (translationMove || randomTranslationMove || rotationMove || randomRotationMove || reinsertionMoveCBMC ||
      swapInsertionMove || swapDeletionMove || swapInsertionMoveCBMC || swapDeletionMoveCBMC || swapLambdaMoveCFCMC ||
      swapLambdaMoveCBCFCMC || GibbsSwapLambdaMoveCFCMC || WidomMoveCBMC || WidomMoveCFCMC || WidomMoveCBCFCMC ||
      volumeMove || GibbsVolumeMove || ParallelTemperingSwap)
  {
    if (translationMove)
    {
      percentage["translation"] = 100.0 * static_cast<double>(translationMove) / static_cast<double>(countTotal);
    }
    if (randomTranslationMove)
    {
      percentage["randomTranslation"] =
          100.0 * static_cast<double>(randomTranslationMove) / static_cast<double>(countTotal);
    }
    if (rotationMove)
    {
      percentage["rotation"] = 100.0 * static_cast<double>(rotationMove) / static_cast<double>(countTotal);
    }
    if (randomRotationMove)
    {
      percentage["randomRotation"] = 100.0 * static_cast<double>(randomRotationMove) / static_cast<double>(countTotal);
    }
    if (reinsertionMoveCBMC)
    {
      percentage["reinsertionCBMC"] =
          100.0 * static_cast<double>(reinsertionMoveCBMC) / static_cast<double>(countTotal);
    }
    if (swapInsertionMove)
    {
      percentage["insertion"] = 100.0 * static_cast<double>(swapInsertionMove) / static_cast<double>(countTotal);
    }
    if (swapDeletionMove)
    {
      percentage["deletion"] = 100.0 * static_cast<double>(swapDeletionMove) / static_cast<double>(countTotal);
    }
    if (swapInsertionMoveCBMC)
    {
      percentage["insertionCBMC"] =
          100.0 * static_cast<double>(swapInsertionMoveCBMC) / static_cast<double>(countTotal);
    }
    if (swapDeletionMoveCBMC)
    {
      percentage["deletionCBMC"] = 100.0 * static_cast<double>(swapDeletionMoveCBMC) / static_cast<double>(countTotal);
    }
    if (swapLambdaMoveCFCMC)
    {
      percentage["insertionDeletionCBMC"] =
          100.0 * static_cast<double>(swapLambdaMoveCFCMC) / static_cast<double>(countTotal);
    }
    if (swapLambdaMoveCBCFCMC)
    {
      percentage["insertionDeletionCBCFCMC"] =
          100.0 * static_cast<double>(swapLambdaMoveCBCFCMC) / static_cast<double>(countTotal);
    }
    if (GibbsSwapMoveCBMC)
    {
      percentage["gibbsSwapCBMC"] = 100.0 * static_cast<double>(GibbsSwapMoveCBMC) / static_cast<double>(countTotal);
    }
    if (GibbsSwapLambdaMoveCFCMC)
    {
      percentage["gibbsSwapCFCMC"] =
          100.0 * static_cast<double>(GibbsSwapLambdaMoveCFCMC) / static_cast<double>(countTotal);
    }
    if (volumeMove)
    {
      percentage["volume"] = 100.0 * static_cast<double>(volumeMove) / static_cast<double>(countTotal);
    }
    if (WidomMoveCBMC)
    {
      percentage["widomCBMC"] = 100.0 * static_cast<double>(WidomMoveCBMC) / static_cast<double>(countTotal);
    }
    if (WidomMoveCFCMC)
    {
      percentage["widomCFCMC"] = 100.0 * static_cast<double>(WidomMoveCFCMC) / static_cast<double>(countTotal);
    }
    if (WidomMoveCBCFCMC)
    {
      percentage["widomCBCFCMC"] = 100.0 * static_cast<double>(WidomMoveCBCFCMC) / static_cast<double>(countTotal);
    }
    if (GibbsVolumeMove)
    {
      percentage["gibbsVolume"] = 100.0 * static_cast<double>(GibbsVolumeMove) / static_cast<double>(countTotal);
    }
    if (ParallelTemperingSwap)
    {
      percentage["parallelTemperingSwap"] =
          100.0 * static_cast<double>(ParallelTemperingSwap) / static_cast<double>(countTotal);
    }

    count["countTotal"] = countTotal;

    count["translation"] = translationMove;
    count["randomTranslation"] = randomTranslationMove;
    count["rotation"] = rotationMove;
    count["randomRotation"] = randomRotationMove;
    count["reinsertionCBMC"] = reinsertionMoveCBMC;
    count["insertion"] = swapInsertionMove;
    count["deletion"] = swapDeletionMove;
    count["insertionCBMC"] = swapInsertionMoveCBMC;
    count["deletionCBMC"] = swapDeletionMoveCBMC;
    count["insertionDeletionCFCMC"] = swapLambdaMoveCFCMC;
    count["InsertionDeletionCBCFCMC"] = swapLambdaMoveCBCFCMC;
    count["gibbsSwapCBMC"] = GibbsSwapMoveCBMC;
    count["gibbsSwapCFCMC"] = GibbsSwapLambdaMoveCFCMC;
    count["widomCBMC"] = WidomMoveCBMC;
    count["widomCFCMC"] = WidomMoveCFCMC;
    count["widomCBCFCMC"] = WidomMoveCBCFCMC;
    count["volume"] = volumeMove;
    count["gibbsVolume"] = GibbsVolumeMove;
    count["parallelTemperingSwap"] = ParallelTemperingSwap;
    count["sum"] = total();
    count["difference"] = countTotal - total();
  }

  status["percentage"] = percentage;
  status["count"] = count;
  return status;
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
  if (versionNumber > c.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
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
