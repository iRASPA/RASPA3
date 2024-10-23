module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_count;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import double3;
import archive;
import json;

/**
 * \brief Struct to keep track of Monte Carlo move counts.
 *
 * The MCMoveCount struct records the number of times each type of Monte Carlo move
 * has been performed during a simulation. It includes counts for various move types
 * such as translation, rotation, insertion, deletion, volume changes, and more.
 * This struct provides methods to clear the counts, write statistics, and accumulate counts.
 */
export struct MCMoveCount
{
  /**
   * \brief Default constructor for MCMoveCount.
   *
   * Initializes all move counts to zero.
   */
  MCMoveCount()
      : translationMove(0),
        randomTranslationMove(0),
        rotationMove(0),
        randomRotationMove(0),
        reinsertionMoveCBMC(0),
        swapInsertionMove(0),
        swapDeletionMove(0),
        swapInsertionMoveCBMC(0),
        swapDeletionMoveCBMC(0),
        swapLambdaMoveCFCMC(0),
        swapLambdaMoveCBCFCMC(0),
        GibbsSwapMoveCBMC(0),
        GibbsSwapLambdaMoveCFCMC(0),
        WidomMoveCBMC(0),
        WidomMoveCFCMC(0),
        WidomMoveCBCFCMC(0),
        volumeMove(0),
        GibbsVolumeMove(0),
        ParallelTemperingSwap(0) {};

  uint64_t versionNumber{1};  ///< Version number for serialization purposes.

  bool operator==(MCMoveCount const &) const = default;

  size_t translationMove;           ///< Count of translation moves performed.
  size_t randomTranslationMove;     ///< Count of random translation moves performed.
  size_t rotationMove;              ///< Count of rotation moves performed.
  size_t randomRotationMove;        ///< Count of random rotation moves performed.
  size_t reinsertionMoveCBMC;       ///< Count of CBMC reinsertion moves performed.
  size_t swapInsertionMove;         ///< Count of insertion moves performed.
  size_t swapDeletionMove;          ///< Count of deletion moves performed.
  size_t swapInsertionMoveCBMC;     ///< Count of CBMC insertion moves performed.
  size_t swapDeletionMoveCBMC;      ///< Count of CBMC deletion moves performed.
  size_t swapLambdaMoveCFCMC;       ///< Count of CFCMC swap lambda moves performed.
  size_t swapLambdaMoveCBCFCMC;     ///< Count of CB/CFCMC swap lambda moves performed.
  size_t GibbsSwapMoveCBMC;         ///< Count of CBMC Gibbs swap moves performed.
  size_t GibbsSwapLambdaMoveCFCMC;  ///< Count of CFCMC Gibbs swap lambda moves performed.
  size_t WidomMoveCBMC;             ///< Count of CBMC Widom moves performed.
  size_t WidomMoveCFCMC;            ///< Count of CFCMC Widom moves performed.
  size_t WidomMoveCBCFCMC;          ///< Count of CB/CFCMC Widom moves performed.
  size_t volumeMove;                ///< Count of volume moves performed.
  size_t GibbsVolumeMove;           ///< Count of Gibbs volume moves performed.
  size_t ParallelTemperingSwap;     ///< Count of parallel tempering swap moves performed.

  /**
   * \brief Calculates the total number of moves performed.
   *
   * Sums all individual move counts to provide the total number of moves.
   *
   * \return The total number of moves performed.
   */
  inline size_t total() const
  {
    return translationMove + randomTranslationMove + rotationMove + randomRotationMove + reinsertionMoveCBMC +
           swapInsertionMove + swapDeletionMove + swapInsertionMoveCBMC + swapDeletionMoveCBMC + swapLambdaMoveCFCMC +
           swapLambdaMoveCBCFCMC + GibbsSwapMoveCBMC + GibbsSwapLambdaMoveCFCMC + WidomMoveCBMC + WidomMoveCFCMC +
           WidomMoveCBCFCMC + volumeMove + GibbsVolumeMove + ParallelTemperingSwap;
  }

  /**
   * \brief Resets all move counts to zero.
   *
   * Clears the statistics by setting all move counts to zero.
   */
  void clearCountStatistics();
  /**
   * \brief Generates a string with formatted move statistics.
   *
   * Provides a string representation of move counts and their percentages
   * based on the total number of moves.
   *
   * \param countTotal The total number of moves performed.
   * \return A string containing the formatted move statistics.
   */
  const std::string writeAllSystemStatistics(size_t countTotal) const;
  /**
   * \brief Generates a JSON object with move statistics.
   *
   * Provides a JSON representation of move counts and their percentages
   * based on the total number of moves.
   *
   * \param countTotal The total number of moves performed.
   * \return A JSON object containing the move statistics.
   */
  const nlohmann::json jsonAllSystemStatistics(size_t countTotal) const;

  /**
   * \brief Accumulates move counts from another MCMoveCount instance.
   *
   * Adds the move counts from another MCMoveCount object to this one.
   *
   * \param b Another MCMoveCount object whose counts will be added.
   * \return A reference to this MCMoveCount object after addition.
   */
  inline MCMoveCount &operator+=(const MCMoveCount &b)
  {
    translationMove += b.translationMove;
    randomTranslationMove += b.randomTranslationMove;
    rotationMove += b.rotationMove;
    randomRotationMove += b.randomRotationMove;
    reinsertionMoveCBMC += b.reinsertionMoveCBMC;
    swapInsertionMove += b.swapInsertionMove;
    swapDeletionMove += b.swapDeletionMove;
    swapInsertionMoveCBMC += b.swapInsertionMoveCBMC;
    swapDeletionMoveCBMC += b.swapDeletionMoveCBMC;
    swapLambdaMoveCFCMC += b.swapLambdaMoveCFCMC;
    swapLambdaMoveCBCFCMC += b.swapLambdaMoveCBCFCMC;
    GibbsSwapMoveCBMC += b.GibbsSwapMoveCBMC;
    GibbsSwapLambdaMoveCFCMC += b.GibbsSwapLambdaMoveCFCMC;
    WidomMoveCBMC += b.WidomMoveCBMC;
    WidomMoveCFCMC += b.WidomMoveCFCMC;
    WidomMoveCBCFCMC += b.WidomMoveCBCFCMC;
    volumeMove += b.volumeMove;
    GibbsVolumeMove += b.GibbsVolumeMove;
    ParallelTemperingSwap += b.ParallelTemperingSwap;

    return *this;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveCount &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveCount &c);
};

export inline MCMoveCount operator+(const MCMoveCount &a, const MCMoveCount &b)
{
  MCMoveCount m;

  m.translationMove = a.translationMove + b.translationMove;
  m.randomTranslationMove = a.randomTranslationMove + b.randomTranslationMove;
  m.rotationMove = a.rotationMove + b.rotationMove;
  m.randomRotationMove = a.randomRotationMove + b.randomRotationMove;
  m.reinsertionMoveCBMC = a.reinsertionMoveCBMC + b.reinsertionMoveCBMC;
  m.swapInsertionMove = a.swapInsertionMove + b.swapInsertionMove;
  m.swapDeletionMove = a.swapDeletionMove + b.swapDeletionMove;
  m.swapInsertionMoveCBMC = a.swapInsertionMoveCBMC + b.swapInsertionMoveCBMC;
  m.swapDeletionMoveCBMC = a.swapDeletionMoveCBMC + b.swapDeletionMoveCBMC;
  m.swapLambdaMoveCFCMC = a.swapLambdaMoveCFCMC + b.swapLambdaMoveCFCMC;
  m.swapLambdaMoveCBCFCMC = a.swapLambdaMoveCBCFCMC + b.swapLambdaMoveCBCFCMC;
  m.GibbsSwapMoveCBMC = a.GibbsSwapMoveCBMC + b.GibbsSwapMoveCBMC;
  m.GibbsSwapLambdaMoveCFCMC = a.GibbsSwapLambdaMoveCFCMC + b.GibbsSwapLambdaMoveCFCMC;
  m.WidomMoveCBMC = a.WidomMoveCBMC + b.WidomMoveCBMC;
  m.WidomMoveCFCMC = a.WidomMoveCFCMC + b.WidomMoveCFCMC;
  m.WidomMoveCBCFCMC = a.WidomMoveCBCFCMC + b.WidomMoveCBCFCMC;
  m.volumeMove = a.volumeMove + b.volumeMove;
  m.GibbsVolumeMove = a.GibbsVolumeMove + b.GibbsVolumeMove;
  m.ParallelTemperingSwap = a.ParallelTemperingSwap + b.ParallelTemperingSwap;

  return m;
}
