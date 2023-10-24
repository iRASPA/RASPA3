export module mc_moves_count;

import <string>;
import <chrono>;
import <fstream>;

import double3;
import archive;

export struct MCMoveCount
{
  MCMoveCount(): translationMove(0), randomTranslationMove(0), rotationMove(0), randomRotationMove(0),
                 reinsertionMoveCBMC(0), swapInsertionMoveCBMC(0), swapDeletionMoveCBMC(0), 
                 swapLambdaMoveCFCMC(0), swapLambdaMoveCBCFCMC(0), 
                 GibbsSwapMoveCBMC(0), GibbsSwapLambdaMoveCFCMC(0),
                 WidomMoveCBMC(0), WidomMoveCFCMC(0), WidomMoveCBCFCMC(0), 
                 volumeMove(0), GibbsVolumeMove(0)
                 {
                 };

  uint64_t versionNumber{ 1 };

  bool operator==(MCMoveCount const&) const = default;

  size_t translationMove;
  size_t randomTranslationMove;
  size_t rotationMove;
  size_t randomRotationMove;
  size_t reinsertionMoveCBMC;
  size_t swapInsertionMoveCBMC;
  size_t swapDeletionMoveCBMC;
  size_t swapLambdaDeletionMove;
  size_t swapLambdaMoveCFCMC;
  size_t swapLambdaMoveCBCFCMC;
  size_t GibbsSwapMoveCBMC;
  size_t GibbsSwapLambdaMoveCFCMC;
  size_t WidomMoveCBMC;
  size_t WidomMoveCFCMC;
  size_t WidomMoveCBCFCMC;
  size_t volumeMove;
  size_t GibbsVolumeMove;

  inline size_t total() const { return translationMove + randomTranslationMove +
                                       rotationMove + randomRotationMove +
                                       reinsertionMoveCBMC + swapInsertionMoveCBMC + swapDeletionMoveCBMC +
                                       swapLambdaMoveCFCMC + swapLambdaMoveCBCFCMC +
                                       GibbsSwapMoveCBMC + GibbsSwapLambdaMoveCFCMC +
                                       WidomMoveCBMC + WidomMoveCFCMC + WidomMoveCBCFCMC +
                                       volumeMove + GibbsVolumeMove;}
 
  void clearCountStatistics();
  const std::string writeMCMoveCountSystemStatistics(size_t countTotal) const;
  const std::string writeMCMoveCountComponentStatistics(size_t countTotal, size_t componentId, const std::string &componentName) const;
  const std::string writeMCMoveCountAllSystemStatistics(size_t countTotal) const;

  inline MCMoveCount& operator+=(const MCMoveCount& b)
  {
    translationMove += b.translationMove;
    randomTranslationMove += b.randomTranslationMove;
    rotationMove += b.rotationMove;
    randomRotationMove += b.randomRotationMove;
    reinsertionMoveCBMC += b.reinsertionMoveCBMC;
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

    return *this;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveCount &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveCount &c);
};

export inline MCMoveCount operator+(const MCMoveCount& a, const MCMoveCount& b)
{
  MCMoveCount m;

  m.translationMove = a.translationMove + b.translationMove;
  m.randomTranslationMove = a.randomTranslationMove + b.randomTranslationMove;
  m.rotationMove = a.rotationMove + b.rotationMove;
  m.randomRotationMove = a.randomRotationMove + b.randomRotationMove;
  m.reinsertionMoveCBMC = a.reinsertionMoveCBMC + b.reinsertionMoveCBMC;
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

  return m;
}

