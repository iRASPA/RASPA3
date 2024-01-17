export module mc_moves_cputime;

import <string>;
import <chrono>;
import <fstream>;

import double3;
import archive;

export struct MCMoveCpuTime
{
  MCMoveCpuTime(): propertySampling(0.0), energyPressureComputation(0.0),
                   translationMove(0.0), translationMoveNonEwald(0.0), translationMoveEwald(0.0),
                   randomTranslationMove(0.0), randomTranslationMoveNonEwald(0.0), randomTranslationMoveEwald(0.0),
                   rotationMove(0.0), rotationMoveNonEwald(0.0), rotationMoveEwald(0.0), randomRotationMove(0.0),
                   randomRotationMoveNonEwald(0.0), randomRotationMoveEwald(0.0), reinsertionMoveCBMC(0.0),
                   reinsertionMoveCBMCNonEwald(0.0), reinsertionMoveCBMCEwald(0.0), swapInsertionMoveCBMC(0.0),
                   swapInsertionMoveCBMCNonEwald(0.0), swapInsertionMoveCBMCEwald(0.0), swapInsertionMoveCBMCTail(0.0),
                   swapDeletionMoveCBMC(0.0), swapDeletionMoveCBMCNonEwald(0.0), swapDeletionMoveCBMCEwald(0.0),
                   swapDeletionMoveCBMCTail(0.0),swapLambdaMoveCFCMC(0.0), swapLambdaInsertionMoveCFCMCNonEwald(0.0), 
                   swapLambdaInsertionMoveCFCMCEwald(0.0), swapLambdaChangeMoveCFCMCNonEwald(0.0), 
                   swapLambdaChangeMoveCFCMCEwald(0.0),  swapLambdaDeletionMoveCFCMCNonEwald(0.0), 
                   swapLambdaDeletionMoveCFCMCEwald(0.0), swapLambdaMoveCBCFCMC(0.0), 
                   swapLambdaInsertionMoveCBCFCMCNonEwald(0.0), swapLambdaInsertionMoveCBCFCMCEwald(0.0),
                   swapLambdaChangeMoveCBCFCMCNonEwald(0.0), swapLambdaChangeMoveCBCFCMCEwald(0.0), 
                   swapLambdaChangeMoveCBCFCMCTail(0.0), swapLambdaDeletionMoveCBCFCMCNonEwald(0.0), 
                   swapLambdaDeletionMoveCBCFCMCEwald(0.0), GibbsSwapMoveCBMC(0.0), GibbsSwapMoveCBMCNonEwald(0.0), 
                   GibbsSwapMoveCBMCEwald(0.0), GibbsSwapMoveCBMCTail(0.0),GibbsSwapLambdaMoveCFCMC(0.0),
                   GibbsSwapLambdaInterChangeMoveCFCMCNonEwald(0.0), GibbsSwapLambdaInterChangeMoveCFCMCEwald(0.0), 
                   GibbsSwapLambdaInterChangeMoveCFCMCTail(0.0), GibbsSwapLambdaChangeMoveCFCMCNonEwald(0.0), 
                   GibbsSwapLambdaChangeMoveCFCMCEwald(0.0), GibbsSwapLambdaChangeMoveCFCMCTail(0.0),
                   GibbsSwapLambdaShuffleMoveCFCMCNonEwald(0.0), GibbsSwapLambdaShuffleMoveCFCMCEwald(0.0), 
                   GibbsSwapLambdaShuffleMoveCFCMCTail(0.0), WidomMoveCBMC(0.0), WidomMoveCBMCNonEwald(0.0), 
                   WidomMoveCBMCEwald(0.0), WidomMoveCFCMC(0.0), WidomMoveCFCMCNonEwald(0.0), WidomMoveCFCMCEwald(0.0), 
                   WidomMoveCBCFCMC(0.0), WidomMoveCBCFCMCNonEwald(0.0), WidomMoveCBCFCMCEwald(0.0),
                   volumeMove(0.0), volumeMoveNonEwald(0.0), volumeMoveEwald(0.0),
                   GibbsVolumeMove(0.0), GibbsVolumeMoveNonEwald(0.0), GibbsVolumeMoveEwald(0.0)
                   {
                   };

  bool operator==(MCMoveCpuTime const&) const = default;

  uint64_t versionNumber{ 1 };

  std::chrono::duration<double> propertySampling;
  std::chrono::duration<double> energyPressureComputation;

  std::chrono::duration<double> translationMove;
  std::chrono::duration<double> translationMoveNonEwald;
  std::chrono::duration<double> translationMoveEwald;

  std::chrono::duration<double> randomTranslationMove;
  std::chrono::duration<double> randomTranslationMoveNonEwald;
  std::chrono::duration<double> randomTranslationMoveEwald;

  std::chrono::duration<double> rotationMove;
  std::chrono::duration<double> rotationMoveNonEwald;
  std::chrono::duration<double> rotationMoveEwald;

  std::chrono::duration<double> randomRotationMove;
  std::chrono::duration<double> randomRotationMoveNonEwald;
  std::chrono::duration<double> randomRotationMoveEwald;

  std::chrono::duration<double> reinsertionMoveCBMC;
  std::chrono::duration<double> reinsertionMoveCBMCNonEwald;
  std::chrono::duration<double> reinsertionMoveCBMCEwald;

  std::chrono::duration<double> swapInsertionMoveCBMC;
  std::chrono::duration<double> swapInsertionMoveCBMCNonEwald;
  std::chrono::duration<double> swapInsertionMoveCBMCEwald;
  std::chrono::duration<double> swapInsertionMoveCBMCTail;

  std::chrono::duration<double> swapDeletionMoveCBMC;
  std::chrono::duration<double> swapDeletionMoveCBMCNonEwald;
  std::chrono::duration<double> swapDeletionMoveCBMCEwald;
  std::chrono::duration<double> swapDeletionMoveCBMCTail;

  std::chrono::duration<double> swapLambdaDeletionMove;
  std::chrono::duration<double> swapLambdaDeletionMoveNonEwald;
  std::chrono::duration<double> swapLambdaDeletionMoveEwald;

  std::chrono::duration<double> swapLambdaMoveCFCMC;
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCNonEwald;
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCEwald;
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCNonEwald;
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCEwald;
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCNonEwald;
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCEwald;

  std::chrono::duration<double> swapLambdaMoveCBCFCMC;
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCNonEwald;
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCEwald;
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCNonEwald;
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCEwald;
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCTail;
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCNonEwald;
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCEwald;

  std::chrono::duration<double> GibbsSwapMoveCBMC;
  std::chrono::duration<double> GibbsSwapMoveCBMCNonEwald;
  std::chrono::duration<double> GibbsSwapMoveCBMCEwald;
  std::chrono::duration<double> GibbsSwapMoveCBMCTail;

  std::chrono::duration<double> GibbsSwapLambdaMoveCFCMC;
  std::chrono::duration<double> GibbsSwapLambdaInterChangeMoveCFCMCNonEwald;
  std::chrono::duration<double> GibbsSwapLambdaInterChangeMoveCFCMCEwald;
  std::chrono::duration<double> GibbsSwapLambdaInterChangeMoveCFCMCTail;
  std::chrono::duration<double> GibbsSwapLambdaChangeMoveCFCMCNonEwald;
  std::chrono::duration<double> GibbsSwapLambdaChangeMoveCFCMCEwald;
  std::chrono::duration<double> GibbsSwapLambdaChangeMoveCFCMCTail;
  std::chrono::duration<double> GibbsSwapLambdaShuffleMoveCFCMCNonEwald;
  std::chrono::duration<double> GibbsSwapLambdaShuffleMoveCFCMCEwald;
  std::chrono::duration<double> GibbsSwapLambdaShuffleMoveCFCMCTail;

  std::chrono::duration<double> WidomMoveCBMC;
  std::chrono::duration<double> WidomMoveCBMCNonEwald;
  std::chrono::duration<double> WidomMoveCBMCEwald;

  std::chrono::duration<double> WidomMoveCFCMC;
  std::chrono::duration<double> WidomMoveCFCMCNonEwald;
  std::chrono::duration<double> WidomMoveCFCMCEwald;

  std::chrono::duration<double> WidomMoveCBCFCMC;
  std::chrono::duration<double> WidomMoveCBCFCMCNonEwald;
  std::chrono::duration<double> WidomMoveCBCFCMCEwald;

  std::chrono::duration<double> volumeMove;
  std::chrono::duration<double> volumeMoveNonEwald;
  std::chrono::duration<double> volumeMoveEwald;

  std::chrono::duration<double> GibbsVolumeMove;
  std::chrono::duration<double> GibbsVolumeMoveNonEwald;
  std::chrono::duration<double> GibbsVolumeMoveEwald;

  inline std::chrono::duration<double> total() const { return propertySampling + energyPressureComputation +
                                                              translationMove + randomTranslationMove +
                                                              rotationMove + randomRotationMove +
                                                              reinsertionMoveCBMC + swapInsertionMoveCBMC + 
                                                              swapDeletionMoveCBMC + swapLambdaMoveCFCMC + 
                                                              swapLambdaMoveCBCFCMC + GibbsSwapMoveCBMC + 
                                                              GibbsSwapLambdaMoveCFCMC + WidomMoveCBMC + 
                                                              WidomMoveCFCMC + WidomMoveCBCFCMC +
                                                              volumeMove + GibbsVolumeMove;}
 
  void clearTimingStatistics();
  const std::string writeMCMoveCPUTimeStatistics() const;
  const std::string writeMCMoveCPUTimeStatistics(size_t componentId, const std::string &componentName) const;
  const std::string writeMCMoveCPUTimeStatistics(std::chrono::duration<double> total) const;

  inline MCMoveCpuTime& operator+=(const MCMoveCpuTime& b)
  {
    propertySampling += b.propertySampling;
    energyPressureComputation += b.energyPressureComputation;

    translationMove += b.translationMove;
    translationMoveNonEwald += b.translationMoveNonEwald;
    translationMoveEwald += b.translationMoveEwald;

    randomTranslationMove += b.randomTranslationMove;
    randomTranslationMoveNonEwald += b.randomTranslationMoveNonEwald;
    randomTranslationMoveEwald += b.randomTranslationMoveEwald;

    rotationMove += b.rotationMove;
    rotationMoveNonEwald += b.rotationMoveNonEwald;
    rotationMoveEwald += b.rotationMoveEwald;

    randomRotationMove += b.randomRotationMove;
    randomRotationMoveNonEwald += b.randomRotationMoveNonEwald;
    randomRotationMoveEwald += b.randomRotationMoveEwald;

    reinsertionMoveCBMC += b.reinsertionMoveCBMC;
    reinsertionMoveCBMCNonEwald += b.reinsertionMoveCBMCNonEwald;
    reinsertionMoveCBMCEwald += b.reinsertionMoveCBMCEwald;

    swapInsertionMoveCBMC += b.swapInsertionMoveCBMC;
    swapInsertionMoveCBMCNonEwald += b.swapInsertionMoveCBMCNonEwald;
    swapInsertionMoveCBMCEwald += b.swapInsertionMoveCBMCEwald;
    swapInsertionMoveCBMCTail += b.swapInsertionMoveCBMCTail;

    swapDeletionMoveCBMC += b.swapDeletionMoveCBMC;
    swapDeletionMoveCBMCNonEwald += b.swapDeletionMoveCBMCNonEwald;
    swapDeletionMoveCBMCEwald += b.swapDeletionMoveCBMCEwald;
    swapDeletionMoveCBMCTail += b.swapDeletionMoveCBMCTail;

    swapLambdaMoveCFCMC += b.swapLambdaMoveCFCMC;
    swapLambdaInsertionMoveCFCMCNonEwald += b.swapLambdaInsertionMoveCFCMCNonEwald;
    swapLambdaInsertionMoveCFCMCEwald += b.swapLambdaInsertionMoveCFCMCEwald;
    swapLambdaChangeMoveCFCMCNonEwald += b.swapLambdaChangeMoveCFCMCNonEwald;
    swapLambdaChangeMoveCFCMCEwald += b.swapLambdaChangeMoveCFCMCEwald;
    swapLambdaDeletionMoveCFCMCNonEwald += b.swapLambdaDeletionMoveCFCMCNonEwald;
    swapLambdaDeletionMoveCFCMCEwald += b.swapLambdaDeletionMoveCFCMCEwald;

    swapLambdaMoveCBCFCMC += b.swapLambdaMoveCBCFCMC;
    swapLambdaInsertionMoveCBCFCMCNonEwald += b.swapLambdaInsertionMoveCBCFCMCNonEwald;
    swapLambdaInsertionMoveCBCFCMCEwald += b.swapLambdaInsertionMoveCBCFCMCEwald;
    swapLambdaChangeMoveCBCFCMCNonEwald += b.swapLambdaChangeMoveCBCFCMCNonEwald;
    swapLambdaChangeMoveCBCFCMCEwald += b.swapLambdaChangeMoveCBCFCMCEwald;
    swapLambdaChangeMoveCBCFCMCTail += b.swapLambdaChangeMoveCBCFCMCTail;
    swapLambdaDeletionMoveCBCFCMCNonEwald += b.swapLambdaDeletionMoveCBCFCMCNonEwald;
    swapLambdaDeletionMoveCBCFCMCEwald += b.swapLambdaDeletionMoveCBCFCMCEwald;

    GibbsSwapMoveCBMC += b.GibbsSwapMoveCBMC;
    GibbsSwapMoveCBMCNonEwald += b.GibbsSwapMoveCBMCNonEwald;
    GibbsSwapMoveCBMCEwald += b.GibbsSwapMoveCBMCEwald;
    GibbsSwapMoveCBMCTail += b.GibbsSwapMoveCBMCTail;

    GibbsSwapLambdaMoveCFCMC += b.GibbsSwapLambdaMoveCFCMC;
    GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += b.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald;
    GibbsSwapLambdaInterChangeMoveCFCMCEwald += b.GibbsSwapLambdaInterChangeMoveCFCMCEwald;
    GibbsSwapLambdaInterChangeMoveCFCMCTail += b.GibbsSwapLambdaInterChangeMoveCFCMCTail;
    GibbsSwapLambdaChangeMoveCFCMCNonEwald += b.GibbsSwapLambdaChangeMoveCFCMCNonEwald;
    GibbsSwapLambdaChangeMoveCFCMCEwald += b.GibbsSwapLambdaChangeMoveCFCMCEwald;
    GibbsSwapLambdaChangeMoveCFCMCTail += b.GibbsSwapLambdaChangeMoveCFCMCTail;
    GibbsSwapLambdaShuffleMoveCFCMCNonEwald += b.GibbsSwapLambdaShuffleMoveCFCMCNonEwald;
    GibbsSwapLambdaShuffleMoveCFCMCEwald += b.GibbsSwapLambdaShuffleMoveCFCMCEwald;
    GibbsSwapLambdaShuffleMoveCFCMCTail += b.GibbsSwapLambdaShuffleMoveCFCMCTail;

    WidomMoveCBMC += b.WidomMoveCBMC;
    WidomMoveCBMCNonEwald += b.WidomMoveCBMCNonEwald;
    WidomMoveCBMCEwald += b.WidomMoveCBMCEwald;

    WidomMoveCFCMC += b.WidomMoveCFCMC;
    WidomMoveCFCMCNonEwald += b.WidomMoveCFCMCNonEwald;
    WidomMoveCFCMCEwald += b.WidomMoveCFCMCEwald;

    WidomMoveCBCFCMC += b.WidomMoveCBCFCMC;
    WidomMoveCBCFCMCNonEwald += b.WidomMoveCBCFCMCNonEwald;
    WidomMoveCBCFCMCEwald += b.WidomMoveCBCFCMCEwald;

    volumeMove += b.volumeMove;
    volumeMoveNonEwald += b.volumeMoveNonEwald;
    volumeMoveEwald += b.volumeMoveEwald;

    GibbsVolumeMove += b.GibbsVolumeMove;
    GibbsVolumeMoveNonEwald += b.GibbsVolumeMoveNonEwald;
    GibbsVolumeMoveEwald += b.GibbsVolumeMoveEwald;

    return *this;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveCpuTime &t);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveCpuTime &t);
};

export inline MCMoveCpuTime operator+(const MCMoveCpuTime& a, const MCMoveCpuTime& b)
{
  MCMoveCpuTime m;

  m.propertySampling = a.propertySampling + b.propertySampling;
  m.energyPressureComputation = a.energyPressureComputation + b.energyPressureComputation;

  m.translationMove = a.translationMove + b.translationMove;
  m.translationMoveNonEwald = a.translationMoveNonEwald + b.translationMoveNonEwald;
  m.translationMoveEwald = a.translationMoveEwald + b.translationMoveEwald;

  m.randomTranslationMove = a.randomTranslationMove + b.randomTranslationMove;
  m.randomTranslationMoveNonEwald = a.randomTranslationMoveNonEwald + b.randomTranslationMoveNonEwald;
  m.randomTranslationMoveEwald = a.randomTranslationMoveEwald + b.randomTranslationMoveEwald;

  m.rotationMove = a.rotationMove + b.rotationMove;
  m.rotationMoveNonEwald = a.rotationMoveNonEwald + b.rotationMoveNonEwald;
  m.rotationMoveEwald = a.rotationMoveEwald + b.rotationMoveEwald;

  m.randomRotationMove = a.randomRotationMove + b.randomRotationMove;
  m.randomRotationMoveNonEwald = a.randomRotationMoveNonEwald + b.randomRotationMoveNonEwald;
  m.randomRotationMoveEwald = a.randomRotationMoveEwald + b.randomRotationMoveEwald;

  m.reinsertionMoveCBMC = a.reinsertionMoveCBMC + b.reinsertionMoveCBMC;
  m.reinsertionMoveCBMCNonEwald = a.reinsertionMoveCBMCNonEwald + b.reinsertionMoveCBMCNonEwald;
  m.reinsertionMoveCBMCEwald = a.reinsertionMoveCBMCEwald + b.reinsertionMoveCBMCEwald;

  m.swapInsertionMoveCBMC = a.swapInsertionMoveCBMC + b.swapInsertionMoveCBMC;
  m.swapInsertionMoveCBMCNonEwald = a.swapInsertionMoveCBMCNonEwald + b.swapInsertionMoveCBMCNonEwald;
  m.swapInsertionMoveCBMCEwald = a.swapInsertionMoveCBMCEwald + b.swapInsertionMoveCBMCEwald;
  m.swapInsertionMoveCBMCTail = a.swapInsertionMoveCBMCEwald + b.swapInsertionMoveCBMCTail;

  m.swapDeletionMoveCBMC = a.swapDeletionMoveCBMC + b.swapDeletionMoveCBMC;
  m.swapDeletionMoveCBMCNonEwald = a.swapDeletionMoveCBMCNonEwald + b.swapDeletionMoveCBMCNonEwald;
  m.swapDeletionMoveCBMCEwald = a.swapDeletionMoveCBMCEwald + b.swapDeletionMoveCBMCEwald;
  m.swapDeletionMoveCBMCTail = a.swapDeletionMoveCBMCEwald + b.swapDeletionMoveCBMCTail;

  m.swapLambdaMoveCFCMC = a.swapLambdaMoveCFCMC + b.swapLambdaMoveCFCMC;
  m.swapLambdaInsertionMoveCFCMCNonEwald = 
    a.swapLambdaInsertionMoveCFCMCNonEwald + b.swapLambdaInsertionMoveCFCMCNonEwald;
  m.swapLambdaInsertionMoveCFCMCEwald = a.swapLambdaInsertionMoveCFCMCEwald + b.swapLambdaInsertionMoveCFCMCEwald;
  m.swapLambdaChangeMoveCFCMCNonEwald = a.swapLambdaChangeMoveCFCMCNonEwald + b.swapLambdaChangeMoveCFCMCNonEwald;
  m.swapLambdaChangeMoveCFCMCEwald = a.swapLambdaChangeMoveCFCMCEwald + b.swapLambdaChangeMoveCFCMCEwald;
  m.swapLambdaDeletionMoveCFCMCNonEwald = 
   a.swapLambdaDeletionMoveCFCMCNonEwald + b.swapLambdaDeletionMoveCFCMCNonEwald;
  m.swapLambdaDeletionMoveCFCMCEwald = a.swapLambdaDeletionMoveCFCMCEwald + b.swapLambdaDeletionMoveCFCMCEwald;

  m.swapLambdaMoveCBCFCMC = a.swapLambdaMoveCBCFCMC + b.swapLambdaMoveCBCFCMC;
  m.swapLambdaInsertionMoveCBCFCMCNonEwald = 
    a.swapLambdaInsertionMoveCBCFCMCNonEwald + b.swapLambdaInsertionMoveCBCFCMCNonEwald;
  m.swapLambdaInsertionMoveCBCFCMCEwald = 
    a.swapLambdaInsertionMoveCBCFCMCEwald + b.swapLambdaInsertionMoveCBCFCMCEwald;
  m.swapLambdaChangeMoveCBCFCMCNonEwald = 
    a.swapLambdaChangeMoveCBCFCMCNonEwald + b.swapLambdaChangeMoveCBCFCMCNonEwald;
  m.swapLambdaChangeMoveCBCFCMCEwald = a.swapLambdaChangeMoveCBCFCMCEwald + b.swapLambdaChangeMoveCBCFCMCEwald;
  m.swapLambdaChangeMoveCBCFCMCTail = a.swapLambdaChangeMoveCBCFCMCTail + b.swapLambdaChangeMoveCBCFCMCTail;
  m.swapLambdaDeletionMoveCBCFCMCNonEwald = 
    a.swapLambdaDeletionMoveCBCFCMCNonEwald + b.swapLambdaDeletionMoveCBCFCMCNonEwald;
  m.swapLambdaDeletionMoveCBCFCMCEwald = a.swapLambdaDeletionMoveCBCFCMCEwald + b.swapLambdaDeletionMoveCBCFCMCEwald;

  m.GibbsSwapMoveCBMC = a.GibbsSwapMoveCBMC + b.GibbsSwapMoveCBMC;
  m.GibbsSwapMoveCBMCNonEwald = a.GibbsSwapMoveCBMCNonEwald + b.GibbsSwapMoveCBMCNonEwald;
  m.GibbsSwapMoveCBMCEwald = a.GibbsSwapMoveCBMCEwald + b.GibbsSwapMoveCBMCEwald;
  m.GibbsSwapMoveCBMCTail = a.GibbsSwapMoveCBMCTail + b.GibbsSwapMoveCBMCTail;

  m.GibbsSwapLambdaMoveCFCMC = a.GibbsSwapLambdaMoveCFCMC + b.GibbsSwapLambdaMoveCFCMC;
  m.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald = 
    a.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald + b.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald;
  m.GibbsSwapLambdaInterChangeMoveCFCMCEwald = 
    a.GibbsSwapLambdaInterChangeMoveCFCMCEwald + b.GibbsSwapLambdaInterChangeMoveCFCMCEwald;
  m.GibbsSwapLambdaInterChangeMoveCFCMCTail = 
    a.GibbsSwapLambdaInterChangeMoveCFCMCTail + b.GibbsSwapLambdaInterChangeMoveCFCMCTail;
  m.GibbsSwapLambdaChangeMoveCFCMCNonEwald = 
    a.GibbsSwapLambdaChangeMoveCFCMCNonEwald + b.GibbsSwapLambdaChangeMoveCFCMCNonEwald;
  m.GibbsSwapLambdaChangeMoveCFCMCEwald = 
    a.GibbsSwapLambdaChangeMoveCFCMCEwald + b.GibbsSwapLambdaChangeMoveCFCMCEwald;
  m.GibbsSwapLambdaChangeMoveCFCMCTail = 
    a.GibbsSwapLambdaChangeMoveCFCMCTail + b.GibbsSwapLambdaChangeMoveCFCMCTail;
  m.GibbsSwapLambdaShuffleMoveCFCMCNonEwald = 
    a.GibbsSwapLambdaShuffleMoveCFCMCNonEwald + b.GibbsSwapLambdaShuffleMoveCFCMCNonEwald;
  m.GibbsSwapLambdaShuffleMoveCFCMCEwald = 
    a.GibbsSwapLambdaShuffleMoveCFCMCEwald + b.GibbsSwapLambdaShuffleMoveCFCMCEwald;
  m.GibbsSwapLambdaShuffleMoveCFCMCTail = 
    a.GibbsSwapLambdaShuffleMoveCFCMCTail + b.GibbsSwapLambdaShuffleMoveCFCMCTail;

  m.WidomMoveCBMC = a.WidomMoveCBMC + b.WidomMoveCBMC;
  m.WidomMoveCBMCNonEwald = a.WidomMoveCBMCNonEwald + b.WidomMoveCBMCNonEwald;
  m.WidomMoveCBMCEwald = a.WidomMoveCBMCEwald + b.WidomMoveCBMCEwald;

  m.WidomMoveCFCMC = a.WidomMoveCFCMC + b.WidomMoveCFCMC;
  m.WidomMoveCFCMCNonEwald = a.WidomMoveCFCMCNonEwald + b.WidomMoveCFCMCNonEwald;
  m.WidomMoveCFCMCEwald = a.WidomMoveCFCMCEwald + b.WidomMoveCFCMCEwald;

  m.WidomMoveCBCFCMC = a.WidomMoveCBCFCMC + b.WidomMoveCBCFCMC;
  m.WidomMoveCBCFCMCNonEwald = a.WidomMoveCBCFCMCNonEwald + b.WidomMoveCBCFCMCNonEwald;
  m.WidomMoveCBCFCMCEwald = a.WidomMoveCBCFCMCEwald + b.WidomMoveCBCFCMCEwald;

  m.volumeMove = a.volumeMove + b.volumeMove;
  m.volumeMoveNonEwald = a.volumeMoveNonEwald + b.volumeMoveNonEwald;
  m.volumeMoveEwald = a.volumeMoveEwald + b.volumeMoveEwald;

  m.GibbsVolumeMove = a.GibbsVolumeMove + b.GibbsVolumeMove;
  m.GibbsVolumeMoveNonEwald = a.GibbsVolumeMoveNonEwald + b.GibbsVolumeMoveNonEwald;
  m.GibbsVolumeMoveEwald = a.GibbsVolumeMoveEwald + b.GibbsVolumeMoveEwald;

  return m;
}

