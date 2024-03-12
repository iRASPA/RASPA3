export module mc_moves_cputime;

import <string>;
import <chrono>;
import <fstream>;

import double3;
import archive;

export struct MCMoveCpuTime
{
  MCMoveCpuTime()
  {
  };

  uint64_t versionNumber{ 1 };

  std::chrono::duration<double> propertySampling{ 0.0 };
  std::chrono::duration<double> energyPressureComputation{ 0.0 };

  std::chrono::duration<double> translationMove{ 0.0 };
  std::chrono::duration<double> translationMoveExternalFieldMolecule{ 0.0 };
  std::chrono::duration<double> translationMoveFrameworkMolecule{ 0.0 };
  std::chrono::duration<double> translationMoveMoleculeMolecule{ 0.0 };
  std::chrono::duration<double> translationMoveEwald{ 0.0 };

  std::chrono::duration<double> randomTranslationMove{ 0.0 };
  std::chrono::duration<double> randomTranslationMoveExternalFieldMolecule{ 0.0 };
  std::chrono::duration<double> randomTranslationMoveFrameworkMolecule{ 0.0 };
  std::chrono::duration<double> randomTranslationMoveMoleculeMolecule{ 0.0 };
  std::chrono::duration<double> randomTranslationMoveEwald{ 0.0 };

  std::chrono::duration<double> rotationMove{ 0.0 };
  std::chrono::duration<double> rotationMoveExternalFieldMolecule{ 0.0 };
  std::chrono::duration<double> rotationMoveFrameworkMolecule{ 0.0 };
  std::chrono::duration<double> rotationMoveMoleculeMolecule{ 0.0 };
  std::chrono::duration<double> rotationMoveEwald{ 0.0 };

  std::chrono::duration<double> randomRotationMove{ 0.0 };
  std::chrono::duration<double> randomRotationMoveExternalFieldMolecule{ 0.0 };
  std::chrono::duration<double> randomRotationMoveFrameworkMolecule{ 0.0 };
  std::chrono::duration<double> randomRotationMoveMoleculeMolecule{ 0.0 };
  std::chrono::duration<double> randomRotationMoveEwald{ 0.0 };

  std::chrono::duration<double> reinsertionMoveCBMC{ 0.0 };
  std::chrono::duration<double> reinsertionMoveCBMCNonEwald{ 0.0 };
  std::chrono::duration<double> reinsertionMoveCBMCEwald{ 0.0 };

  std::chrono::duration<double> swapInsertionMoveCBMC{ 0.0 };
  std::chrono::duration<double> swapInsertionMoveCBMCNonEwald{ 0.0 };
  std::chrono::duration<double> swapInsertionMoveCBMCEwald{ 0.0 };
  std::chrono::duration<double> swapInsertionMoveCBMCTail{ 0.0 };

  std::chrono::duration<double> swapDeletionMoveCBMC{ 0.0 };
  std::chrono::duration<double> swapDeletionMoveCBMCNonEwald{ 0.0 };
  std::chrono::duration<double> swapDeletionMoveCBMCEwald{ 0.0 };
  std::chrono::duration<double> swapDeletionMoveCBMCTail{ 0.0 };

  std::chrono::duration<double> swapLambdaDeletionMove{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveNonEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveEwald{ 0.0 };

  std::chrono::duration<double> swapLambdaMoveCFCMC{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCFramework{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCFCMCTail{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCFramework{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCFCMCTail{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCFramework{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCFCMCTail{ 0.0 };

  std::chrono::duration<double> swapLambdaMoveCBCFCMC{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCFramework{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCNonEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaInsertionMoveCBCFCMCTail{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCFramework{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaChangeMoveCBCFCMCTail{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCFramework{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCNonEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCEwald{ 0.0 };
  std::chrono::duration<double> swapLambdaDeletionMoveCBCFCMCTail{ 0.0 };

  std::chrono::duration<double> GibbsSwapMoveCBMC{ 0.0 };
  std::chrono::duration<double> GibbsSwapMoveCBMCNonEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapMoveCBMCEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapMoveCBMCTail{ 0.0 };

  std::chrono::duration<double> GibbsSwapLambdaMoveCFCMC{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaInterChangeMoveCFCMCNonEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaInterChangeMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaInterChangeMoveCFCMCTail{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaChangeMoveCFCMCNonEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaChangeMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaChangeMoveCFCMCTail{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaShuffleMoveCFCMCNonEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaShuffleMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> GibbsSwapLambdaShuffleMoveCFCMCTail{ 0.0 };

  std::chrono::duration<double> WidomMoveCBMC{ 0.0 };
  std::chrono::duration<double> WidomMoveCBMCNonEwald{ 0.0 };
  std::chrono::duration<double> WidomMoveCBMCEwald{ 0.0 };
  std::chrono::duration<double> WidomMoveCBMCTail{ 0.0 };

  std::chrono::duration<double> WidomMoveCFCMC{ 0.0 };
  std::chrono::duration<double> WidomMoveCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> WidomMoveCFCMCFramework{ 0.0 };
  std::chrono::duration<double> WidomMoveCFCMCMolecule{ 0.0 };
  std::chrono::duration<double> WidomMoveCFCMCEwald{ 0.0 };
  std::chrono::duration<double> WidomMoveCFCMCTail{ 0.0 };

  std::chrono::duration<double> WidomMoveCBCFCMC{ 0.0 };
  std::chrono::duration<double> WidomMoveCBCFCMCExternalField{ 0.0 };
  std::chrono::duration<double> WidomMoveCBCFCMCFramework{ 0.0 };
  std::chrono::duration<double> WidomMoveCBCFCMCMolecule{0.0};
  std::chrono::duration<double> WidomMoveCBCFCMCNonEwald{0.0};
  std::chrono::duration<double> WidomMoveCBCFCMCEwald{ 0.0 };
  std::chrono::duration<double> WidomMoveCBCFCMCTail{ 0.0 };

  std::chrono::duration<double> volumeMove{ 0.0 };
  std::chrono::duration<double> volumeMoveNonEwald{ 0.0 };
  std::chrono::duration<double> volumeMoveEwald{ 0.0 };
  std::chrono::duration<double> volumeMoveTail{ 0.0 };

  std::chrono::duration<double> GibbsVolumeMove{ 0.0 };
  std::chrono::duration<double> GibbsVolumeMoveNonEwald{ 0.0 };
  std::chrono::duration<double> GibbsVolumeMoveEwald{ 0.0 };
  std::chrono::duration<double> GibbsVolumeMoveTail{ 0.0 };

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
    translationMoveExternalFieldMolecule += b.translationMoveExternalFieldMolecule;
    translationMoveFrameworkMolecule += b.translationMoveFrameworkMolecule;
    translationMoveMoleculeMolecule += b.translationMoveMoleculeMolecule;
    translationMoveEwald += b.translationMoveEwald;

    randomTranslationMove += b.randomTranslationMove;
    randomTranslationMoveExternalFieldMolecule += b.randomTranslationMoveExternalFieldMolecule;
    randomTranslationMoveFrameworkMolecule += b.randomTranslationMoveFrameworkMolecule;
    randomTranslationMoveMoleculeMolecule += b.randomTranslationMoveMoleculeMolecule;
    randomTranslationMoveEwald += b.randomTranslationMoveEwald;

    rotationMove += b.rotationMove;
    rotationMoveExternalFieldMolecule += b.rotationMoveExternalFieldMolecule;
    rotationMoveFrameworkMolecule += b.rotationMoveFrameworkMolecule;
    rotationMoveMoleculeMolecule += b.rotationMoveMoleculeMolecule;
    rotationMoveEwald += b.rotationMoveEwald;

    randomRotationMove += b.randomRotationMove;
    randomRotationMoveExternalFieldMolecule += b.randomRotationMoveExternalFieldMolecule;
    randomRotationMoveFrameworkMolecule += b.randomRotationMoveFrameworkMolecule;
    randomRotationMoveMoleculeMolecule += b.randomRotationMoveMoleculeMolecule;
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
    swapLambdaInsertionMoveCFCMCExternalField += b.swapLambdaInsertionMoveCFCMCExternalField;
    swapLambdaInsertionMoveCFCMCFramework += b.swapLambdaInsertionMoveCFCMCFramework;
    swapLambdaInsertionMoveCFCMCMolecule += b.swapLambdaInsertionMoveCFCMCMolecule;
    swapLambdaInsertionMoveCFCMCEwald += b.swapLambdaInsertionMoveCFCMCEwald;
    swapLambdaInsertionMoveCFCMCTail += b.swapLambdaInsertionMoveCFCMCTail;
    swapLambdaChangeMoveCFCMCExternalField += b.swapLambdaChangeMoveCFCMCExternalField;
    swapLambdaChangeMoveCFCMCFramework += b.swapLambdaChangeMoveCFCMCFramework;
    swapLambdaChangeMoveCFCMCMolecule += b.swapLambdaChangeMoveCFCMCMolecule;
    swapLambdaChangeMoveCFCMCEwald += b.swapLambdaChangeMoveCFCMCEwald;
    swapLambdaChangeMoveCFCMCTail += b.swapLambdaChangeMoveCFCMCTail;
    swapLambdaDeletionMoveCFCMCExternalField += b.swapLambdaDeletionMoveCFCMCExternalField;
    swapLambdaDeletionMoveCFCMCFramework += b.swapLambdaDeletionMoveCFCMCFramework;
    swapLambdaDeletionMoveCFCMCMolecule += b.swapLambdaDeletionMoveCFCMCMolecule;
    swapLambdaDeletionMoveCFCMCEwald += b.swapLambdaDeletionMoveCFCMCEwald;
    swapLambdaDeletionMoveCFCMCTail += b.swapLambdaDeletionMoveCFCMCTail;

    swapLambdaMoveCBCFCMC += b.swapLambdaMoveCBCFCMC;
    swapLambdaInsertionMoveCBCFCMCExternalField += b.swapLambdaInsertionMoveCBCFCMCExternalField;
    swapLambdaInsertionMoveCBCFCMCFramework += b.swapLambdaInsertionMoveCBCFCMCFramework;
    swapLambdaInsertionMoveCBCFCMCMolecule += b.swapLambdaInsertionMoveCBCFCMCMolecule;
    swapLambdaInsertionMoveCBCFCMCNonEwald += b.swapLambdaInsertionMoveCBCFCMCNonEwald;
    swapLambdaInsertionMoveCBCFCMCEwald += b.swapLambdaInsertionMoveCBCFCMCEwald;
    swapLambdaInsertionMoveCBCFCMCTail += b.swapLambdaInsertionMoveCBCFCMCTail;
    swapLambdaChangeMoveCBCFCMCExternalField += b.swapLambdaChangeMoveCBCFCMCExternalField;
    swapLambdaChangeMoveCBCFCMCFramework += b.swapLambdaChangeMoveCBCFCMCFramework;
    swapLambdaChangeMoveCBCFCMCMolecule += b.swapLambdaChangeMoveCBCFCMCMolecule;
    swapLambdaChangeMoveCBCFCMCEwald += b.swapLambdaChangeMoveCBCFCMCEwald;
    swapLambdaChangeMoveCBCFCMCTail += b.swapLambdaChangeMoveCBCFCMCTail;
    swapLambdaDeletionMoveCBCFCMCExternalField += b.swapLambdaDeletionMoveCBCFCMCExternalField;
    swapLambdaDeletionMoveCBCFCMCFramework += b.swapLambdaDeletionMoveCBCFCMCFramework;
    swapLambdaDeletionMoveCBCFCMCMolecule += b.swapLambdaDeletionMoveCBCFCMCMolecule;
    swapLambdaDeletionMoveCBCFCMCNonEwald += b.swapLambdaDeletionMoveCBCFCMCNonEwald;
    swapLambdaDeletionMoveCBCFCMCEwald += b.swapLambdaDeletionMoveCBCFCMCEwald;
    swapLambdaDeletionMoveCBCFCMCTail += b.swapLambdaDeletionMoveCBCFCMCTail;

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
    WidomMoveCBMCTail += b.WidomMoveCBMCTail;

    WidomMoveCFCMC += b.WidomMoveCFCMC;
    WidomMoveCFCMCExternalField += b.WidomMoveCFCMCExternalField;
    WidomMoveCFCMCFramework += b.WidomMoveCFCMCFramework;
    WidomMoveCFCMCMolecule += b.WidomMoveCFCMCMolecule;
    WidomMoveCFCMCEwald += b.WidomMoveCFCMCEwald;
    WidomMoveCFCMCTail += b.WidomMoveCFCMCTail;

    WidomMoveCBCFCMC += b.WidomMoveCBCFCMC;
    WidomMoveCBCFCMCExternalField += b.WidomMoveCBCFCMCExternalField;
    WidomMoveCBCFCMCFramework += b.WidomMoveCBCFCMCFramework;
    WidomMoveCBCFCMCMolecule += b.WidomMoveCBCFCMCMolecule;
    WidomMoveCBCFCMCNonEwald += b.WidomMoveCBCFCMCNonEwald;
    WidomMoveCBCFCMCEwald += b.WidomMoveCBCFCMCEwald;
    WidomMoveCBCFCMCTail += b.WidomMoveCBCFCMCTail;

    volumeMove += b.volumeMove;
    volumeMoveNonEwald += b.volumeMoveNonEwald;
    volumeMoveEwald += b.volumeMoveEwald;
    volumeMoveTail += b.volumeMoveTail;

    GibbsVolumeMove += b.GibbsVolumeMove;
    GibbsVolumeMoveNonEwald += b.GibbsVolumeMoveNonEwald;
    GibbsVolumeMoveEwald += b.GibbsVolumeMoveEwald;
    GibbsVolumeMoveTail += b.GibbsVolumeMoveTail;

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
  m.translationMoveExternalFieldMolecule = a.translationMoveExternalFieldMolecule + b.translationMoveExternalFieldMolecule;
  m.translationMoveFrameworkMolecule = a.translationMoveFrameworkMolecule + b.translationMoveFrameworkMolecule;
  m.translationMoveMoleculeMolecule = a.translationMoveMoleculeMolecule + b.translationMoveMoleculeMolecule;
  m.translationMoveEwald = a.translationMoveEwald + b.translationMoveEwald;

  m.randomTranslationMove = a.randomTranslationMove + b.randomTranslationMove;
  m.randomTranslationMoveExternalFieldMolecule = a.randomTranslationMoveExternalFieldMolecule + b.randomTranslationMoveExternalFieldMolecule;
  m.randomTranslationMoveFrameworkMolecule = a.randomTranslationMoveFrameworkMolecule + b.randomTranslationMoveFrameworkMolecule;
  m.randomTranslationMoveMoleculeMolecule = a.randomTranslationMoveMoleculeMolecule + b.randomTranslationMoveMoleculeMolecule;
  m.randomTranslationMoveEwald = a.randomTranslationMoveEwald + b.randomTranslationMoveEwald;

  m.rotationMove = a.rotationMove + b.rotationMove;
  m.rotationMoveExternalFieldMolecule = a.rotationMoveExternalFieldMolecule + b.rotationMoveExternalFieldMolecule;
  m.rotationMoveFrameworkMolecule = a.rotationMoveFrameworkMolecule + b.rotationMoveFrameworkMolecule;
  m.rotationMoveMoleculeMolecule = a.rotationMoveMoleculeMolecule + b.rotationMoveMoleculeMolecule;
  m.rotationMoveEwald = a.rotationMoveEwald + b.rotationMoveEwald;

  m.randomRotationMove = a.randomRotationMove + b.randomRotationMove;
  m.randomRotationMoveExternalFieldMolecule = a.randomRotationMoveExternalFieldMolecule + b.randomRotationMoveExternalFieldMolecule;
  m.randomRotationMoveFrameworkMolecule = a.randomRotationMoveFrameworkMolecule + b.randomRotationMoveFrameworkMolecule;
  m.randomRotationMoveMoleculeMolecule = a.randomRotationMoveMoleculeMolecule + b.randomRotationMoveMoleculeMolecule;
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
  m.swapLambdaInsertionMoveCFCMCExternalField = a.swapLambdaInsertionMoveCFCMCExternalField + b.swapLambdaInsertionMoveCFCMCExternalField;
  m.swapLambdaInsertionMoveCFCMCFramework = a.swapLambdaInsertionMoveCFCMCFramework + b.swapLambdaInsertionMoveCFCMCFramework;
  m.swapLambdaInsertionMoveCFCMCMolecule = a.swapLambdaInsertionMoveCFCMCMolecule + b.swapLambdaInsertionMoveCFCMCMolecule;
  m.swapLambdaInsertionMoveCFCMCEwald = a.swapLambdaInsertionMoveCFCMCEwald + b.swapLambdaInsertionMoveCFCMCEwald;
  m.swapLambdaInsertionMoveCFCMCTail = a.swapLambdaInsertionMoveCFCMCTail + b.swapLambdaInsertionMoveCFCMCTail;
  m.swapLambdaChangeMoveCFCMCExternalField = a.swapLambdaChangeMoveCFCMCExternalField + b.swapLambdaChangeMoveCFCMCExternalField;
  m.swapLambdaChangeMoveCFCMCFramework = a.swapLambdaChangeMoveCFCMCFramework + b.swapLambdaChangeMoveCFCMCFramework;
  m.swapLambdaChangeMoveCFCMCMolecule = a.swapLambdaChangeMoveCFCMCMolecule + b.swapLambdaChangeMoveCFCMCMolecule;
  m.swapLambdaChangeMoveCFCMCEwald = a.swapLambdaChangeMoveCFCMCEwald + b.swapLambdaChangeMoveCFCMCEwald;
  m.swapLambdaChangeMoveCFCMCTail = a.swapLambdaChangeMoveCFCMCTail + b.swapLambdaChangeMoveCFCMCTail;
  m.swapLambdaDeletionMoveCFCMCExternalField = a.swapLambdaDeletionMoveCFCMCExternalField + b.swapLambdaDeletionMoveCFCMCExternalField;
  m.swapLambdaDeletionMoveCFCMCFramework = a.swapLambdaDeletionMoveCFCMCFramework + b.swapLambdaDeletionMoveCFCMCFramework;
  m.swapLambdaDeletionMoveCFCMCMolecule = a.swapLambdaDeletionMoveCFCMCMolecule + b.swapLambdaDeletionMoveCFCMCMolecule;
  m.swapLambdaDeletionMoveCFCMCEwald = a.swapLambdaDeletionMoveCFCMCEwald + b.swapLambdaDeletionMoveCFCMCEwald;
  m.swapLambdaDeletionMoveCFCMCTail = a.swapLambdaDeletionMoveCFCMCTail + b.swapLambdaDeletionMoveCFCMCTail;

  m.swapLambdaMoveCBCFCMC = a.swapLambdaMoveCBCFCMC + b.swapLambdaMoveCBCFCMC;
  m.swapLambdaInsertionMoveCBCFCMCExternalField = a.swapLambdaInsertionMoveCBCFCMCExternalField + b.swapLambdaInsertionMoveCBCFCMCExternalField;
  m.swapLambdaInsertionMoveCBCFCMCFramework = a.swapLambdaInsertionMoveCBCFCMCFramework + b.swapLambdaInsertionMoveCBCFCMCFramework;
  m.swapLambdaInsertionMoveCBCFCMCMolecule = a.swapLambdaInsertionMoveCBCFCMCMolecule + b.swapLambdaInsertionMoveCBCFCMCMolecule;
  m.swapLambdaInsertionMoveCBCFCMCNonEwald = a.swapLambdaInsertionMoveCBCFCMCNonEwald + b.swapLambdaInsertionMoveCBCFCMCNonEwald;
  m.swapLambdaInsertionMoveCBCFCMCEwald = a.swapLambdaInsertionMoveCBCFCMCEwald + b.swapLambdaInsertionMoveCBCFCMCEwald;
  m.swapLambdaInsertionMoveCBCFCMCTail = a.swapLambdaInsertionMoveCBCFCMCTail + b.swapLambdaInsertionMoveCBCFCMCTail;
  m.swapLambdaChangeMoveCBCFCMCExternalField = a.swapLambdaChangeMoveCBCFCMCExternalField + b.swapLambdaChangeMoveCBCFCMCExternalField;
  m.swapLambdaChangeMoveCBCFCMCFramework = a.swapLambdaChangeMoveCBCFCMCFramework + b.swapLambdaChangeMoveCBCFCMCFramework;
  m.swapLambdaChangeMoveCBCFCMCMolecule = a.swapLambdaChangeMoveCBCFCMCMolecule + b.swapLambdaChangeMoveCBCFCMCMolecule;
  m.swapLambdaChangeMoveCBCFCMCEwald = a.swapLambdaChangeMoveCBCFCMCEwald + b.swapLambdaChangeMoveCBCFCMCEwald;
  m.swapLambdaChangeMoveCBCFCMCTail = a.swapLambdaChangeMoveCBCFCMCTail + b.swapLambdaChangeMoveCBCFCMCTail;
  m.swapLambdaDeletionMoveCBCFCMCExternalField = a.swapLambdaDeletionMoveCBCFCMCExternalField + b.swapLambdaDeletionMoveCBCFCMCExternalField;
  m.swapLambdaDeletionMoveCBCFCMCFramework = a.swapLambdaDeletionMoveCBCFCMCFramework + b.swapLambdaDeletionMoveCBCFCMCFramework;
  m.swapLambdaDeletionMoveCBCFCMCMolecule = a.swapLambdaDeletionMoveCBCFCMCMolecule + b.swapLambdaDeletionMoveCBCFCMCMolecule;
  m.swapLambdaDeletionMoveCBCFCMCNonEwald = a.swapLambdaDeletionMoveCBCFCMCNonEwald + b.swapLambdaDeletionMoveCBCFCMCNonEwald;
  m.swapLambdaDeletionMoveCBCFCMCEwald = a.swapLambdaDeletionMoveCBCFCMCEwald + b.swapLambdaDeletionMoveCBCFCMCEwald;
  m.swapLambdaDeletionMoveCBCFCMCTail = a.swapLambdaDeletionMoveCBCFCMCTail + b.swapLambdaDeletionMoveCBCFCMCTail;

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
  m.WidomMoveCBMCTail = a.WidomMoveCBMCTail + b.WidomMoveCBMCTail;

  m.WidomMoveCFCMC = a.WidomMoveCFCMC + b.WidomMoveCFCMC;
  m.WidomMoveCFCMCExternalField = a.WidomMoveCFCMCExternalField + b.WidomMoveCFCMCExternalField;
  m.WidomMoveCFCMCFramework = a.WidomMoveCFCMCFramework + b.WidomMoveCFCMCFramework;
  m.WidomMoveCFCMCMolecule = a.WidomMoveCFCMCMolecule + b.WidomMoveCFCMCMolecule;
  m.WidomMoveCFCMCEwald = a.WidomMoveCFCMCEwald + b.WidomMoveCFCMCEwald;
  m.WidomMoveCFCMCTail = a.WidomMoveCFCMCTail + b.WidomMoveCFCMCTail;

  m.WidomMoveCBCFCMC = a.WidomMoveCBCFCMC + b.WidomMoveCBCFCMC;
  m.WidomMoveCBCFCMCExternalField = a.WidomMoveCBCFCMCExternalField + b.WidomMoveCBCFCMCExternalField;
  m.WidomMoveCBCFCMCFramework = a.WidomMoveCBCFCMCFramework + b.WidomMoveCBCFCMCFramework;
  m.WidomMoveCBCFCMCMolecule = a.WidomMoveCBCFCMCMolecule + b.WidomMoveCBCFCMCMolecule;
  m.WidomMoveCBCFCMCNonEwald = a.WidomMoveCBCFCMCNonEwald + b.WidomMoveCBCFCMCNonEwald;
  m.WidomMoveCBCFCMCEwald = a.WidomMoveCBCFCMCEwald + b.WidomMoveCBCFCMCEwald;
  m.WidomMoveCBCFCMCTail = a.WidomMoveCBCFCMCTail + b.WidomMoveCBCFCMCTail;

  m.volumeMove = a.volumeMove + b.volumeMove;
  m.volumeMoveNonEwald = a.volumeMoveNonEwald + b.volumeMoveNonEwald;
  m.volumeMoveEwald = a.volumeMoveEwald + b.volumeMoveEwald;
  m.volumeMoveTail = a.volumeMoveTail + b.volumeMoveTail;

  m.GibbsVolumeMove = a.GibbsVolumeMove + b.GibbsVolumeMove;
  m.GibbsVolumeMoveNonEwald = a.GibbsVolumeMoveNonEwald + b.GibbsVolumeMoveNonEwald;
  m.GibbsVolumeMoveEwald = a.GibbsVolumeMoveEwald + b.GibbsVolumeMoveEwald;
  m.GibbsVolumeMoveTail = a.GibbsVolumeMoveTail + b.GibbsVolumeMoveTail;

  return m;
}

