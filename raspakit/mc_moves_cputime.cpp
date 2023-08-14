module mc_moves_cputime;

import <chrono>;
import <string>;
import <sstream>;

import double3;

import print;

void MCMoveCpuTime::clearTimingStatistics()
{
  translationMove = std::chrono::duration<double>::zero();
  translationMoveNonEwald  = std::chrono::duration<double>::zero();
  translationMoveEwald = std::chrono::duration<double>::zero();

  randomTranslationMove = std::chrono::duration<double>::zero();
  randomTranslationMoveNonEwald = std::chrono::duration<double>::zero();
  randomTranslationMoveEwald = std::chrono::duration<double>::zero();

  rotationMove = std::chrono::duration<double>::zero();
  rotationMoveNonEwald = std::chrono::duration<double>::zero();
  rotationMoveEwald = std::chrono::duration<double>::zero();

  randomRotationMove = std::chrono::duration<double>::zero();
  randomRotationMoveNonEwald = std::chrono::duration<double>::zero();
  randomRotationMoveEwald = std::chrono::duration<double>::zero();

  reinsertionMoveCBMC = std::chrono::duration<double>::zero();
  reinsertionMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  reinsertionMoveCBMCEwald = std::chrono::duration<double>::zero();

  swapInsertionMoveCBMC = std::chrono::duration<double>::zero();
  swapInsertionMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  swapInsertionMoveCBMCEwald = std::chrono::duration<double>::zero();
  swapInsertionMoveCBMCTail = std::chrono::duration<double>::zero();

  swapDeletionMoveCBMC = std::chrono::duration<double>::zero();
  swapDeletionMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  swapDeletionMoveCBMCEwald = std::chrono::duration<double>::zero();
  swapDeletionMoveCBMCTail = std::chrono::duration<double>::zero();

  swapLambdaMoveCFCMC = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCEwald = std::chrono::duration<double>::zero();

  swapLambdaMoveCBCFCMC = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCTail = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCEwald = std::chrono::duration<double>::zero();

  GibbsSwapMoveCBMC = std::chrono::duration<double>::zero();
  GibbsSwapMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapMoveCBMCEwald = std::chrono::duration<double>::zero();

  GibbsSwapLambdaMoveCFCMC = std::chrono::duration<double>::zero();
  GibbsSwapLambdaInterChangeMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaInterChangeMoveCFCMCEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaChangeMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaChangeMoveCFCMCEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaShuffleMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaShuffleMoveCFCMCEwald = std::chrono::duration<double>::zero();

  WidomMoveCBMC = std::chrono::duration<double>::zero();
  WidomMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  WidomMoveCBMCEwald = std::chrono::duration<double>::zero();

  WidomMoveCFCMC = std::chrono::duration<double>::zero();
  WidomMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  WidomMoveCFCMCEwald = std::chrono::duration<double>::zero();

  WidomMoveCBCFCMC = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCEwald = std::chrono::duration<double>::zero();

  volumeMove = std::chrono::duration<double>::zero();
  volumeMoveNonEwald = std::chrono::duration<double>::zero();
  volumeMoveEwald = std::chrono::duration<double>::zero();

  GibbsVolumeMove = std::chrono::duration<double>::zero();
  GibbsVolumeMoveNonEwald = std::chrono::duration<double>::zero();
  GibbsVolumeMoveEwald = std::chrono::duration<double>::zero();
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics() const
{
  std::ostringstream stream;

  if(volumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Volume:                          {:14f} [s]\n", volumeMove.count());
    std::print(stream, "    Non-Ewald:                   {:14f} [s]\n", volumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                       {:14f} [s]\n", volumeMoveEwald.count());
    std::print(stream, "    Overhead:                    {:14f} [s]\n",
                 volumeMove.count() - volumeMoveNonEwald.count() - volumeMoveEwald.count());
  }

  if(GibbsVolumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs Volume:                    {:14f} [s]\n", GibbsVolumeMove.count());
    std::print(stream, "    Non-Ewald:                   {:14f} [s]\n", GibbsVolumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                       {:14f} [s]\n", GibbsVolumeMoveEwald.count());
    std::print(stream, "    Overhead:                    {:14f} [s]\n",
                 GibbsVolumeMove.count() - GibbsVolumeMoveNonEwald.count() - GibbsVolumeMoveEwald.count());
  }

  std::print(stream, "\n");

  std::print(stream, "Property sampling                {:14f} [s]\n", propertySampling.count());
  std::print(stream, "Energy/pressure sampling:        {:14f} [s]\n", energyPressureComputation.count());

  std::print(stream, "\n\n");

  return stream.str();
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics(size_t componentId, const std::string &componentName) const
{
  std::ostringstream stream;

  std::print(stream, "Component {} {}\n", componentId, componentName);

  if(translationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Translation:                 {:14f} [s]\n", translationMove.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", translationMoveNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", translationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 translationMove.count() - translationMoveNonEwald.count() - translationMoveEwald.count());
  }

  if(randomTranslationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Random Translation:          {:14f} [s]\n", randomTranslationMove.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", randomTranslationMoveNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", randomTranslationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 randomTranslationMove.count() - randomTranslationMoveNonEwald.count() - randomTranslationMoveEwald.count());
  }

  if(rotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Rotation:                    {:14f} [s]\n", rotationMove.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", rotationMoveNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", rotationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 rotationMove.count() - rotationMoveNonEwald.count() - rotationMoveEwald.count());
  }

  if(randomRotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Random Rotation:             {:14f} [s]\n", randomRotationMove.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", randomRotationMoveNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", randomRotationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 randomRotationMove.count() - randomRotationMoveNonEwald.count() - randomRotationMoveEwald.count());
  }

  if(reinsertionMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Reinsertion CBMC:            {:14f} [s]\n", reinsertionMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", reinsertionMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", reinsertionMoveCBMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 reinsertionMoveCBMC.count() - reinsertionMoveCBMCNonEwald.count() - reinsertionMoveCBMCEwald.count());
  }

  if(swapInsertionMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Insertion CBMC:              {:14f} [s]\n", swapInsertionMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", swapInsertionMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", swapInsertionMoveCBMCEwald.count());
    std::print(stream, "        Tail:                    {:14f} [s]\n", swapInsertionMoveCBMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 swapInsertionMoveCBMC.count() - swapInsertionMoveCBMCNonEwald.count() - 
                 swapInsertionMoveCBMCEwald.count() - swapInsertionMoveCBMCTail.count());
  }

  if(swapDeletionMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Deletion CBMC:               {:14f} [s]\n", swapDeletionMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", swapDeletionMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", swapDeletionMoveCBMCEwald.count());
    std::print(stream, "        Tail:                    {:14f} [s]\n", swapDeletionMoveCBMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 swapDeletionMoveCBMC.count() - swapDeletionMoveCBMCNonEwald.count() - 
                 swapDeletionMoveCBMCEwald.count() - swapDeletionMoveCBMCTail.count());
  }

  if(swapLambdaMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Insertion/deletion CFCMC:    {:14f} [s]\n", swapLambdaMoveCFCMC.count());
    std::print(stream, "        Insertion Non-Ewald:     {:14f} [s]\n", swapLambdaInsertionMoveCFCMCNonEwald.count());
    std::print(stream, "        Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCFCMCEwald.count());
    std::print(stream, "        Lambda-change Non-Ewald: {:14f} [s]\n", swapLambdaChangeMoveCFCMCNonEwald.count());
    std::print(stream, "        Lambda-change Ewald:     {:14f} [s]\n", swapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "        Deletion Non-Ewald:      {:14f} [s]\n", swapLambdaDeletionMoveCFCMCNonEwald.count());
    std::print(stream, "        Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCFCMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n", swapLambdaMoveCFCMC.count()
                 - swapLambdaInsertionMoveCFCMCNonEwald.count() - swapLambdaInsertionMoveCFCMCEwald.count()
                 - swapLambdaChangeMoveCFCMCNonEwald.count() - swapLambdaChangeMoveCFCMCEwald.count()
                 - swapLambdaDeletionMoveCFCMCNonEwald.count() - swapLambdaDeletionMoveCFCMCEwald.count());
  }

  if(swapLambdaMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Insertion/deletion CB/CFCMC: {:14f} [s]\n", swapLambdaMoveCBCFCMC.count());
    std::print(stream, "        Insertion Non-Ewald:     {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCEwald.count());
    std::print(stream, "        Lambda-change Non-Ewald: {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Lambda-change Ewald:     {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCEwald.count());
    std::print(stream, "        Lambda-change Tail:      {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCTail.count());
    std::print(stream, "        Deletion Non-Ewald:      {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n", swapLambdaMoveCBCFCMC.count()
                 - swapLambdaInsertionMoveCBCFCMCNonEwald.count() - swapLambdaInsertionMoveCBCFCMCEwald.count()
                 - swapLambdaChangeMoveCBCFCMCNonEwald.count() - swapLambdaChangeMoveCBCFCMCEwald.count() - swapLambdaChangeMoveCBCFCMCTail.count()
                 - swapLambdaDeletionMoveCBCFCMCNonEwald.count() - swapLambdaDeletionMoveCBCFCMCEwald.count());
  }

  if(GibbsSwapMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Gibbs-swap CBMC:             {:14f} [s]\n", GibbsSwapMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", GibbsSwapMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", GibbsSwapMoveCBMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 GibbsSwapMoveCBMC.count() - GibbsSwapMoveCBMCNonEwald.count() - GibbsSwapMoveCBMCEwald.count());
  }

  if(GibbsSwapLambdaMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Gibbs-swap CFCMC:            {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count());
    std::print(stream, "        Inter-change Non-Ewald:  {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count());
    std::print(stream, "        Inter-change Ewald:      {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCEwald.count());
    std::print(stream, "        Lambda-change Non-Ewald: {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCNonEwald.count());
    std::print(stream, "        Lambda-change Ewald:     {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "        Shuffle Non-Ewald:       {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count());
    std::print(stream, "        Shuffle Ewald:           {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count()
                 - GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaInterChangeMoveCFCMCEwald.count()
                 - GibbsSwapLambdaChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaChangeMoveCFCMCEwald.count()
                 - GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count() - GibbsSwapLambdaShuffleMoveCFCMCEwald.count());
  }

  if(WidomMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Widom CBMC:                  {:14f} [s]\n", WidomMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", WidomMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", WidomMoveCBMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 WidomMoveCBMC.count() - WidomMoveCBMCNonEwald.count() - WidomMoveCBMCEwald.count());
  }

  if(WidomMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Widom CFCMC:                 {:14f} [s]\n", WidomMoveCFCMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", WidomMoveCFCMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", WidomMoveCFCMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 WidomMoveCFCMC.count() - WidomMoveCFCMCNonEwald.count() - WidomMoveCFCMCEwald.count());
  }

  if(WidomMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Widom CB/CFCMC:              {:14f} [s]\n", WidomMoveCBCFCMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", WidomMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", WidomMoveCBCFCMCEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 WidomMoveCBCFCMC.count() - WidomMoveCBCFCMCNonEwald.count() - WidomMoveCBCFCMCEwald.count());
  }

  std::print(stream, "\n\n");

  return stream.str();
}

const std::string MCMoveCpuTime::writeMCMoveCPUTimeStatistics(std::chrono::duration<double> totalSimulation) const
{
  std::ostringstream stream;

  if(translationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Translation:                 {:14f} [s]\n", translationMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", translationMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", translationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 translationMove.count() - translationMoveNonEwald.count() - translationMoveEwald.count());
  }

  if(randomTranslationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Random Translation:          {:14f} [s]\n", randomTranslationMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", randomTranslationMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", randomTranslationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 randomTranslationMove.count() - randomTranslationMoveNonEwald.count() - randomTranslationMoveEwald.count());
  }

  if(rotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Rotation:                    {:14f} [s]\n", rotationMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", rotationMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", rotationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 rotationMove.count() - rotationMoveNonEwald.count() - rotationMoveEwald.count());
  }

  if(randomRotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Random rotation:             {:14f} [s]\n", randomRotationMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", randomRotationMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", randomRotationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 randomRotationMove.count() - randomRotationMoveNonEwald.count() - randomRotationMoveEwald.count());
  }

  if(reinsertionMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Reinsertion CBMC:            {:14f} [s]\n", reinsertionMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", reinsertionMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", reinsertionMoveCBMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 reinsertionMoveCBMC.count() - reinsertionMoveCBMCNonEwald.count() - reinsertionMoveCBMCEwald.count());
  }

  if(swapInsertionMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Insertion CBMC:              {:14f} [s]\n", swapInsertionMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", swapInsertionMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", swapInsertionMoveCBMCEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", swapInsertionMoveCBMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 swapInsertionMoveCBMC.count() - swapInsertionMoveCBMCNonEwald.count() - 
                 swapInsertionMoveCBMCEwald.count() - swapInsertionMoveCBMCTail.count());
  }

  if(swapDeletionMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Deletion CBMC:               {:14f} [s]\n", swapDeletionMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", swapDeletionMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", swapDeletionMoveCBMCEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", swapDeletionMoveCBMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 swapDeletionMoveCBMC.count() - swapDeletionMoveCBMCNonEwald.count() - 
                 swapDeletionMoveCBMCEwald.count() - swapDeletionMoveCBMCTail.count());
  }

  if(swapLambdaMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Insertion/deletion CFCMC:    {:14f} [s]\n", swapLambdaMoveCFCMC.count());
    std::print(stream, "    Insertion Non-Ewald:     {:14f} [s]\n", swapLambdaInsertionMoveCFCMCNonEwald.count());
    std::print(stream, "    Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCFCMCEwald.count());
    std::print(stream, "    Lambda-change Non-Ewald: {:14f} [s]\n", swapLambdaChangeMoveCFCMCNonEwald.count());
    std::print(stream, "    Lambda-change Ewald:     {:14f} [s]\n", swapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "    Deletion Non-Ewald:      {:14f} [s]\n", swapLambdaDeletionMoveCFCMCNonEwald.count());
    std::print(stream, "    Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCFCMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n", swapLambdaMoveCFCMC.count()
                 - swapLambdaInsertionMoveCFCMCNonEwald.count() - swapLambdaInsertionMoveCFCMCEwald.count()
                 - swapLambdaChangeMoveCFCMCNonEwald.count() - swapLambdaChangeMoveCFCMCEwald.count()
                 - swapLambdaDeletionMoveCFCMCNonEwald.count() - swapLambdaDeletionMoveCFCMCEwald.count());
  }

  if(swapLambdaMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Insertion/deletion CB/CFCMC: {:14f} [s]\n", swapLambdaMoveCBCFCMC.count());
    std::print(stream, "    Insertion Non-Ewald:     {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCEwald.count());
    std::print(stream, "    Lambda-change Non-Ewald: {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Lambda-change Ewald:     {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCEwald.count());
    std::print(stream, "    Lambda-change Tail:      {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCTail.count());
    std::print(stream, "    Deletion Non-Ewald:      {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n", swapLambdaMoveCBCFCMC.count()
                 - swapLambdaInsertionMoveCBCFCMCNonEwald.count() - swapLambdaInsertionMoveCBCFCMCEwald.count()
                 - swapLambdaChangeMoveCBCFCMCNonEwald.count() - swapLambdaChangeMoveCBCFCMCEwald.count() - swapLambdaChangeMoveCBCFCMCTail.count()
                 - swapLambdaDeletionMoveCBCFCMCNonEwald.count() - swapLambdaDeletionMoveCBCFCMCEwald.count());
  }

  if(GibbsSwapMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs-swap CBMC:             {:14f} [s]\n", GibbsSwapMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", GibbsSwapMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", GibbsSwapMoveCBMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 GibbsSwapMoveCBMC.count() - GibbsSwapMoveCBMCNonEwald.count() - GibbsSwapMoveCBMCEwald.count());
  }

  if(GibbsSwapLambdaMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs-swap CFCMC:            {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count());
    std::print(stream, "    Inter-change Non-Ewald:  {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count());
    std::print(stream, "    Inter-change Ewald:      {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCEwald.count());
    std::print(stream, "    Lambda-change Non-Ewald: {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCNonEwald.count());
    std::print(stream, "    Lambda-change Ewald:     {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "    Shuffle Non-Ewald:       {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count());
    std::print(stream, "    Shuffle Ewald:           {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count()
                 - GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaInterChangeMoveCFCMCEwald.count()
                 - GibbsSwapLambdaChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaChangeMoveCFCMCEwald.count()
                 - GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count() - GibbsSwapLambdaShuffleMoveCFCMCEwald.count());
  }

  if(volumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Volume:                      {:14f} [s]\n", volumeMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", volumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", volumeMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 volumeMove.count() - volumeMoveNonEwald.count() - volumeMoveEwald.count());
  }

  if(WidomMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Widom CBMC:                  {:14f} [s]\n", WidomMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", WidomMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", WidomMoveCBMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 WidomMoveCBMC.count() - WidomMoveCBMCNonEwald.count() - WidomMoveCBMCEwald.count());
  }

  if(WidomMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Widom CFCMC:                 {:14f} [s]\n", WidomMoveCFCMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", WidomMoveCFCMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", WidomMoveCFCMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 WidomMoveCFCMC.count() - WidomMoveCFCMCNonEwald.count() - WidomMoveCFCMCEwald.count());
  }

  if(WidomMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Widom CB/CFCMC:              {:14f} [s]\n", WidomMoveCBCFCMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", WidomMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", WidomMoveCBCFCMCEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 WidomMoveCBCFCMC.count() - WidomMoveCBCFCMCNonEwald.count() - WidomMoveCBCFCMCEwald.count());
  }

  if(GibbsVolumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs Volume:                {:14f} [s]\n", GibbsVolumeMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", GibbsVolumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", GibbsVolumeMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 GibbsVolumeMove.count() - GibbsVolumeMoveNonEwald.count() - GibbsVolumeMoveEwald.count());
  }

  std::print(stream, "\n");
  std::print(stream, "Property sampling:           {:14f} [s]\n", propertySampling.count());
  std::print(stream, "Energy/pressure sampling:    {:14f} [s]\n\n", energyPressureComputation.count());

  std::print(stream, "Production simulation time:  {:14f} [s]\n", totalSimulation.count());
  std::print(stream, "                All summed:  {:14f} [s]\n", total().count());
  std::print(stream, "                  Overhead:  {:14f} [s]\n", totalSimulation.count() - total().count());

  std::print(stream, "\n\n");

  return stream.str();
}
