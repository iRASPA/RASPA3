module mc_moves_cputime;

import <chrono>;
import <string>;
import <sstream>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <fstream>;
import <complex>;

import double3;
import stringutils;
import archive;


void MCMoveCpuTime::clearTimingStatistics()
{
  translationMove = std::chrono::duration<double>::zero();
  translationMoveExternalFieldMolecule = std::chrono::duration<double>::zero();
  translationMoveFrameworkMolecule = std::chrono::duration<double>::zero();
  translationMoveMoleculeMolecule = std::chrono::duration<double>::zero();
  translationMoveEwald = std::chrono::duration<double>::zero();

  randomTranslationMove = std::chrono::duration<double>::zero();
  randomTranslationMoveExternalFieldMolecule = std::chrono::duration<double>::zero();
  randomTranslationMoveFrameworkMolecule = std::chrono::duration<double>::zero();
  randomTranslationMoveMoleculeMolecule = std::chrono::duration<double>::zero();
  randomTranslationMoveEwald = std::chrono::duration<double>::zero();

  rotationMove = std::chrono::duration<double>::zero();
  rotationMoveExternalFieldMolecule = std::chrono::duration<double>::zero();
  rotationMoveFrameworkMolecule = std::chrono::duration<double>::zero();
  rotationMoveMoleculeMolecule = std::chrono::duration<double>::zero();
  rotationMoveEwald = std::chrono::duration<double>::zero();

  randomRotationMove = std::chrono::duration<double>::zero();
  randomRotationMoveExternalFieldMolecule = std::chrono::duration<double>::zero();
  randomRotationMoveFrameworkMolecule = std::chrono::duration<double>::zero();
  randomRotationMoveMoleculeMolecule = std::chrono::duration<double>::zero();
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
  swapLambdaInsertionMoveCFCMCExternalField = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCFCMCFramework = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCFCMCMolecule = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCFCMCTail = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCExternalField = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCFramework = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCMolecule = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCFCMCTail = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCExternalField = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCFramework = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCMolecule = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCFCMCTail = std::chrono::duration<double>::zero();

  swapLambdaMoveCBCFCMC = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCExternalField = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCFramework = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCMolecule = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaInsertionMoveCBCFCMCTail = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCExternalField = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCFramework = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCMolecule = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaChangeMoveCBCFCMCTail = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCExternalField = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCFramework = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCMolecule = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCEwald = std::chrono::duration<double>::zero();
  swapLambdaDeletionMoveCBCFCMCTail = std::chrono::duration<double>::zero();

  GibbsSwapMoveCBMC = std::chrono::duration<double>::zero();
  GibbsSwapMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapMoveCBMCEwald = std::chrono::duration<double>::zero();
  GibbsSwapMoveCBMCTail = std::chrono::duration<double>::zero();

  GibbsSwapLambdaMoveCFCMC = std::chrono::duration<double>::zero();
  GibbsSwapLambdaInterChangeMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaInterChangeMoveCFCMCEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaInterChangeMoveCFCMCTail = std::chrono::duration<double>::zero();
  GibbsSwapLambdaChangeMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaChangeMoveCFCMCEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaChangeMoveCFCMCTail = std::chrono::duration<double>::zero();
  GibbsSwapLambdaShuffleMoveCFCMCNonEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaShuffleMoveCFCMCEwald = std::chrono::duration<double>::zero();
  GibbsSwapLambdaShuffleMoveCFCMCTail = std::chrono::duration<double>::zero();

  WidomMoveCBMC = std::chrono::duration<double>::zero();
  WidomMoveCBMCNonEwald = std::chrono::duration<double>::zero();
  WidomMoveCBMCEwald = std::chrono::duration<double>::zero();
  WidomMoveCBMCTail = std::chrono::duration<double>::zero();

  WidomMoveCFCMC = std::chrono::duration<double>::zero();
  WidomMoveCFCMCExternalField = std::chrono::duration<double>::zero();
  WidomMoveCFCMCFramework = std::chrono::duration<double>::zero();
  WidomMoveCFCMCMolecule = std::chrono::duration<double>::zero();
  WidomMoveCFCMCEwald = std::chrono::duration<double>::zero();
  WidomMoveCFCMCTail = std::chrono::duration<double>::zero();

  WidomMoveCBCFCMC = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCExternalField = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCFramework = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCMolecule = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCNonEwald = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCEwald = std::chrono::duration<double>::zero();
  WidomMoveCBCFCMCTail = std::chrono::duration<double>::zero();

  volumeMove = std::chrono::duration<double>::zero();
  volumeMoveNonEwald = std::chrono::duration<double>::zero();
  volumeMoveEwald = std::chrono::duration<double>::zero();
  volumeMoveTail = std::chrono::duration<double>::zero();

  GibbsVolumeMove = std::chrono::duration<double>::zero();
  GibbsVolumeMoveNonEwald = std::chrono::duration<double>::zero();
  GibbsVolumeMoveEwald = std::chrono::duration<double>::zero();
  GibbsVolumeMoveTail = std::chrono::duration<double>::zero();
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
    std::print(stream, "    Tail:                        {:14f} [s]\n", volumeMoveTail.count());
    std::print(stream, "    Overhead:                    {:14f} [s]\n",
                 volumeMove.count() - volumeMoveNonEwald.count() - volumeMoveEwald.count() - volumeMoveTail.count());
  }

  if(GibbsVolumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs Volume:                    {:14f} [s]\n", GibbsVolumeMove.count());
    std::print(stream, "    Non-Ewald:                   {:14f} [s]\n", GibbsVolumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                       {:14f} [s]\n", GibbsVolumeMoveEwald.count());
    std::print(stream, "    Tail:                        {:14f} [s]\n", GibbsVolumeMoveTail.count());
    std::print(stream, "    Overhead:                    {:14f} [s]\n",
                 GibbsVolumeMove.count() - GibbsVolumeMoveNonEwald.count() - GibbsVolumeMoveEwald.count() - GibbsVolumeMoveTail.count());
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
    std::print(stream, "        ExternalField-Molecule:  {:14f} [s]\n", translationMoveExternalFieldMolecule.count());
    std::print(stream, "        Framework-Molecule:      {:14f} [s]\n", translationMoveFrameworkMolecule.count());
    std::print(stream, "        Molecule-Molecule:       {:14f} [s]\n", translationMoveMoleculeMolecule.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", translationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 translationMove.count() - translationMoveExternalFieldMolecule.count() - translationMoveFrameworkMolecule.count() -
                 translationMoveMoleculeMolecule.count() - translationMoveEwald.count());
  }

  if(randomTranslationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Random Translation:          {:14f} [s]\n", randomTranslationMove.count());
    std::print(stream, "        ExternalField-Molecule:  {:14f} [s]\n", randomTranslationMoveExternalFieldMolecule.count());
    std::print(stream, "        Framework-Molecule:      {:14f} [s]\n", randomTranslationMoveMoleculeMolecule.count());
    std::print(stream, "        Molecule-Molecul:        {:14f} [s]\n", randomTranslationMoveMoleculeMolecule.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", randomTranslationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 randomTranslationMove.count() - randomTranslationMoveExternalFieldMolecule.count() - 
                 randomTranslationMoveFrameworkMolecule.count() - randomTranslationMoveMoleculeMolecule.count() -
                 randomTranslationMoveEwald.count());
  }

  if(rotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Rotation:                    {:14f} [s]\n", rotationMove.count());
    std::print(stream, "        ExternalField-Molecule:  {:14f} [s]\n", rotationMoveExternalFieldMolecule.count());
    std::print(stream, "        Framework-Molecule:      {:14f} [s]\n", rotationMoveMoleculeMolecule.count());
    std::print(stream, "        Molecule-Molecul:        {:14f} [s]\n", rotationMoveMoleculeMolecule.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", rotationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 rotationMove.count() - rotationMoveExternalFieldMolecule.count() - rotationMoveFrameworkMolecule.count() -
                 rotationMoveMoleculeMolecule.count() - rotationMoveEwald.count());
  }

  if(randomRotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Random Rotation:             {:14f} [s]\n", randomRotationMove.count());
    std::print(stream, "        ExternalField-Molecule:  {:14f} [s]\n", randomRotationMoveExternalFieldMolecule.count());
    std::print(stream, "        Framework-Molecule:      {:14f} [s]\n", randomRotationMoveMoleculeMolecule.count());
    std::print(stream, "        Molecule-Molecul:        {:14f} [s]\n", randomRotationMoveMoleculeMolecule.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", randomRotationMoveEwald.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 randomRotationMove.count() - randomRotationMoveExternalFieldMolecule.count() - 
                 randomRotationMoveFrameworkMolecule.count() - randomRotationMoveMoleculeMolecule.count() -
                 randomRotationMoveEwald.count());
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
    std::print(stream, "        Insertion ExternalField: {:14f} [s]\n", swapLambdaInsertionMoveCFCMCExternalField.count());
    std::print(stream, "        Insertion Framework:     {:14f} [s]\n", swapLambdaInsertionMoveCFCMCFramework.count());
    std::print(stream, "        Insertion Molecule:      {:14f} [s]\n", swapLambdaInsertionMoveCFCMCMolecule.count());
    std::print(stream, "        Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCFCMCEwald.count());
    std::print(stream, "        Insertion Tail:          {:14f} [s]\n", swapLambdaInsertionMoveCFCMCTail.count());
    std::print(stream, "        Lambda ExternalField:    {:14f} [s]\n", swapLambdaChangeMoveCFCMCExternalField.count());
    std::print(stream, "        Lambda Framework:        {:14f} [s]\n", swapLambdaChangeMoveCFCMCFramework.count());
    std::print(stream, "        Lambda Molecule:         {:14f} [s]\n", swapLambdaChangeMoveCFCMCMolecule.count());
    std::print(stream, "        Lambda Ewald:            {:14f} [s]\n", swapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "        Lambda Tail:             {:14f} [s]\n", swapLambdaChangeMoveCFCMCTail.count());
    std::print(stream, "        Deletion ExternalField:  {:14f} [s]\n", swapLambdaDeletionMoveCFCMCExternalField.count());
    std::print(stream, "        Deletion Framework:      {:14f} [s]\n", swapLambdaDeletionMoveCFCMCFramework.count());
    std::print(stream, "        Deletion Molecule:    :  {:14f} [s]\n", swapLambdaDeletionMoveCFCMCMolecule.count());
    std::print(stream, "        Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCFCMCEwald.count());
    std::print(stream, "        Deletion Tail:           {:14f} [s]\n", swapLambdaDeletionMoveCFCMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n", swapLambdaMoveCFCMC.count()
                 - swapLambdaInsertionMoveCFCMCExternalField.count() - swapLambdaInsertionMoveCFCMCFramework.count() 
                 - swapLambdaInsertionMoveCFCMCMolecule.count() - swapLambdaInsertionMoveCFCMCEwald.count()
                 - swapLambdaInsertionMoveCFCMCTail.count()
                 - swapLambdaChangeMoveCFCMCExternalField.count() - swapLambdaChangeMoveCFCMCFramework.count()
                 - swapLambdaChangeMoveCFCMCMolecule.count() - swapLambdaChangeMoveCFCMCEwald.count()
                 - swapLambdaChangeMoveCFCMCTail.count()
                 - swapLambdaDeletionMoveCFCMCExternalField.count() - swapLambdaDeletionMoveCFCMCFramework.count()
                 - swapLambdaDeletionMoveCFCMCMolecule.count() - swapLambdaDeletionMoveCFCMCEwald.count()
                 - swapLambdaDeletionMoveCFCMCTail.count());
  }

  if(swapLambdaMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Insertion/deletion CB/CFCMC: {:14f} [s]\n", swapLambdaMoveCBCFCMC.count());
    std::print(stream, "        Insertion ExternalField: {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCExternalField.count());
    std::print(stream, "        Insertion Framework:     {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCFramework.count());
    std::print(stream, "        Insertion Molecule:      {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCMolecule.count());
    std::print(stream, "        Insertion Non-Ewald:     {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCEwald.count());
    std::print(stream, "        Insertion Tail:          {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCTail.count());
    std::print(stream, "        Lambda ExternalField:    {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCExternalField.count());
    std::print(stream, "        Lambda Framework:        {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCFramework.count());
    std::print(stream, "        Lambda Molecule:         {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCMolecule.count());
    std::print(stream, "        Lambda Ewald:            {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCEwald.count());
    std::print(stream, "        Lambda Tail:             {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCTail.count());
    std::print(stream, "        Deletion ExternalField:  {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCExternalField.count());
    std::print(stream, "        Deletion Framework:      {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCFramework.count());
    std::print(stream, "        Deletion Molecule:       {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCMolecule.count());
    std::print(stream, "        Deletion Non-Ewald:      {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCEwald.count());
    std::print(stream, "        Deletion Tail:           {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n", swapLambdaMoveCBCFCMC.count()
                 - swapLambdaInsertionMoveCBCFCMCExternalField.count() - swapLambdaInsertionMoveCBCFCMCFramework.count() - swapLambdaInsertionMoveCBCFCMCMolecule.count()
                 - swapLambdaInsertionMoveCBCFCMCNonEwald.count() - swapLambdaInsertionMoveCBCFCMCEwald.count() - swapLambdaInsertionMoveCBCFCMCTail.count()
                 - swapLambdaChangeMoveCBCFCMCExternalField.count() - swapLambdaChangeMoveCBCFCMCFramework.count() - swapLambdaChangeMoveCBCFCMCMolecule.count()
                 - swapLambdaChangeMoveCBCFCMCEwald.count() - swapLambdaChangeMoveCBCFCMCTail.count()
                 - swapLambdaDeletionMoveCBCFCMCExternalField.count() - swapLambdaDeletionMoveCBCFCMCFramework.count() - swapLambdaDeletionMoveCBCFCMCMolecule.count()
                 - swapLambdaDeletionMoveCBCFCMCNonEwald.count() - swapLambdaDeletionMoveCBCFCMCEwald.count() - swapLambdaDeletionMoveCBCFCMCTail.count());
  }

  if(GibbsSwapMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Gibbs-swap CBMC:             {:14f} [s]\n", GibbsSwapMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", GibbsSwapMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", GibbsSwapMoveCBMCEwald.count());
    std::print(stream, "        Tail:                    {:14f} [s]\n", GibbsSwapMoveCBMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 GibbsSwapMoveCBMC.count() - GibbsSwapMoveCBMCNonEwald.count() - 
                 GibbsSwapMoveCBMCEwald.count() - GibbsSwapMoveCBMCTail.count());
  }

  if(GibbsSwapLambdaMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Gibbs-swap CFCMC:            {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count());
    std::print(stream, "        Inter-change Non-Ewald:  {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count());
    std::print(stream, "        Inter-change Ewald:      {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCEwald.count());
    std::print(stream, "        Inter-change Tail:       {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCTail.count());
    std::print(stream, "        Lambda-change Non-Ewald: {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCNonEwald.count());
    std::print(stream, "        Lambda-change Ewald:     {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "        Lambda-change Tail:      {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCTail.count());
    std::print(stream, "        Shuffle Non-Ewald:       {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count());
    std::print(stream, "        Shuffle Ewald:           {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCEwald.count());
    std::print(stream, "        Shuffle Tail:            {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count()
                 - GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaInterChangeMoveCFCMCEwald.count() -GibbsSwapLambdaInterChangeMoveCFCMCTail.count()
                 - GibbsSwapLambdaChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaChangeMoveCFCMCEwald.count() - GibbsSwapLambdaChangeMoveCFCMCTail.count()
                 - GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count() - GibbsSwapLambdaShuffleMoveCFCMCEwald.count() - GibbsSwapLambdaShuffleMoveCFCMCTail.count());
  }

  if(WidomMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Widom CBMC:                  {:14f} [s]\n", WidomMoveCBMC.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", WidomMoveCBMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", WidomMoveCBMCEwald.count());
    std::print(stream, "        Tail:                    {:14f} [s]\n", WidomMoveCBMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 WidomMoveCBMC.count() - WidomMoveCBMCNonEwald.count() - WidomMoveCBMCEwald.count() - WidomMoveCBMCTail.count());
  }

  if(WidomMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Widom CFCMC:                 {:14f} [s]\n", WidomMoveCFCMC.count());
    std::print(stream, "        ExternalField:           {:14f} [s]\n", WidomMoveCFCMCExternalField.count());
    std::print(stream, "        Framework:               {:14f} [s]\n", WidomMoveCFCMCFramework.count());
    std::print(stream, "        Molecule:                {:14f} [s]\n", WidomMoveCFCMCMolecule.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", WidomMoveCFCMCEwald.count());
    std::print(stream, "        Tail:                    {:14f} [s]\n", WidomMoveCFCMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 WidomMoveCFCMC.count() - WidomMoveCFCMCExternalField.count() 
                 - WidomMoveCFCMCFramework.count() - WidomMoveCFCMCMolecule.count()
                 - WidomMoveCFCMCEwald.count() - WidomMoveCFCMCTail.count());
  }

  if(WidomMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "    Widom CB/CFCMC:              {:14f} [s]\n", WidomMoveCBCFCMC.count());
    std::print(stream, "        ExternalField:           {:14f} [s]\n", WidomMoveCBCFCMCExternalField.count());
    std::print(stream, "        Framework:               {:14f} [s]\n", WidomMoveCBCFCMCFramework.count());
    std::print(stream, "        Molecule:                {:14f} [s]\n", WidomMoveCBCFCMCMolecule.count());
    std::print(stream, "        Non-Ewald:               {:14f} [s]\n", WidomMoveCBCFCMCNonEwald.count());
    std::print(stream, "        Ewald:                   {:14f} [s]\n", WidomMoveCBCFCMCEwald.count());
    std::print(stream, "        Tail:                    {:14f} [s]\n", WidomMoveCBCFCMCTail.count());
    std::print(stream, "        Overhead:                {:14f} [s]\n",
                 WidomMoveCBCFCMC.count() - WidomMoveCBCFCMCExternalField.count() - WidomMoveCBCFCMCFramework.count() - WidomMoveCBCFCMCMolecule.count()
                 - WidomMoveCBCFCMCNonEwald.count() - WidomMoveCBCFCMCEwald.count() - WidomMoveCBCFCMCTail.count());
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
    std::print(stream, "Translation:                   {:14f} [s]\n", translationMove.count());
    std::print(stream, "    ExternalField-Molecule:    {:14f} [s]\n", translationMoveExternalFieldMolecule.count());
    std::print(stream, "    Framework-Molecule:        {:14f} [s]\n", translationMoveFrameworkMolecule.count());
    std::print(stream, "    Molecule-Molecule:         {:14f} [s]\n", translationMoveMoleculeMolecule.count());
    std::print(stream, "    Ewald:                     {:14f} [s]\n", translationMoveEwald.count());
    std::print(stream, "    Overhead:                  {:14f} [s]\n",
                 translationMove.count() - translationMoveExternalFieldMolecule.count() - translationMoveFrameworkMolecule.count() -
                 translationMoveMoleculeMolecule.count() - translationMoveEwald.count());
  }

  if(randomTranslationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Random Translation:          {:14f} [s]\n", randomTranslationMove.count());
    std::print(stream, "    ExternalField-Molecule:  {:14f} [s]\n", randomTranslationMoveExternalFieldMolecule.count());
    std::print(stream, "    Framework-Molecule:      {:14f} [s]\n", randomTranslationMoveMoleculeMolecule.count());
    std::print(stream, "    Molecule-Molecul:        {:14f} [s]\n", randomTranslationMoveMoleculeMolecule.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", randomTranslationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 randomTranslationMove.count() - randomTranslationMoveExternalFieldMolecule.count() - 
                 randomTranslationMoveFrameworkMolecule.count() - randomTranslationMoveMoleculeMolecule.count() -
                 randomTranslationMoveEwald.count());
  }

  if(rotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Rotation:                    {:14f} [s]\n", rotationMove.count());
    std::print(stream, "    ExternalField-Molecule:  {:14f} [s]\n", rotationMoveExternalFieldMolecule.count());
    std::print(stream, "    Framework-Molecule:      {:14f} [s]\n", rotationMoveMoleculeMolecule.count());
    std::print(stream, "    Molecule-Molecul:        {:14f} [s]\n", rotationMoveMoleculeMolecule.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", rotationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 rotationMove.count() - rotationMoveExternalFieldMolecule.count() - rotationMoveFrameworkMolecule.count() -
                 rotationMoveMoleculeMolecule.count() - rotationMoveEwald.count());
  }

  if(randomRotationMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Random rotation:             {:14f} [s]\n", randomRotationMove.count());
    std::print(stream, "    ExternalField-Molecule:  {:14f} [s]\n", randomRotationMoveExternalFieldMolecule.count());
    std::print(stream, "    Framework-Molecule:      {:14f} [s]\n", randomRotationMoveMoleculeMolecule.count());
    std::print(stream, "    Molecule-Molecul:        {:14f} [s]\n", randomRotationMoveMoleculeMolecule.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", randomRotationMoveEwald.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 randomRotationMove.count() - randomRotationMoveExternalFieldMolecule.count() -
                 randomRotationMoveFrameworkMolecule.count() - randomRotationMoveMoleculeMolecule.count() -
                 randomRotationMoveEwald.count());
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
    std::print(stream, "    Insertion ExternalField: {:14f} [s]\n", swapLambdaInsertionMoveCFCMCExternalField.count());
    std::print(stream, "    Insertion Framework:     {:14f} [s]\n", swapLambdaInsertionMoveCFCMCFramework.count());
    std::print(stream, "    Insertion Molecule:      {:14f} [s]\n", swapLambdaInsertionMoveCFCMCMolecule.count());
    std::print(stream, "    Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCFCMCEwald.count());
    std::print(stream, "    Insertion Tail:          {:14f} [s]\n", swapLambdaInsertionMoveCFCMCTail.count());
    std::print(stream, "    Lambda ExternalField:    {:14f} [s]\n", swapLambdaChangeMoveCFCMCExternalField.count());
    std::print(stream, "    Lambda Framework:        {:14f} [s]\n", swapLambdaChangeMoveCFCMCFramework.count());
    std::print(stream, "    Lambda Molecule:         {:14f} [s]\n", swapLambdaChangeMoveCFCMCMolecule.count());
    std::print(stream, "    Lambda Ewald:            {:14f} [s]\n", swapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "    Lambda Tail:             {:14f} [s]\n", swapLambdaChangeMoveCFCMCTail.count());
    std::print(stream, "    Deletion ExternalField:  {:14f} [s]\n", swapLambdaDeletionMoveCFCMCExternalField.count());
    std::print(stream, "    Deletion Framework:      {:14f} [s]\n", swapLambdaDeletionMoveCFCMCFramework.count());
    std::print(stream, "    Deletion Molecule:    :  {:14f} [s]\n", swapLambdaDeletionMoveCFCMCMolecule.count());
    std::print(stream, "    Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCFCMCEwald.count());
    std::print(stream, "    Deletion Tail:           {:14f} [s]\n", swapLambdaDeletionMoveCFCMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n", swapLambdaMoveCFCMC.count()
                 - swapLambdaInsertionMoveCFCMCExternalField.count() - swapLambdaInsertionMoveCFCMCFramework.count()
                 - swapLambdaInsertionMoveCFCMCMolecule.count() - swapLambdaInsertionMoveCFCMCEwald.count()
                 - swapLambdaInsertionMoveCFCMCTail.count()
                 - swapLambdaChangeMoveCFCMCExternalField.count() - swapLambdaChangeMoveCFCMCFramework.count()
                 - swapLambdaChangeMoveCFCMCMolecule.count() - swapLambdaChangeMoveCFCMCEwald.count()
                 - swapLambdaChangeMoveCFCMCTail.count()
                 - swapLambdaDeletionMoveCFCMCExternalField.count() - swapLambdaDeletionMoveCFCMCFramework.count()
                 - swapLambdaDeletionMoveCFCMCMolecule.count() - swapLambdaDeletionMoveCFCMCEwald.count()
                 - swapLambdaDeletionMoveCFCMCTail.count());
  }

  if(swapLambdaMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Insertion/deletion CB/CFCMC: {:14f} [s]\n", swapLambdaMoveCBCFCMC.count());
    std::print(stream, "    Insertion ExternalField: {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCExternalField.count());
    std::print(stream, "    Insertion Framework:     {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCFramework.count());
    std::print(stream, "    Insertion Molecule:      {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCMolecule.count());
    std::print(stream, "    Insertion Non-Ewald:     {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Insertion Ewald:         {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCEwald.count());
    std::print(stream, "    Insertion Tail:          {:14f} [s]\n", swapLambdaInsertionMoveCBCFCMCTail.count());
    std::print(stream, "    Lambda ExternalField:    {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCExternalField.count());
    std::print(stream, "    Lambda Framework:        {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCFramework.count());
    std::print(stream, "    Lambda Molecule:         {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCMolecule.count());
    std::print(stream, "    Lambda Ewald:            {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCEwald.count());
    std::print(stream, "    Lambda Tail:             {:14f} [s]\n", swapLambdaChangeMoveCBCFCMCTail.count());
    std::print(stream, "    Deletion ExternalField:  {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCExternalField.count());
    std::print(stream, "    Deletion Framework:      {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCFramework.count());
    std::print(stream, "    Deletion Molecule:       {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCMolecule.count());
    std::print(stream, "    Deletion Non-Ewald:      {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Deletion Ewald:          {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCEwald.count());
    std::print(stream, "    Deletion Tail:           {:14f} [s]\n", swapLambdaDeletionMoveCBCFCMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n", swapLambdaMoveCBCFCMC.count()
                 - swapLambdaInsertionMoveCBCFCMCExternalField.count() - swapLambdaInsertionMoveCBCFCMCFramework.count() - swapLambdaInsertionMoveCBCFCMCMolecule.count()
                 - swapLambdaInsertionMoveCBCFCMCNonEwald.count() - swapLambdaInsertionMoveCBCFCMCEwald.count() - swapLambdaInsertionMoveCBCFCMCTail.count()
                 - swapLambdaChangeMoveCBCFCMCExternalField.count() - swapLambdaChangeMoveCBCFCMCFramework.count() - swapLambdaChangeMoveCBCFCMCMolecule.count()
                 - swapLambdaChangeMoveCBCFCMCEwald.count() - swapLambdaChangeMoveCBCFCMCTail.count()
                 - swapLambdaDeletionMoveCBCFCMCExternalField.count() - swapLambdaDeletionMoveCBCFCMCFramework.count() - swapLambdaDeletionMoveCBCFCMCMolecule.count()
                 - swapLambdaDeletionMoveCBCFCMCNonEwald.count() - swapLambdaDeletionMoveCBCFCMCEwald.count() - swapLambdaDeletionMoveCBCFCMCTail.count());
  }

  if(GibbsSwapMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs-swap CBMC:             {:14f} [s]\n", GibbsSwapMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", GibbsSwapMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", GibbsSwapMoveCBMCEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", GibbsSwapMoveCBMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 GibbsSwapMoveCBMC.count() - GibbsSwapMoveCBMCNonEwald.count() - 
                 GibbsSwapMoveCBMCEwald.count() - GibbsSwapMoveCBMCTail.count());
  }

  if(GibbsSwapLambdaMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs-swap CFCMC:            {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count());
    std::print(stream, "    Inter-change Non-Ewald:  {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count());
    std::print(stream, "    Inter-change Ewald:      {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCEwald.count());
    std::print(stream, "    Inter-change Tail:       {:14f} [s]\n", GibbsSwapLambdaInterChangeMoveCFCMCTail.count());
    std::print(stream, "    Lambda-change Non-Ewald: {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCNonEwald.count());
    std::print(stream, "    Lambda-change Ewald:     {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCEwald.count());
    std::print(stream, "    Lambda-change Tail:      {:14f} [s]\n", GibbsSwapLambdaChangeMoveCFCMCTail.count());
    std::print(stream, "    Shuffle Non-Ewald:       {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count());
    std::print(stream, "    Shuffle Ewald:           {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCEwald.count());
    std::print(stream, "    Shuffle Tail:            {:14f} [s]\n", GibbsSwapLambdaShuffleMoveCFCMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n", GibbsSwapLambdaMoveCFCMC.count()
                 - GibbsSwapLambdaInterChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaInterChangeMoveCFCMCEwald.count() - GibbsSwapLambdaInterChangeMoveCFCMCTail.count()
                 - GibbsSwapLambdaChangeMoveCFCMCNonEwald.count() - GibbsSwapLambdaChangeMoveCFCMCEwald.count() - GibbsSwapLambdaChangeMoveCFCMCTail.count()
                 - GibbsSwapLambdaShuffleMoveCFCMCNonEwald.count() - GibbsSwapLambdaShuffleMoveCFCMCEwald.count() - GibbsSwapLambdaShuffleMoveCFCMCTail.count());
  }

  if(volumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Volume:                      {:14f} [s]\n", volumeMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", volumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", volumeMoveEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", volumeMoveTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 volumeMove.count() - volumeMoveNonEwald.count() - volumeMoveEwald.count() - volumeMoveTail.count());
  }

  if(WidomMoveCBMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Widom CBMC:                  {:14f} [s]\n", WidomMoveCBMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", WidomMoveCBMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", WidomMoveCBMCEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", WidomMoveCBMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 WidomMoveCBMC.count() - WidomMoveCBMCNonEwald.count() - WidomMoveCBMCEwald.count() - WidomMoveCBMCTail.count());
  }

  if(WidomMoveCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Widom CFCMC:                 {:14f} [s]\n", WidomMoveCFCMC.count());
    std::print(stream, "    ExternalField:           {:14f} [s]\n", WidomMoveCFCMCExternalField.count());
    std::print(stream, "    Framework:               {:14f} [s]\n", WidomMoveCFCMCFramework.count());
    std::print(stream, "    Molecule:                {:14f} [s]\n", WidomMoveCFCMCMolecule.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", WidomMoveCFCMCEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", WidomMoveCFCMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 WidomMoveCFCMC.count() - WidomMoveCFCMCExternalField.count()
                 - WidomMoveCFCMCFramework.count() - WidomMoveCFCMCMolecule.count()
                 - WidomMoveCFCMCEwald.count() - WidomMoveCFCMCTail.count());
  }

  if(WidomMoveCBCFCMC > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Widom CB/CFCMC:              {:14f} [s]\n", WidomMoveCBCFCMC.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", WidomMoveCBCFCMCNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", WidomMoveCBCFCMCEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", WidomMoveCBCFCMCTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 WidomMoveCBCFCMC.count() - WidomMoveCBCFCMCNonEwald.count() - WidomMoveCBCFCMCEwald.count() - WidomMoveCBCFCMCTail.count());
  }

  if(GibbsVolumeMove > std::chrono::duration<double>::zero())
  {
    std::print(stream, "\n");
    std::print(stream, "Gibbs Volume:                {:14f} [s]\n", GibbsVolumeMove.count());
    std::print(stream, "    Non-Ewald:               {:14f} [s]\n", GibbsVolumeMoveNonEwald.count());
    std::print(stream, "    Ewald:                   {:14f} [s]\n", GibbsVolumeMoveEwald.count());
    std::print(stream, "    Tail:                    {:14f} [s]\n", GibbsVolumeMoveTail.count());
    std::print(stream, "    Overhead:                {:14f} [s]\n",
                 GibbsVolumeMove.count() - GibbsVolumeMoveNonEwald.count() - GibbsVolumeMoveEwald.count() - GibbsVolumeMoveTail.count());
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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveCpuTime &t)
{
  archive << t.versionNumber;

  archive << t.propertySampling;
  archive << t.energyPressureComputation;

  archive << t.translationMove;
  archive << t.translationMoveExternalFieldMolecule;
  archive << t.translationMoveFrameworkMolecule;
  archive << t.translationMoveMoleculeMolecule;
  archive << t.translationMoveEwald;

  archive << t.randomTranslationMove;
  archive << t.randomTranslationMoveExternalFieldMolecule;
  archive << t.randomTranslationMoveFrameworkMolecule;
  archive << t.randomTranslationMoveMoleculeMolecule;
  archive << t.randomTranslationMoveEwald;

  archive << t.rotationMove;
  archive << t.rotationMoveExternalFieldMolecule;
  archive << t.rotationMoveFrameworkMolecule;
  archive << t.rotationMoveMoleculeMolecule;
  archive << t.rotationMoveEwald;

  archive << t.randomRotationMove;
  archive << t.randomRotationMoveExternalFieldMolecule;
  archive << t.randomRotationMoveFrameworkMolecule;
  archive << t.randomRotationMoveMoleculeMolecule;
  archive << t.randomRotationMoveEwald;

  archive << t.reinsertionMoveCBMC;
  archive << t.reinsertionMoveCBMCNonEwald;
  archive << t.reinsertionMoveCBMCEwald;

  archive << t.swapInsertionMoveCBMC;
  archive << t.swapInsertionMoveCBMCNonEwald;
  archive << t.swapInsertionMoveCBMCEwald;
  archive << t.swapInsertionMoveCBMCTail;

  archive << t.swapDeletionMoveCBMC;
  archive << t.swapDeletionMoveCBMCNonEwald;
  archive << t.swapDeletionMoveCBMCEwald;
  archive << t.swapDeletionMoveCBMCTail;

  archive << t.swapLambdaDeletionMove;
  archive << t.swapLambdaDeletionMoveNonEwald;
  archive << t.swapLambdaDeletionMoveEwald;

  archive << t.swapLambdaMoveCFCMC;
  archive << t.swapLambdaInsertionMoveCFCMCExternalField;
  archive << t.swapLambdaInsertionMoveCFCMCFramework;
  archive << t.swapLambdaInsertionMoveCFCMCMolecule;
  archive << t.swapLambdaInsertionMoveCFCMCEwald;
  archive << t.swapLambdaInsertionMoveCFCMCTail;
  archive << t.swapLambdaChangeMoveCFCMCExternalField;
  archive << t.swapLambdaChangeMoveCFCMCFramework;
  archive << t.swapLambdaChangeMoveCFCMCMolecule;
  archive << t.swapLambdaChangeMoveCFCMCEwald;
  archive << t.swapLambdaChangeMoveCFCMCTail;
  archive << t.swapLambdaDeletionMoveCFCMCExternalField;
  archive << t.swapLambdaDeletionMoveCFCMCFramework;
  archive << t.swapLambdaDeletionMoveCFCMCMolecule;
  archive << t.swapLambdaDeletionMoveCFCMCEwald;
  archive << t.swapLambdaDeletionMoveCFCMCTail;

  archive << t.swapLambdaMoveCBCFCMC;
  archive << t.swapLambdaInsertionMoveCBCFCMCExternalField;
  archive << t.swapLambdaInsertionMoveCBCFCMCFramework;
  archive << t.swapLambdaInsertionMoveCBCFCMCMolecule;
  archive << t.swapLambdaInsertionMoveCBCFCMCNonEwald;
  archive << t.swapLambdaInsertionMoveCBCFCMCEwald;
  archive << t.swapLambdaInsertionMoveCBCFCMCTail;
  archive << t.swapLambdaChangeMoveCBCFCMCExternalField;
  archive << t.swapLambdaChangeMoveCBCFCMCFramework;
  archive << t.swapLambdaChangeMoveCBCFCMCMolecule;
  archive << t.swapLambdaChangeMoveCBCFCMCEwald;
  archive << t.swapLambdaChangeMoveCBCFCMCTail;
  archive << t.swapLambdaDeletionMoveCBCFCMCExternalField;
  archive << t.swapLambdaDeletionMoveCBCFCMCFramework;
  archive << t.swapLambdaDeletionMoveCBCFCMCMolecule;
  archive << t.swapLambdaDeletionMoveCBCFCMCNonEwald;
  archive << t.swapLambdaDeletionMoveCBCFCMCEwald;
  archive << t.swapLambdaDeletionMoveCBCFCMCTail;

  archive << t.GibbsSwapMoveCBMC;
  archive << t.GibbsSwapMoveCBMCNonEwald;
  archive << t.GibbsSwapMoveCBMCEwald;
  archive << t.GibbsSwapMoveCBMCTail;

  archive << t.GibbsSwapLambdaMoveCFCMC;
  archive << t.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald;
  archive << t.GibbsSwapLambdaInterChangeMoveCFCMCEwald;
  archive << t.GibbsSwapLambdaInterChangeMoveCFCMCTail;
  archive << t.GibbsSwapLambdaChangeMoveCFCMCNonEwald;
  archive << t.GibbsSwapLambdaChangeMoveCFCMCEwald;
  archive << t.GibbsSwapLambdaChangeMoveCFCMCTail;
  archive << t.GibbsSwapLambdaShuffleMoveCFCMCNonEwald;
  archive << t.GibbsSwapLambdaShuffleMoveCFCMCEwald;
  archive << t.GibbsSwapLambdaShuffleMoveCFCMCTail;

  archive << t.WidomMoveCBMC;
  archive << t.WidomMoveCBMCNonEwald;
  archive << t.WidomMoveCBMCEwald;
  archive << t.WidomMoveCBMCTail;

  archive << t.WidomMoveCFCMC;
  archive << t.WidomMoveCFCMCExternalField;
  archive << t.WidomMoveCFCMCFramework;
  archive << t.WidomMoveCFCMCMolecule;
  archive << t.WidomMoveCFCMCEwald;
  archive << t.WidomMoveCFCMCTail;

  archive << t.WidomMoveCBCFCMC;
  archive << t.WidomMoveCBCFCMCExternalField;
  archive << t.WidomMoveCBCFCMCFramework;
  archive << t.WidomMoveCBCFCMCMolecule;
  archive << t.WidomMoveCBCFCMCNonEwald;
  archive << t.WidomMoveCBCFCMCEwald;
  archive << t.WidomMoveCBCFCMCTail;

  archive << t.volumeMove;
  archive << t.volumeMoveNonEwald;
  archive << t.volumeMoveEwald;
  archive << t.volumeMoveTail;

  archive << t.GibbsVolumeMove;
  archive << t.GibbsVolumeMoveNonEwald;
  archive << t.GibbsVolumeMoveEwald;
  archive << t.GibbsVolumeMoveTail;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveCpuTime &t)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > t.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveCpuTime' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> t.propertySampling;
  archive >> t.energyPressureComputation;

  archive >> t.translationMove;
  archive >> t.translationMoveExternalFieldMolecule;
  archive >> t.translationMoveFrameworkMolecule;
  archive >> t.translationMoveMoleculeMolecule;
  archive >> t.translationMoveEwald;

  archive >> t.randomTranslationMove;
  archive >> t.randomTranslationMoveExternalFieldMolecule;
  archive >> t.randomTranslationMoveFrameworkMolecule;
  archive >> t.randomTranslationMoveMoleculeMolecule;
  archive >> t.randomTranslationMoveEwald;

  archive >> t.rotationMove;
  archive >> t.rotationMoveExternalFieldMolecule;
  archive >> t.rotationMoveFrameworkMolecule;
  archive >> t.rotationMoveMoleculeMolecule;
  archive >> t.rotationMoveEwald;

  archive >> t.randomRotationMove;
  archive >> t.randomRotationMoveExternalFieldMolecule;
  archive >> t.randomRotationMoveFrameworkMolecule;
  archive >> t.randomRotationMoveMoleculeMolecule;
  archive >> t.randomRotationMoveEwald;

  archive >> t.reinsertionMoveCBMC;
  archive >> t.reinsertionMoveCBMCNonEwald;
  archive >> t.reinsertionMoveCBMCEwald;

  archive >> t.swapInsertionMoveCBMC;
  archive >> t.swapInsertionMoveCBMCNonEwald;
  archive >> t.swapInsertionMoveCBMCEwald;
  archive >> t.swapInsertionMoveCBMCTail;

  archive >> t.swapDeletionMoveCBMC;
  archive >> t.swapDeletionMoveCBMCNonEwald;
  archive >> t.swapDeletionMoveCBMCEwald;
  archive >> t.swapDeletionMoveCBMCTail;

  archive >> t.swapLambdaDeletionMove;
  archive >> t.swapLambdaDeletionMoveNonEwald;
  archive >> t.swapLambdaDeletionMoveEwald;

  archive >> t.swapLambdaMoveCFCMC;
  archive >> t.swapLambdaInsertionMoveCFCMCExternalField;
  archive >> t.swapLambdaInsertionMoveCFCMCFramework;
  archive >> t.swapLambdaInsertionMoveCFCMCMolecule;
  archive >> t.swapLambdaInsertionMoveCFCMCEwald;
  archive >> t.swapLambdaInsertionMoveCFCMCTail;
  archive >> t.swapLambdaChangeMoveCFCMCExternalField;
  archive >> t.swapLambdaChangeMoveCFCMCFramework;
  archive >> t.swapLambdaChangeMoveCFCMCMolecule;
  archive >> t.swapLambdaChangeMoveCFCMCEwald;
  archive >> t.swapLambdaChangeMoveCFCMCTail;
  archive >> t.swapLambdaDeletionMoveCFCMCExternalField;
  archive >> t.swapLambdaDeletionMoveCFCMCFramework;
  archive >> t.swapLambdaDeletionMoveCFCMCMolecule;
  archive >> t.swapLambdaDeletionMoveCFCMCEwald;
  archive >> t.swapLambdaDeletionMoveCFCMCTail;

  archive >> t.swapLambdaMoveCBCFCMC;
  archive >> t.swapLambdaInsertionMoveCBCFCMCExternalField;
  archive >> t.swapLambdaInsertionMoveCBCFCMCFramework;
  archive >> t.swapLambdaInsertionMoveCBCFCMCMolecule;
  archive >> t.swapLambdaInsertionMoveCBCFCMCNonEwald;
  archive >> t.swapLambdaInsertionMoveCBCFCMCEwald;
  archive >> t.swapLambdaInsertionMoveCBCFCMCTail;
  archive >> t.swapLambdaChangeMoveCBCFCMCExternalField;
  archive >> t.swapLambdaChangeMoveCBCFCMCFramework;
  archive >> t.swapLambdaChangeMoveCBCFCMCMolecule;
  archive >> t.swapLambdaChangeMoveCBCFCMCEwald;
  archive >> t.swapLambdaChangeMoveCBCFCMCTail;
  archive >> t.swapLambdaDeletionMoveCBCFCMCExternalField;
  archive >> t.swapLambdaDeletionMoveCBCFCMCFramework;
  archive >> t.swapLambdaDeletionMoveCBCFCMCMolecule;
  archive >> t.swapLambdaDeletionMoveCBCFCMCNonEwald;
  archive >> t.swapLambdaDeletionMoveCBCFCMCEwald;
  archive >> t.swapLambdaDeletionMoveCBCFCMCTail;

  archive >> t.GibbsSwapMoveCBMC;
  archive >> t.GibbsSwapMoveCBMCNonEwald;
  archive >> t.GibbsSwapMoveCBMCEwald;
  archive >> t.GibbsSwapMoveCBMCTail;

  archive >> t.GibbsSwapLambdaMoveCFCMC;
  archive >> t.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald;
  archive >> t.GibbsSwapLambdaInterChangeMoveCFCMCEwald;
  archive >> t.GibbsSwapLambdaInterChangeMoveCFCMCTail;
  archive >> t.GibbsSwapLambdaChangeMoveCFCMCNonEwald;
  archive >> t.GibbsSwapLambdaChangeMoveCFCMCEwald;
  archive >> t.GibbsSwapLambdaChangeMoveCFCMCTail;
  archive >> t.GibbsSwapLambdaShuffleMoveCFCMCNonEwald;
  archive >> t.GibbsSwapLambdaShuffleMoveCFCMCEwald;
  archive >> t.GibbsSwapLambdaShuffleMoveCFCMCTail;

  archive >> t.WidomMoveCBMC;
  archive >> t.WidomMoveCBMCNonEwald;
  archive >> t.WidomMoveCBMCEwald;
  archive >> t.WidomMoveCBMCTail;

  archive >> t.WidomMoveCFCMC;
  archive >> t.WidomMoveCFCMCExternalField;
  archive >> t.WidomMoveCFCMCFramework;
  archive >> t.WidomMoveCFCMCMolecule;
  archive >> t.WidomMoveCFCMCEwald;
  archive >> t.WidomMoveCFCMCTail;

  archive >> t.WidomMoveCBCFCMC;
  archive >> t.WidomMoveCBCFCMCExternalField;
  archive >> t.WidomMoveCBCFCMCFramework;
  archive >> t.WidomMoveCBCFCMCMolecule;
  archive >> t.WidomMoveCBCFCMCNonEwald;
  archive >> t.WidomMoveCBCFCMCEwald;
  archive >> t.WidomMoveCBCFCMCTail;

  archive >> t.volumeMove;
  archive >> t.volumeMoveNonEwald;
  archive >> t.volumeMoveEwald;
  archive >> t.volumeMoveTail;

  archive >> t.GibbsVolumeMove;
  archive >> t.GibbsVolumeMoveNonEwald;
  archive >> t.GibbsVolumeMoveEwald;
  archive >> t.GibbsVolumeMoveTail;

  return archive;
}

