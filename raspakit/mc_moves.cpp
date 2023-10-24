module;

module mc_moves;

import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
import <fstream>;
import <format>;
import <print>;
import <exception>;
import <source_location>;
import <complex>;

import archive;
import double3;
import double3x3;
import simd_quatd;
import randomnumbers;

import component;
import atom;
import simulationbox;
import cbmc;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import mc_moves_probabilities_particles;
import mc_moves_probabilities_system;
import mc_moves_cputime;
import transition_matrix;


void MC_Moves::performRandomMove(RandomNumber &random, System& selectedSystem, System& selectedSecondSystem, size_t selectedComponent, size_t &fractionalMoleculeSystem)
{
  double randomNumber = random.uniform();

  MCMoveProbabilitiesParticles &mc_moves_probabilities = selectedSystem.components[selectedComponent].mc_moves_probabilities;

  size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  if (randomNumber < mc_moves_probabilities.accumulatedProbabilityTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = translationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityRandomTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomTranslationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = rotationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityRandomRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomRotationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityVolumeMove)
  {
    std::optional<RunningEnergy> energy = volumeMove(random, selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityReinsertionMove_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = reinsertionMove(random, selectedSystem, selectedComponent, selectedMolecule, molecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }

      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityIdentityChangeMove_CBMC)
  {
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilitySwapMove_CBMC)
  {
    if (random.uniform() < 0.5)
    {
      const auto [energyDifference, Pacc] = insertionMove(random, selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

      const auto [energyDifference, Pacc] = deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    const auto [energyDifference, Pacc] = swapMove_CFCMC(random, selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    selectedSystem.tmmc.updateMatrix(Pacc, oldN);
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    const auto [energyDifference, Pacc] = swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    selectedSystem.tmmc.updateMatrix(Pacc, oldN);
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityGibbsVolumeMove)
  {
    std::optional<std::pair<RunningEnergy,RunningEnergy>> energy = GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CBMC)
  {
    if(random.uniform() < 0.5)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(random, selectedSystem, selectedSecondSystem, fractionalMoleculeSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else if(selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(random, selectedSecondSystem, selectedSystem, fractionalMoleculeSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityWidomMove)
  {
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC)
  {
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC_CBMC)
  {
    const auto [energyDifference, Pacc] = swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
  }
}

void MC_Moves::performRandomMoveProduction(RandomNumber &random, System& selectedSystem, System& selectedSecondSystem, size_t selectedComponent, size_t &fractionalMoleculeSystem, size_t currentBlock)
{
  double randomNumber = random.uniform();

  MCMoveProbabilitiesParticles &mc_moves_probabilities = selectedSystem.components[selectedComponent].mc_moves_probabilities;
  size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  if (randomNumber < mc_moves_probabilities.accumulatedProbabilityTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = translationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.translationMove += (t2 - t1);
    selectedSystem.mc_moves_cputime.translationMove += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.translationMove++;
    selectedSystem.mc_moves_count.translationMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityRandomTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomTranslationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.randomTranslationMove += (t2 - t1);
    selectedSystem.mc_moves_cputime.randomTranslationMove += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.randomTranslationMove++;
    selectedSystem.mc_moves_count.randomTranslationMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = rotationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.rotationMove += (t2 - t1);
    selectedSystem.mc_moves_cputime.rotationMove += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.rotationMove++;
    selectedSystem.mc_moves_count.rotationMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityRandomRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomRotationMove(random, selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.randomRotationMove += (t2 - t1);
    selectedSystem.mc_moves_cputime.randomRotationMove += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.randomRotationMove++;
    selectedSystem.mc_moves_count.randomRotationMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityVolumeMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> energy = volumeMove(random, selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.mc_moves_cputime.volumeMove += (t2 - t1);
    selectedSystem.mc_moves_count.volumeMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityReinsertionMove_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = reinsertionMove(random, selectedSystem, selectedComponent, selectedMolecule, molecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.reinsertionMoveCBMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.reinsertionMoveCBMC++;
    selectedSystem.mc_moves_count.reinsertionMoveCBMC++;
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityIdentityChangeMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CBMC)
  {
    if (random.uniform() < 0.5)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] = insertionMove(random, selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.swapInsertionMoveCBMC += (t2 - t1);

      selectedSystem.components[selectedComponent].mc_moves_count.swapInsertionMoveCBMC++;
      selectedSystem.mc_moves_count.swapInsertionMoveCBMC++;
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] = deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.swapDeletionMoveCBMC += (t2 - t1);

      selectedSystem.components[selectedComponent].mc_moves_count.swapDeletionMoveCBMC++;
      selectedSystem.mc_moves_count.swapDeletionMoveCBMC++;
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] = swapMove_CFCMC(random, selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    selectedSystem.tmmc.updateMatrix(Pacc, oldN);

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.swapLambdaMoveCFCMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.swapLambdaMoveCFCMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.swapLambdaMoveCFCMC++;
    selectedSystem.mc_moves_count.swapLambdaMoveCFCMC++;
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] = swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    selectedSystem.tmmc.updateMatrix(Pacc, oldN);

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.swapLambdaMoveCBCFCMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.swapLambdaMoveCBCFCMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.swapLambdaMoveCBCFCMC++;
    selectedSystem.mc_moves_count.swapLambdaMoveCBCFCMC++;
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsVolumeMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.mc_moves_cputime.GibbsVolumeMove += (t2 - t1);
    selectedSystem.mc_moves_count.GibbsVolumeMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CBMC)
  {
    if (random.uniform() < 0.5)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);

      selectedSystem.components[selectedComponent].mc_moves_count.GibbsSwapMoveCBMC++;
      selectedSystem.mc_moves_count.GibbsSwapMoveCBMC++;
    }
    else
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSecondSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);
      selectedSecondSystem.mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);

      selectedSecondSystem.components[selectedComponent].mc_moves_count.GibbsSwapMoveCBMC++;
      selectedSecondSystem.mc_moves_count.GibbsSwapMoveCBMC++;
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);

      selectedSystem.components[selectedComponent].mc_moves_count.GibbsSwapLambdaMoveCFCMC++;
      selectedSystem.mc_moves_count.GibbsSwapLambdaMoveCFCMC++;
    }
    else if (selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSecondSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);
      selectedSecondSystem.mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);

      selectedSecondSystem.components[selectedComponent].mc_moves_count.GibbsSwapLambdaMoveCFCMC++;
      selectedSecondSystem.mc_moves_count.GibbsSwapLambdaMoveCFCMC++;
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityWidomMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::optional<double> RosenbluthWeight = WidomMove(random, selectedSystem, selectedComponent);

    selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, RosenbluthWeight.value_or(0.0), selectedSystem.weight());

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCBMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.WidomMoveCBMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.WidomMoveCBMC++;
    selectedSystem.mc_moves_count.WidomMoveCBMC++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] = swapMove_CFCMC(random, selectedSystem, selectedComponent, 0, true, true);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCFCMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.WidomMoveCFCMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.WidomMoveCFCMC++;
    selectedSystem.mc_moves_count.WidomMoveCFCMC++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC_CBMC)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] = swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCBCFCMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.WidomMoveCBCFCMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.WidomMoveCBCFCMC++;
    selectedSystem.mc_moves_count.WidomMoveCBCFCMC++;
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MC_Moves &mc)
{
  archive << mc.versionNumber;

  archive << mc.energyOverlapCriteria;
  archive << mc.useDualCutOff;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MC_Moves &mc)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MC_Moves' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> mc.energyOverlapCriteria;
  archive >> mc.useDualCutOff;

  return archive;
}
