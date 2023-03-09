module;

module mc_moves;

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import lambda;
import property_widom;
import averages;

import <complex>;
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


void MC_Moves::performRandomMove(System& selectedSystem, size_t selectedComponent)
{
  double randomNumber = RandomNumber::Uniform();

  if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = translationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityRandomTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomTranslationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = rotationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityRandomRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomRotationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityVolumeMove)
  {
    std::optional<RunningEnergy> energy = volumeMove(selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityReinsertionMove_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = reinsertionMove(selectedSystem, selectedComponent, selectedMolecule, molecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityIdentityChangeMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CBMC)
  {
    if (RandomNumber::Uniform() < 0.5)
    {
      std::optional<RunningEnergy> energyDifference = insertionMove(selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(selectedComponent);

      std::optional<RunningEnergy> energyDifference = deletionMove(selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    std::optional<RunningEnergy> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsVolumeMove)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC_CBMC)
  {
    std::optional<RunningEnergy> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, 0, true, true);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
  }
}

void MC_Moves::performRandomMoveProduction(System& selectedSystem, size_t selectedComponent, size_t currentBlock)
{
  double randomNumber = RandomNumber::Uniform();

  if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = translationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_TranslationMove += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityRandomTranslationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomTranslationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_RandomTranslationMove += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = rotationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_RotationMove += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityRandomRotationMove)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = randomRotationMove(selectedSystem, selectedComponent, molecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_RandomRotationMove += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityVolumeMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> energy = volumeMove(selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.cpuTime_VolumeMove += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityReinsertionMove_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      std::optional<RunningEnergy> energyDifference = reinsertionMove(selectedSystem, selectedComponent, selectedMolecule, molecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_ReinsertionMove_CBMC += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityIdentityChangeMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CBMC)
  {
    if (RandomNumber::Uniform() < 0.5)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      std::optional<RunningEnergy> energyDifference = insertionMove(selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_SwapInsertionMove_CBMC += (t2 - t1);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      std::optional<RunningEnergy> energyDifference = deletionMove(selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_SwapDeletionMove_CBMC += (t2 - t1);
    }
      
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC_CBMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::optional<RunningEnergy> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_SwapMove_CFCMC_CBMC += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsVolumeMove)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::optional<double> RosenbluthWeight = WidomMove(selectedSystem, selectedComponent);
    selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, RosenbluthWeight.value_or(0.0), selectedSystem.weight());

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_WidomMove_CBMC += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove_CFCMC_CBMC)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::optional<RunningEnergy> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, 0, true, true);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_timings.cpuTime_WidomMove_CFCMC_CBMC += (t2 - t1);
  }
}

