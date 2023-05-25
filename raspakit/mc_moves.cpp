module;

module mc_moves;

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
import transition_matrix;


void MC_Moves::performRandomMove(System& selectedSystem, System& selectedSecondSystem, size_t selectedComponent, size_t &fractionalMoleculeSystem)
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Translation);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Translation);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Rotation);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Rotation);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Reinsertion);
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityIdentityChangeMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CBMC)
  {
    size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];
    if (RandomNumber::Uniform() < 0.5)
    {
      const auto [energyDifference, Pacc] = insertionMove(selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.update(Pacc, oldN, TransitionMatrix::MoveType::Insertion);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(selectedComponent);

      const auto [energyDifference, Pacc] = deletionMove(selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.update(Pacc, oldN, TransitionMatrix::MoveType::Deletion);
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

    std::optional<RunningEnergy> energyDifference = swapMove_CFCMC(selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
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
    std::optional<std::pair<RunningEnergy,RunningEnergy>> energy = GibbsVolumeMove(selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CBMC)
  {
    if(RandomNumber::Uniform() < 0.5)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(selectedSystem, selectedSecondSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(selectedSecondSystem, selectedSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(selectedSystem, selectedSecondSystem, fractionalMoleculeSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else if(selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(selectedSecondSystem, selectedSystem, fractionalMoleculeSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC_CBMC(selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else if(selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC_CBMC(selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
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

void MC_Moves::performRandomMoveProduction(System& selectedSystem, System& selectedSecondSystem, size_t selectedComponent, size_t &fractionalMoleculeSystem, size_t currentBlock)
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Translation);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_TranslationMove += (t2 - t1);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Translation);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_RandomTranslationMove += (t2 - t1);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Rotation);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_RotationMove += (t2 - t1);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Rotation);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_RandomRotationMove += (t2 - t1);
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
    selectedSystem.mc_moves_probabilities.cpuTime_VolumeMove += (t2 - t1);
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
      selectedSystem.tmmc.update(1.0, selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent], TransitionMatrix::MoveType::Reinsertion);
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_ReinsertionMove_CBMC += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityIdentityChangeMove_CBMC)
  {
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CBMC)
  {
    size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];
    if (RandomNumber::Uniform() < 0.5)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] = insertionMove(selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.update(Pacc, oldN, TransitionMatrix::MoveType::Insertion);

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CBMC += (t2 - t1);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] = deletionMove(selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.update(Pacc, oldN, TransitionMatrix::MoveType::Deletion);
      
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapDeletionMove_CBMC += (t2 - t1);
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilitySwapMove_CFCMC)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::optional<RunningEnergy> energyDifference = swapMove_CFCMC(selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapMove_CFCMC += (t2 - t1);
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
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapMove_CFCMC_CBMC += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsVolumeMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsVolumeMove(selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.mc_moves_probabilities.cpuTime_GibbsVolumeMove += (t2 - t1);
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CBMC)
  {
    if (RandomNumber::Uniform() < 0.5)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(selectedSystem, selectedSecondSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapMove_CBMC += (t2 - t1);
    }
    else
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CBMC(selectedSystem, selectedSecondSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapMove_CBMC += (t2 - t1);
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC += (t2 - t1);
    }
    else if (selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC(selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC += (t2 - t1);
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC_CBMC(selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC += (t2 - t1);
    }
    else if (selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = GibbsSwapMove_CFCMC_CBMC(selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC += (t2 - t1);
    }
  }
  else if (randomNumber < selectedSystem.components[selectedComponent].mc_moves_probabilities.accumulatedProbabilityWidomMove)
  {
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::optional<double> RosenbluthWeight = WidomMove(selectedSystem, selectedComponent);

    selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, RosenbluthWeight.value_or(0.0), selectedSystem.weight());

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CBMC += (t2 - t1);
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
    selectedSystem.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_CBMC += (t2 - t1);
  }
}

