module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <print>
#include <source_location>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves;

#ifndef USE_LEGACY_HEADERS
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
import <exception>;
import <source_location>;
import <complex>;
import <print>;
#endif

import archive;
import double3;
import double3x3;
import simd_quatd;
import randomnumbers;

import component;
import atom;
import molecule;
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
import mc_moves_translation;
import mc_moves_random_translation;
import mc_moves_rotation;
import mc_moves_random_rotation;
import mc_moves_reinsertion;
import mc_moves_insertion;
import mc_moves_deletion;
import mc_moves_insertion_cbmc;
import mc_moves_deletion_cbmc;
import mc_moves_swap_cfcmc;
import mc_moves_swap_cfcmc_cbmc;
import mc_moves_gibbs_swap_cbmc;
import mc_moves_volume;
import mc_moves_gibbs_volume;
import mc_moves_identity_change;
import mc_moves_swap_cfcmc;
import mc_moves_swap_cfcmc_cbmc;
import mc_moves_gibbs_swap_cbmc;
import mc_moves_gibbs_swap_cfcmc;
import mc_moves_reaction;
import mc_moves_reaction_cfcmc_cbmc;
import mc_moves_widom;
import mc_moves_parallel_tempering_swap;
import mc_moves_hybridmc;

void MC_Moves::performRandomMove(RandomNumber &random, System &selectedSystem, System &selectedSecondSystem,
                                 size_t selectedComponent, size_t &fractionalMoleculeSystem)
{
  double randomNumber = random.uniform();

  MCMoveProbabilitiesParticles &mc_moves_probabilities =
      selectedSystem.components[selectedComponent].mc_moves_probabilities;

  size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  if (randomNumber < mc_moves_probabilities.accumulatedTranslationProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference =
          MC_Moves::translationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                    selectedSystem.components, molecule, molecule_atoms);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedRandomTranslationProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::randomTranslationMove(
          random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedRotationProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::rotationMove(
          random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedRandomRotationProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::randomRotationMove(
          random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedVolumeChangeProbability)
  {
    std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedReinsertionCBMCProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::reinsertionMove(
          random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }

      selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedIdentityChangeCBMCProbability)
  {
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapProbability)
  {
    if (random.uniform() < 0.5)
    {
      const auto [energyDifference, Pacc] = MC_Moves::insertionMove(random, selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

      const auto [energyDifference, Pacc] =
          MC_Moves::deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapCBMCProbability)
  {
    if (random.uniform() < 0.5)
    {
      const auto [energyDifference, Pacc] = MC_Moves::insertionMoveCBMC(random, selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
    }
    else
    {
      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

      const auto [energyDifference, Pacc] =
          MC_Moves::deletionMoveCBMC(random, selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapCFCMCProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    selectedSystem.tmmc.updateMatrix(Pacc, oldN);
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapCBCFCMCProbability)
  {
    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, selectedMolecule);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
    selectedSystem.tmmc.updateMatrix(Pacc, oldN);
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedGibbsVolumeChangeProbability)
  {
    std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
        MC_Moves::GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedGibbsSwapCBMCProbability)
  {
    if (random.uniform() < 0.5)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedGibbsSwapCFCMCProbability)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
          random, selectedSystem, selectedSecondSystem, fractionalMoleculeSystem, selectedComponent);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
    }
    else if (selectedSecondSystem.containsTheFractionalMolecule)
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
          random, selectedSecondSystem, selectedSystem, fractionalMoleculeSystem, selectedComponent);
      if (energy)
      {
        selectedSecondSystem.runningEnergies += energy.value().first;
        selectedSystem.runningEnergies += energy.value().second;
      }
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedWidomProbability)
  {
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedWidomCFCMCProbability)
  {
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedWidomCBCFCMCProbability)
  {
    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
    if (energyDifference)
    {
      selectedSystem.runningEnergies += energyDifference.value();
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedParallelTemperingProbability)
  {
    std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
        MC_Moves::ParallelTemperingSwap(random, selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedHybridMCProbability)
  {
    std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
  }
}

void MC_Moves::performRandomMoveProduction(RandomNumber &random, System &selectedSystem, System &selectedSecondSystem,
                                           size_t selectedComponent, size_t &fractionalMoleculeSystem,
                                           size_t currentBlock)
{
  double randomNumber = random.uniform();

  MCMoveProbabilitiesParticles &mc_moves_probabilities =
      selectedSystem.components[selectedComponent].mc_moves_probabilities;
  size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  if (randomNumber < mc_moves_probabilities.accumulatedTranslationProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.translationMove.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference =
          MC_Moves::translationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                    selectedSystem.components, molecule, molecule_atoms);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedRandomTranslationProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.randomTranslationMove.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::randomTranslationMove(
          random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedRotationProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.rotationMove.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::rotationMove(
          random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedRandomRotationProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.randomRotationMove.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::randomRotationMove(
          random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedVolumeChangeProbability)
  {
    selectedSystem.mc_moves_statistics.volumeMove.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.mc_moves_cputime.volumeMove += (t2 - t1);
    selectedSystem.mc_moves_count.volumeMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedReinsertionCBMCProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
    {
      std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
      size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
      Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
      std::optional<RunningEnergy> energyDifference = MC_Moves::reinsertionMove(
          random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

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
  else if (randomNumber < mc_moves_probabilities.accumulatedIdentityChangeCBMCProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.identityChangeMove_CBMC.allCounts += 1uz;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapProbability)
  {
    if (random.uniform() < 0.5)
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.swapInsertionMove.allCounts += 1uz;

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] = MC_Moves::insertionMove(random, selectedSystem, selectedComponent);

      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.swapInsertionMove += (t2 - t1);
      selectedSystem.mc_moves_cputime.swapInsertionMove += (t2 - t1);

      selectedSystem.components[selectedComponent].mc_moves_count.swapInsertionMove++;
      selectedSystem.mc_moves_count.swapInsertionMove++;
    }
    else
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.swapDeletionMove.allCounts += 1uz;

      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] =
          MC_Moves::deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

      if (energyDifference)
      {
        selectedSystem.runningEnergies -= energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.swapDeletionMove += (t2 - t1);
      selectedSystem.mc_moves_cputime.swapDeletionMove += (t2 - t1);

      selectedSystem.components[selectedComponent].mc_moves_count.swapDeletionMove++;
      selectedSystem.mc_moves_count.swapDeletionMove++;
    }
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapCBMCProbability)
  {
    if (random.uniform() < 0.5)
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.allCounts += 1uz;

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] = MC_Moves::insertionMoveCBMC(random, selectedSystem, selectedComponent);

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
      selectedSystem.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.allCounts += 1uz;

      size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] =
          MC_Moves::deletionMoveCBMC(random, selectedSystem, selectedComponent, selectedMolecule);

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
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapCFCMCProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.swapMove_CFCMC.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, selectedMolecule);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedSwapCBCFCMCProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.swapMove_CFCMC_CBMC.allCounts += 1uz;

    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, selectedMolecule);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedGibbsVolumeChangeProbability)
  {
    selectedSystem.mc_moves_statistics.GibbsVolumeMove.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
        MC_Moves::GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.mc_moves_cputime.GibbsVolumeMove += (t2 - t1);
    selectedSystem.mc_moves_count.GibbsVolumeMove++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedGibbsSwapCBMCProbability)
  {
    if (random.uniform() < 0.5)
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CBMC.allCounts += 1uz;

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
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
      selectedSecondSystem.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CBMC.allCounts += 1uz;

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedGibbsSwapCFCMCProbability)
  {
    if (selectedSystem.containsTheFractionalMolecule)
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.allCounts += 1uz;

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
          random, selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
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
      selectedSecondSystem.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.allCounts += 1uz;
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
          random, selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedWidomProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.WidomMove_CBMC.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    std::pair<double, double> value = MC_Moves::WidomMove(random, selectedSystem, selectedComponent);

    selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(
        currentBlock, value.first, 2.0 * value.second, selectedSystem.weight());

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCBMC += (t2 - t1);
    selectedSystem.mc_moves_cputime.WidomMoveCBMC += (t2 - t1);

    selectedSystem.components[selectedComponent].mc_moves_count.WidomMoveCBMC++;
    selectedSystem.mc_moves_count.WidomMoveCBMC++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedWidomCFCMCProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.WidomMove_CFCMC.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, 0, true, true);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedWidomCBCFCMCProbability)
  {
    selectedSystem.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC_CBMC.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    const auto [energyDifference, Pacc] =
        MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
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
  else if (randomNumber < mc_moves_probabilities.accumulatedParallelTemperingProbability)
  {
    selectedSystem.mc_moves_statistics.ParallelTemperingSwap.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
        MC_Moves::ParallelTemperingSwap(random, selectedSystem, selectedSecondSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value().first;
      selectedSecondSystem.runningEnergies = energy.value().second;
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.mc_moves_cputime.ParallelTemperingSwap += (t2 - t1);
    selectedSystem.mc_moves_count.ParallelTemperingSwap++;
  }
  else if (randomNumber < mc_moves_probabilities.accumulatedHybridMCProbability)
  {
    selectedSystem.mc_moves_statistics.hybridMC.allCounts += 1uz;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
    if (energy)
    {
      selectedSystem.runningEnergies = energy.value();
    }
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    selectedSystem.mc_moves_cputime.hybridMC += (t2 - t1);
    selectedSystem.mc_moves_count.hybridMC++;
  }
}
